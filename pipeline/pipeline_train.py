#!/usr/bin/env python3

import argparse
import gzip
import pickle
import pdb
import warnings

import numpy
import pandas
import pyaging
import sklearn
import sklearn.decomposition
import sklearn.linear_model
import sklearn.model_selection


class NonNegInt(int):
	def __new__(cls, value):
		value = int(value)
		if value < 0:
			raise ValueError(f"Value must be 0 or a positive integer, got {value}")
		return int.__new__(cls, value)


class CommaSeparatedList(argparse.Action):
	def __call__(self, parser, namespace, values, option_string=None):
		setattr(namespace, self.dest, values.split(","))
		return


def get_args():
	ap = argparse.ArgumentParser()
	ap.add_argument("-m", "--meta", required=True,
		metavar="tsv",
		help="metadata table for samples (required)",
	)
	ap.add_argument("-b", "--beta", required=True,
		metavar="tsv",
		help="methylation beta matrix (required)",
	)
	ap.add_argument("-d", "--depth", required=True,
		metavar="tsv",
		help="depth matrix (required)",
	)
	ap.add_argument("-o", "--output", required=True,
		metavar="pkl.gz",
		help="output adaptation results (required)",
	)

	# below are adaptation arguments
	ap.add_argument("--aggregate-epicv2-probes", action="store_true",
		help="aggregate EPICv2 probes (default: %(default)s)"
	)
	ap.add_argument("--pca-model", type=str, required=True,
		metavar="pkl.gz",
		help="path to PCA model file (required)"
	)
	ap.add_argument("--depth-filter-thres", type=NonNegInt, default=10,
		metavar="int",
		help="minimum depth threshold for filtering (default: %(default)s)",
	)
	ap.add_argument("--imputation-depth-thres", type=NonNegInt, default=20,
		metavar="int",
		help="maximum depth threshold for imputation (default: %(default)s)",
	)
	ap.add_argument("--imputation-method", type=str, default="knn",
		choices=["knn", "mean", "median", "constant"],
		help="imputation method for missing values (default: %(default)s)"
	)
	ap.add_argument("--clock-names", action=CommaSeparatedList, required=True,
		metavar="clock[,clock[,...]]",
		help="comma-separated list of clock names, e.g. horvath2013,hannum (required); check pyaging for available clocks",
	)
	ap.add_argument("--pyaging-data", type=str, default="pyaging_data",
		metavar="dir",
		help="path to pyaging data directory (default: %(default)s)",
	)
	ap.add_argument("--cv-fold", type=NonNegInt, default=5,
		metavar="int",
		help="number of folds for K-Fold cross-validation (default: %(default)s)",
	)
	ap.add_argument("--cv-seed", type=int, default=0,
		metavar="int",
		help="random seed for K-Fold cross-validation (default: %(default)s)",
	)

	# parse and refine args
	args = ap.parse_args()

	return args


def main():
	args = get_args()

	meta = pandas.read_csv(args.meta, sep="\t", index_col=0)
	beta = pandas.read_csv(args.beta, sep="\t", index_col=0).T
	depth = pandas.read_csv(args.depth, sep="\t", index_col=0).T

	if not meta.index.is_unique:
		raise ValueError("metadata has non-unique sample IDs")
	if not meta.columns.is_unique:
		raise ValueError("metadata matrix has non-unique column IDs")

	meta_samples = set(meta.index)
	beta_samples = set(beta.index)
	depth_samples = set(depth.index)
	if not (beta_samples == meta_samples == depth_samples):
		warnings.warn(
			"sample ids in beta/depth matrix do not match metadata, non-matching samples will be dropped"
		)
	shared_samples = meta.index.intersection(beta.index).intersection(depth.index)
	meta = meta.loc[shared_samples]
	beta = beta.loc[shared_samples]
	depth = depth.loc[shared_samples]
	if "age" not in meta.columns:
		raise ValueError("metadata does not contain 'age' column")
	meta_cols = meta.columns.tolist()

	beta_cpgs = set(beta.columns)
	depth_cpgs = set(depth.columns)
	if beta_cpgs != depth_cpgs:
		warnings.warn(
			"CpG ids in beta/depth matrix do not match, non-matching CpGs will be dropped"
		)
	shared_cpgs = beta.columns.intersection(depth.columns)
	beta = beta.loc[:, shared_cpgs]
	depth = depth.loc[:, shared_cpgs]

	df_filt = depth < args.depth_filter_thres
	depth[df_filt] = numpy.nan

	im_filt = ((beta == 0) | (beta == 1)) & (depth <= args.imputation_depth_thres)
	depth[im_filt] = numpy.nan

	if args.aggregate_epicv2_probes:
		beta = pyaging.pp.epicv2_probe_aggregation(beta, verbose=False)

	concat_df = pandas.concat([meta, beta], axis=1)
	adata = pyaging.pp.df_to_adata(concat_df,
		metadata_cols=meta_cols,
		imputer_strategy=args.imputation_method,
		verbose=False,
	)
	pyaging.pred.predict_age(adata,
		clock_names=args.clock_names,
		dir=args.pyaging_data,
		verbose=False,
	)
	base_pred = adata.obs

	with gzip.open(args.pca_model, "rb") as fp:
		pca_bundle = pickle.load(fp)
	pca_cpgs = pca_bundle["cpgs"]
	pca_obj = pca_bundle["pca_obj"]

	beta_pca = beta.reindex(columns=pca_cpgs, fill_value=0).values
	numpy.nan_to_num(beta_pca, copy=False)
	beta_pca = pca_obj.transform(beta_pca)

	adapt_pred = meta.copy()
	adapt_models = dict()
	for clock in args.clock_names:
		clock_pred = numpy.full(len(adapt_pred), numpy.nan)
		# concatenate base predictions with PCA features
		clock_x = numpy.hstack([base_pred[clock].values.reshape(-1, 1), beta_pca])
		self_cv_pred = numpy.full(len(adapt_pred), numpy.nan)
		model_list = list()
		for train_idx, test_idx in sklearn.model_selection.KFold(
			n_splits=args.cv_fold, shuffle=True, random_state=args.cv_seed
		).split(beta_pca):
			train_x = clock_x[train_idx]
			train_y = meta["age"].values[train_idx]
			test_x = clock_x[test_idx]
			_model = sklearn.linear_model.ElasticNetCV(
				alphas=numpy.logspace(-4, 4, 100),
				l1_ratio=numpy.linspace(0, 1, 11),
				cv=5,
				n_jobs=-1,
				random_state=args.cv_seed,
			)
			_model.fit(train_x, train_y)
			model_list.append(_model)
			self_cv_pred[test_idx] = _model.predict(test_x)
		adapt_pred[clock] = self_cv_pred
		_final_alpha = numpy.median([m.alpha_ for m in model_list])
		_final_l1_ratio = numpy.median([m.l1_ratio_ for m in model_list])
		final_model = sklearn.linear_model.ElasticNet(
			alpha=_final_alpha,
			l1_ratio=_final_l1_ratio,
			random_state=args.cv_seed,
		)
		final_model.fit(clock_x, meta["age"].values)
		adapt_models[clock] = final_model

	out_obj = {
		"adapt_pred": adapt_pred,
		"adapt_models": adapt_models,
		"pca_model": args.pca_model,
		"depth_filter_thres": args.depth_filter_thres,
		"imputation_depth_thres": args.imputation_depth_thres,
		"imputation_method": args.imputation_method,
		"clock_names": args.clock_names,
	}
	with gzip.open(args.output, "wb") as fp:
		pickle.dump(out_obj, fp)

	return


if __name__ == "__main__":
	main()
