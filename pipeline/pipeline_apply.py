#!/usr/bin/env python3

import argparse
import gzip
import pickle
import warnings

import numpy
import pandas
import pyaging


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
	ap.add_argument("-b", "--beta", required=True,
		metavar="tsv",
		help="methylation beta matrix (required)",
	)
	ap.add_argument("-d", "--depth", required=True,
		metavar="tsv",
		help="depth matrix (required)",
	)
	ap.add_argument("-p", "--trained-pipeline", type=str, required=True,
		metavar="pkl.gz",
		help="path to trained pipeline file (required)"
	)
	ap.add_argument("-o", "--output", required=True,
		metavar="tsv",
		help="output adaptation prediction results (required)",
	)

	# below are adaptation arguments
	ap.add_argument("--aggregate-epicv2-probes", action="store_true",
		help="aggregate EPICv2 probes (default: %(default)s)"
	)
	ap.add_argument("--pyaging-data", type=str, default="pyaging_data",
		metavar="dir",
		help="path to pyaging data directory (default: %(default)s)",
	)

	# parse and refine args
	args = ap.parse_args()

	return args


def main():
	args = get_args()

	beta = pandas.read_csv(args.beta, sep="\t", index_col=0).T
	depth = pandas.read_csv(args.depth, sep="\t", index_col=0).T

	beta_samples = set(beta.index)
	depth_samples = set(depth.index)
	if not (beta_samples == depth_samples):
		warnings.warn(
			"sample ids in beta/depth matrix do not match metadata, non-matching samples will be dropped"
		)
	shared_samples = beta.index.intersection(depth.index)
	beta = beta.loc[shared_samples]
	depth = depth.loc[shared_samples]

	with gzip.open(args.trained_pipeline, "rb") as fp:
		pipeline_bundle = pickle.load(fp)
	adapt_models = pipeline_bundle["adapt_models"]

	with gzip.open(pipeline_bundle["pca_model"], "rb") as fp:
		pca_bundle = pickle.load(fp)
	cpgs = pca_bundle["cpgs"]
	pca_obj = pca_bundle["pca_obj"]

	beta = beta.reindex(columns=cpgs)
	depth = depth.reindex(columns=cpgs, fill_value=0)

	df_filt = depth < pipeline_bundle["depth_filter_thres"]
	depth[df_filt] = numpy.nan

	im_filt = ((beta == 0) | (beta == 1)) & (
		depth <= pipeline_bundle["imputation_depth_thres"])
	depth[im_filt] = numpy.nan

	if args.aggregate_epicv2_probes:
		beta = pyaging.pp.epicv2_probe_aggregation(beta, verbose=False)

	adata = pyaging.pp.df_to_adata(beta,
		metadata_cols=[],
		imputer_strategy=pipeline_bundle["imputation_method"],
		verbose=False,
	)
	pyaging.pred.predict_age(adata,
		clock_names=pipeline_bundle["clock_names"],
		dir=args.pyaging_data,
		verbose=False,
	)
	base_pred = adata.obs

	beta_pca = beta.values
	numpy.nan_to_num(beta_pca, copy=False)
	beta_pca = pca_obj.transform(beta_pca)

	adapt_pred = pandas.DataFrame(index=beta.index)
	for clock in pipeline_bundle["clock_names"]:
		# concatenate base predictions with PCA features
		clock_x = numpy.hstack([base_pred[clock].values.reshape(-1, 1), beta_pca])
		adapt_pred[clock] = adapt_models[clock].predict(clock_x)

	adapt_pred.to_csv(args.output, sep="\t", index=True, header=True)
	return


if __name__ == "__main__":
	main()
