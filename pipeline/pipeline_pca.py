#!/usr/bin/env python3

import argparse
import gzip
import pickle

import pandas
import numpy
import sklearn
import sklearn.decomposition


class NonNegInt(int):
	def __new__(cls, value):
		value = int(value)
		if value < 0:
			raise ValueError(f"Value must be 0 or a positive integer, got {value}")
		return int.__new__(cls, value)


def get_args():
	ap = argparse.ArgumentParser()
	ap.add_argument("-b", "--beta", required=True,
		metavar="tsv",
		help="methylation beta matrix (required)",
	)
	ap.add_argument("-l", "--cpgs-list", type=str,
		metavar="txt",
		help="list of CpGs to subset the beta matrix [use all]",
	)
	ap.add_argument("-o", "--output", required=True,
		metavar="pkl.gz",
		help="output pca results (required)",
	)

	ap.add_argument("--n-components", type=NonNegInt, required=True,
		metavar="int",
		help="number of PCA components to retain (required)"
	)

	# parse and refine args
	args = ap.parse_args()

	return args


def main():
	args = get_args()

	beta = pandas.read_csv(args.beta, sep="\t", index_col=0).T
	if args.cpgs_list is not None:
		with open(args.cpgs_list) as fp:
			cpgs = [v.strip() for v in fp]
	else:
		cpgs = beta.columns.tolist()
	beta = beta.reindex(columns=cpgs)
	beta = numpy.nan_to_num(beta.values, nan=0.0)

	pca_obj = sklearn.decomposition.PCA(n_components=args.n_components)
	pca_obj.fit(beta)

	out_obj = {
		"cpgs": cpgs,
		"pca_obj": pca_obj,
		"n_components": args.n_components,
	}
	with gzip.open(args.output, "wb") as fp:
		pickle.dump(out_obj, fp)

	return


if __name__ == "__main__":
	main()
