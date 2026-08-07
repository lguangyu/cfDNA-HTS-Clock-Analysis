#!/usr/bin/env python3

import argparse
import pandas
import pdb


def get_args():
	ap = argparse.ArgumentParser()
	ap.add_argument("-i", "--input", required=True,
		metavar="tsv",
		help="input clinical data table (required)",
	)
	ap.add_argument("-o", "--output", required=True,
		metavar="tsv",
		help="output clinical age predictions (required)",
	)

	# parse and refine args
	args = ap.parse_args()

	return args


def main():
	args = get_args()

	in_df = pandas.read_csv(args.input, sep="\t", index_col=0)

	coef_df = pandas.read_csv("result/coef.txt", sep="\t")

	intercept_row = coef_df.iloc[0]
	intercept = intercept_row["COEF"]

	coef_df = coef_df.iloc[1:, :]
	features = coef_df["FEAT_CODE"]
	scaler_mean = coef_df["SCALER_MEAN"]
	scaler_std = coef_df["SCALER_STD"]
	coef = coef_df["COEF"]

	curated_df = pandas.DataFrame(index=in_df.index, columns=features)
	for feat in features:
		if feat not in in_df.columns:
			raise ValueError(f"Feature {feat} not found in input data")
		curated_df[feat] = in_df[feat]
	curated_df.fillna(0, inplace=True)

	# predict
	pred_values = (curated_df.values - scaler_mean.values.reshape(1, -1)
	    ) / scaler_std.values.reshape(1, -1)
	pred_values = pred_values @ coef.values.reshape(-1, 1) + intercept
	pred_df = pandas.DataFrame(pred_values, index=in_df.index,
		columns=["clinical_age"])

	pred_df.to_csv(args.output, sep="\t", header=True)
	return


if __name__ == "__main__":
	main()
