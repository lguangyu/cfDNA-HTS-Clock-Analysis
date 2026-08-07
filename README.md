# cfDNA HTS Clock Analysis

This repository contains all code to generate the result figures of paper "A robust and integrated framework for cross-platform adaptation of epigenetic clocks in cell-free DNA sequencing"

## The DF-IM-TL pipeline

### 1. Train a PCA model

The examples can be found under `pipeline/`. To run this example, a trained PCA model is required. I'd prefer a PCA model that is trained with a sufficiently large dataset such as GSE55763, though it may take sometime to do so. The user may need to prepare a beta matrix file, and pre-determine the number of components. To obtain such a model, call run the following command:

```bash
# in this case, the gse55763 dataset contains 2711 samples, and the number of
# compoents is chosen as 2500
python3 pipeline_pca.py \
	-b gse55763.beta.tsv \
	-o gse55763.pca.pkl.gz \
	-n 2500
```

Users may package their own PCA model instead of training. The data structure of the PCA model bundle is:

```python3
obj = {
	"cpgs": List[str],
	"pca_obj": sklearn.decomposition.PCA,
	"n_components": int,
}
```

Where:
"cpgs" is a list of CpG reference IDs that are compatible with the beta/depth matrices.
"pca_obj" is a fitted instance of sklearn.decomposition.PCA class
"n_components" is a redundant record of n_components used in the "pca_obj" (optional).

Note:  the "cpgs" must be in element-wise correspondance to the input dimensions of the "pca_obj", both in order and in size.

### 2. Fit the pipeline

The next step is to fit the TL step of the adaptation pipeline, despite that the DF and IM step doesn't require training. The input to the fitting script requires:

* the metadata table of the input dataset (sample x features)
* the beta matrix of the input dataset (cpgs x samples)
* the depth matrix of the input dataset (cpgs x samples)
* a PCA model bundle (as one obtained in the previous step)

With these inputs ready, the fitting command can be executed by:

```bash
python3 pipeline_train.py \
	-m input.metadata.tsv \
	-b input.beta.tsv \
	-d input.depth.tsv \
	--pca-model gse_55763.pca.pkl.gz \
	--depth-filter-thres 10 \
	--imputation-depth-thres 20 \
	--imputation-method knn \
	--clock-names horvath2013,hannum,han,zhangblup \
	-o adapt_test.pkl.gz
```

Check for a complete list of command-line options by `python3 pipeline_train.py -h`. When it finishes, the fitted pipeline data bundle is saved in file `piepline_fitted.pkl.gz`.


### 3. Apply the pipeline

To apply the fitted pipeline on a target dataset, the following data is required:

* the fitted pipeline file (pkl.gz)
* the beta matrix of the input dataset (cpgs x samples)
* the depth matrix of the input dataset (cpgs x samples)
* a PCA model bundle (should be in the same path as when the pipeline is trained)

Then call with command:

```bash
python3 pipeline_apply.py \
	-b target.beta.tsv \
	-d target.depth.tsv \
	-p adapt_test.pkl.gz \
	-o adapt_test.pred.tsv
```

The generated file contains predictions of adapted models.


## The ClinicalAge

ClinicalAge is an in-house trained model to predict age from clinical measurements without using chronological age as an input. The related script can be found under `clinical_age/`.

A tiny run example of here:

```bash
```