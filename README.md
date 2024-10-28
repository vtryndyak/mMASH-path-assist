# Predicting MASH score using machine learning algorithm

This script is intended to use pre-trained machine learning model to predict MASH score using gene expression data. To make the prediction, please follow these steps:

1. Create conda environment.

At this step, we will ensure that all dependencies are present. As the model was trained with h2o, it requires a specific version. The r-h2o used in training was version 3.44.0.3, and it comes with this conda environment. To create the environment:
    
```bash
conda env create -f predict_mash_score.R
conda activate h2o-env
```

**_NOTE:_** to learn more about conda environments (specifying path, change name, etc), please follow official [conda documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)

2. Prepare the input data.

The input gene expression data is expected to match certain criteria. Please follow these guidelines to prepare the data for MASH score prediction:
    
- The format is expected to be a comma-separated file (`.csv`).
- The first column is expected to be a gene indentifier. The identifier column is handled by position and therefore the name of the column is not taken into consideration and can be any value.
- Each subsequent column is expected to represent an individual sample. The name of the column will be used as sample id.
- No other columns besides gene identifier and sample data are allowed.
- The expression of all genes should be provided as there is an internal step of per-sample scaling.
- Please refer to the supplied `test.csv` for an example of the properly formatted data.

3. Run the script to make the prediction.

The script accepts the following positional unnamed arguments:
    
- <path/to/ge_data.csv> Path to the input file with gene expression data formatted as described in step 2.
- <model> String indicating which model to use. Accepted values are `lm` for linear model, `nn` for neural network, `rf` for random forest, and `gbm` for gradient boosting machine. The `lm` is a recommended value.
- <path/to/output> Path to the file where predictions will be saved.

To make predictions, run the provided R script in the terminal vindow. For example:
    
```bash
Rscript predict_mash_score.R test.csv lm predictions.csv
```
To use random forest model, simply modify the value of model argument:

```bash
Rscript predict_mash_score.R test.csv rf predictions.csv
```

The output will contain the following columns:

- predict: The predicted MASH score
- sample_id: The sample identifier inferred from the column name in the input data
- trained_on: A string with `MF` value indicating model was originally developed on data from both males and females. Used for development purpose and can be ignored by end user.
- gene_set: A panel of genes used for the model development. Used for development purpose and can be ignored by end user.
- Model: a string indicating which model was used for the prediction.
