### Predict MESH score using GE data ### 
# Uses pre-trained h2o models to predict the MASH score. 

#!/usr/bin/env Rscript
#
# Usage: 
#   Rscript predict_mash_score.R <path/to/ge_data.csv> <model> <path/to/output>
#
# Notes: 
#   This script expects to have gene expression matrix at the input path, following 
#   the structure where first column represent gene identifier and all subsequent
#   columns represent individual samples. No other columns are allowed in the
#   input data. The input data is expected to be in the csv format.
#
#   The gene expression matrix is supposed to contain raw counts for all genes
#   without prior sub-setting in order to perform proper scaling. If the data is
#   pre-filtered to a particular gene set, the performance of the model is not
#   guaranteed.
#

# Load packages -----------------------------------------------------------
message("Loading packages...")
suppressPackageStartupMessages({
    library(h2o)
    library(dplyr)
    library(tibble)
    library(readr)
})

# Determine arguments -----------------------------------------------------

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE) %>% as.list()
arg_names <- c("ge_path", "model", "output_path")
args <- setNames(args, arg_names[1:length(args)])

# Check arguments

if(! endsWith(args$ge_path, '.csv')){
    stop(
        "The input data is not in csv format (does not end in .csv)."
    )
}

args$model <- tolower(args$model)
available_models <- c(
    "rf", "lm", "nn", "gbm"
)

if(! args$model %in% available_models){
    message("The requested model is not awailable!")
    stop(
        "Please choose one of lm, nn, rf, or gbm."
    )
}

# The set of 74 predictor genes
genes74 <- readRDS(file = "predictor_genes.rds")

# Import raw counts of GE data
message("Reading input data...")
data <- suppressMessages(read_csv(args$ge_path))

colnames(data)[1] <- "gene_symbol"

# Perform scaling
message("Performing per-sample scaling...")
data <- data %>%
    mutate(Total = rowSums(select_if(., is.numeric), na.rm = TRUE)) %>%
    group_by(gene_symbol) %>%
    slice_head() %>%
    ungroup %>%
    filter(!is.na(gene_symbol)) %>%
    column_to_rownames("gene_symbol")

data[data < 0] <- 0

# Scale per sample and only subset to 74 predictors
data_scaled <- sweep(data,2,colSums(data),`/`)
data_scaled <- t(data_scaled[genes74,-ncol(data_scaled)])

# Make predictions
message("Making predictions...")

model_path <- paste0(
    "combined_genes74_",
    args$model
)

suppressWarnings(h2o.init())

model <- h2o.loadModel(
    model_path
)


predictions <- predict(
    model,
    newdata = as.h2o(data_scaled)
) %>%
    round() %>%
    as.data.frame %>%
    mutate(
        sample_id = rownames(data_scaled),
        trained_on = "MF",
        gene_set = "genes74",
        Model = args$model
    )

# Save predicitons
message("Saving predictions...")

predictions %>%
    arrange(trained_on, sample_id) %>%
    mutate(predict = ifelse(predict<0, 0, predict)) %>%
    write_csv(args$output_path)

message("Done!")
