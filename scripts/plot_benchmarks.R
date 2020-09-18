#!/usr/bin/env Rscript
library(tidyverse)

# --- Functions to get file paths ---

#return the data directory appended to path
get_data_dir <- function(path){
  return(paste0("../data/", path))
}

#return file name from stems where subpath
#is a dir in the data directory and suffix
#is the part after the stem.
#[../data/subpath/stem1suffix, ../data/subpath/stem2suffix, ...]
get_file <- function(stem, subpath, suffix){
  dir <- paste0(subpath, paste0('/', stem))
  return(paste0(get_data_dir(dir), suffix))
}

#return ../data/subpath/stem.recal.after.csv,
# the default is ../data/stem.recal.after.csv
get_gatk_csv_file <- function(stem, subpath = '.'){
  return(get_file(stem, subpath, ".recal.after.csv"))
}

# --- Utility Functions ---

#turn a probability to a phred-scaled quality score
p_to_q <- function(p){
  return(floor(-10*log10(p)))
}

# --- Dataframe Manipulation ---

#Gets only QualityScore and After rows in the dataframe
#and changes CovariateValue title to PredictedQuality
get_quality_rows <- function(df){
  df %>%
    filter(CovariateName == 'QualityScore' & Recalibration == 'After') %>%
    rename(PredictedQuality = CovariateValue)
}

#Returns a df with covariatevalue and reportedquality columns
average_over_rgs_fnr <- function(df){
  df %>%
    get_quality_rows() %>%
    mutate_at(vars(FalseNegativeRate, PredictedQuality), as.numeric) %>%
    group_by(FalseNegativeRate, PredictedQuality) %>%
    summarize(ActualQuality = p_to_q(sum(Errors)/sum(Observations)))
}

average_over_rgs_fpr <- function(df){
  df %>%
    get_quality_rows() %>%
    mutate_at(vars(FalsePositiveRate, PredictedQuality), as.numeric) %>%
    group_by(FalsePositiveRate, PredictedQuality) %>%
    summarize(ActualQuality = p_to_q(sum(Errors)/sum(Observations)))
}

average_over_rgs_fnrfpr <- function(df){
  df %>%
    get_quality_rows() %>%
    mutate_at(vars(FalsePositiveRate, FalseNegativeRate, PredictedQuality),
              as.numeric) %>%
    group_by(FalseNegativeRate, FalsePositiveRate, PredictedQuality) %>%
    summarize(ActualQuality = p_to_q(sum(Errors)/sum(Observations)))
}

# --- Importing CSVs ---

import_fnr_csvs <- function(fnr){
  csv_files <- get_gatk_csv_file(fnr)
  dfs <- map(csv_files, read_csv)
  names(dfs) <- fnr
  df <- bind_rows(dfs, .id = 'FalseNegativeRate')
  return(df)
}

import_fpr_csvs <- function(fpr){
  csv_files <- get_gatk_csv_file(fpr, subpath = 'fpr')
  dfs <- map(csv_files, read_csv)
  names(dfs) <- fpr
  df <- bind_rows(dfs, .id = 'FalsePositiveRate')
  return(df)
}

import_fnrfpr_csvs <- function(fnr, fpr){
  stems <- unlist(map(fnr, paste0, '_', fpr))
  csv_files <- get_gatk_csv_file(stems, subpath = 'fnrfpr')
  dfs <- map(csv_files, read_csv)
  names(dfs) <- stems
  df <- bind_rows(dfs, .id = 'FNR_FPR') %>%
    separate(FNR_FPR, into = c("FalseNegativeRate","FalsePositiveRate"),
             sep = "_", convert = TRUE)
}

# --- Plotting Data ---

plot_fnr_csvs <- function(gatk_df){
  ggplot(gatk_df, aes(PredictedQuality, ActualQuality)) +
    geom_abline(slope = 1, intercept = 0) +
    geom_point(aes(color = factor(FalseNegativeRate)), size = 2) +
    geom_line(aes(color = factor(FalseNegativeRate)), size = 1) +
    scale_color_brewer('Known Sites\nFalse Negative\nRate', palette = 'PRGn') +
    scale_x_continuous("Predicted Quality") +
    scale_y_continuous("Actual Quality") +
    ggtitle('False Negative Variants Cause Underconfidence') +
    coord_fixed(ratio = 1) +
    theme_minimal(base_size = 18) + 
    theme(plot.margin = margin(0,0,0,0))
}

plot_fpr_csvs <- function(gatk_df){
  ggplot(gatk_df, aes(PredictedQuality, ActualQuality)) +
    geom_abline(slope = 1, intercept = 0) +
    geom_point(aes(color = factor(FalsePositiveRate)), size = 2) +
    geom_line(aes(color = factor(FalsePositiveRate)), size = 1) +
    scale_color_brewer('Known Sites\nFalse Positive\nRate', palette = 'PRGn') +
    scale_x_continuous("Predicted Quality") +
    scale_y_continuous("Actual Quality") +
    ggtitle('False Positive Variants Have Little Effect') +
    coord_fixed(ratio = 1) +
    theme_minimal(base_size = 18) + 
    theme(plot.margin = margin(0,0,0,0))
}

plot_comparison_dfs <- function(comparison_dfs){
  ggplot(comparison_dfs, aes(PredictedQuality, ActualQuality)) +
    geom_abline(slope = 1, intercept = 0) +
    geom_point(aes(color = CalibrationMethod), size = 2) +
    geom_line(aes(color = CalibrationMethod), size = 1) +
    scale_color_brewer('Calibration\nMethod', palette = 'Dark2') +
    scale_x_continuous("Predicted Quality") +
    scale_y_continuous('Actual Quality') +
    ggtitle('GATK and KBBQ Perform Similarly') +
    coord_fixed(ratio = 1) +
    theme_minimal(base_size = 18) +
    theme(plot.margin = margin(0,0,0,0))
}



# --- Main --- 

fnr <- c(0, 20, 40, 60, 80, 100)
fpr <- c(0, 20, 40, 60, 80, 100)

fnr_df <- average_over_rgs(import_fnr_csvs(fnr))

pdf('../figures/fnr.pdf', width = 9, height = 7)
plot_fnr_csvs(fnr_df)
dev.off()

kbbq_comparison_df <- get_gatk_csv_file('kbbq') %>%
  read_csv() %>%
  filter(CovariateName == 'QualityScore' & Recalibration == 'After')
gatk_comparison_df <- get_gatk_csv_file('0') %>%
  read_csv() %>%
  filter(CovariateName == 'QualityScore' & Recalibration == 'After')
raw_comparison_df <- get_gatk_csv_file('0') %>%
  read_csv() %>%
  filter(CovariateName == 'QualityScore' & Recalibration == 'Before')

comparison_dfs <- list( KBBQ = kbbq_comparison_df, GATK = gatk_comparison_df, Raw = raw_comparison_df)
comparison_df <- bind_rows(comparison_dfs, .id = 'CalibrationMethod') %>%
  rename(PredictedQuality = CovariateValue) %>%
  mutate(PredictedQuality = as.numeric(PredictedQuality)) %>%
  group_by(CalibrationMethod, PredictedQuality) %>%
  summarize(ActualQuality = p_to_q(sum(Errors)/sum(Observations)))

pdf('../figures/comparison.pdf', width = 9, height = 7)
plot_comparison_dfs(comparison_df)
dev.off()

fpr_df <- average_over_rgs_fpr(import_fpr_csvs(fpr))

pdf('../figures/fpr.pdf', width = 9, height = 7)
plot_fpr_csvs(fpr_df)
dev.off()
    
fnrfpr_df <- average_over_rgs_fnrfpr(import)



