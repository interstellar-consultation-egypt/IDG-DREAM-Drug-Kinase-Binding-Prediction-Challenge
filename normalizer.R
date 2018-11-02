library(tidyverse)
library(BBmisc)


## read compound features
train_compound_features <- read_csv("data/train/compoundFeatures.csv")
test_compound_features <- read_csv("data/test/compoundFeatures.csv")

## add compound_IDs column to test data
tmpDf <- data.frame(compound_IDs = rep(NA, nrow(test_compound_features)))
test_compound_features <- cbind(tmpDf, test_compound_features)

## combine train and test data sets together
compound_features <- rbind(train_compound_features, test_compound_features)
# compound_features <- full_join(train_compound_features, test_compound_features)

## normalization + write to file
compound_features_normalized <- normalize(compound_features, method = "range")
write_csv(subset(compound_features_normalized, !is.na((compound_IDs))), 
          "data/train/compoundFeaturesNormalized.csv")
write_csv(subset(compound_features_normalized, is.na((compound_IDs))), 
          "data/test/compoundFeaturesNormalized.csv")


## clear environment
rm(list = ls())


## read target features
train_target_features <- read_csv("data/train/targetFeatures.csv")
test_target_features <- read_csv("data/test/targetFeatures.csv")

## add target_IDs column to test data
tmpDf <- data.frame(target_IDs = rep(NA, nrow(test_target_features)))
test_target_features <- cbind(tmpDf, test_target_features)

## combine train and test data sets together
target_features <- rbind(train_target_features, test_target_features)
# target_features <- full_join(train_target_features, test_target_features)

## normalization
target_features_normalized <- target_features
target_features_normalized[,2:420]    <- 
    target_features_normalized[,2:420]    / 100
target_features_normalized[,691:1189] <- 
    target_features_normalized[,691:1189] / 100
remainingFeatures <- 1:ncol(target_features_normalized)
remainingFeatures <- setdiff(remainingFeatures, c(2:420, 691:1189))
target_features_normalized[,remainingFeatures] <- 
    normalize(target_features_normalized[,remainingFeatures], method = "range")

## write to file
write_csv(subset(target_features_normalized, !is.na((target_IDs))), 
          "data/train/targetFeaturesNormalized.csv")
write_csv(subset(target_features_normalized, is.na((target_IDs))), 
          "data/test/targetFeaturesNormalized.csv")
