

## path of the training and test datasets to be created
path <- 'data/dataset2/'


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


library(tidyverse)
library(BBmisc)


## read compound features
colTypesSpecs <- paste('c', paste(rep('d',1426), collapse=''), sep='')
train_compound_features <- 
    read_csv(paste(path, 'train/compoundFeatures.csv', sep = ''), 
             col_types = colTypesSpecs)
colTypesSpecs <- paste(rep('d',1426), collapse='')
test_compound_features <- 
    read_csv(paste(path, 'test/compoundFeatures.csv', sep = ''), 
             col_types = colTypesSpecs)

## add compound_IDs column to test data
tmpDf <- data.frame(compound_id = rep(NA, nrow(test_compound_features)))
test_compound_features <- cbind(tmpDf, test_compound_features)

## combine train and test data sets together
compound_features <- rbind(train_compound_features, test_compound_features)
# compound_features <- 
#     full_join(train_compound_features, test_compound_features)

## fill NA's with column mean (if any)
replaceWithMean <- function(x) {
    x <- replace_na(x, mean(unique(x), na.rm = TRUE))
}
compound_features[,-1] <- map_df(compound_features[,-1], replaceWithMean)

## normalization + write to file
## Features post-processing: scale to range [0,1] using min-max normalization
compound_features_normalized <- normalize(compound_features, method = 'range')
write_csv(subset(compound_features_normalized, !is.na((compound_id))), 
          paste(path, 'train/compoundFeaturesNormalized.csv', sep = ''))
write_csv(subset(compound_features_normalized, is.na((compound_id))), 
          paste(path, 'test/compoundFeaturesNormalized.csv', sep = ''))


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


## read target features
train_target_features <- 
    read_csv(paste(path, 'train/targetFeatures.csv', sep = ''),
             na = character())
test_target_features <- 
    read_csv(paste(path, 'test/targetFeatures.csv', sep = ''),
             na = character())

## add target_IDs column to test data
tmpDf <- data.frame(target_IDs = rep(NA, nrow(test_target_features)))
test_target_features <- cbind(tmpDf, test_target_features)

## combine train and test data sets together
target_features <- rbind(train_target_features, test_target_features)
# target_features <- full_join(train_target_features, test_target_features)

## fill NA's with column mean (if any)
replaceWithMean <- function(x) {
    x <- replace_na(x, mean(unique(x), na.rm = TRUE))
}
target_features[,-1] <- map_df(target_features[,-1], replaceWithMean)

## normalization
## Features post-processing: percentage-style variables to be divided by 100
## These are the variables starting with '[G1.', '[G2.' and '[G4.'
target_features_normalized <- target_features
target_features_normalized[,2:568]    <- 
    target_features_normalized[,2:568] / 100

## normalization (remaining variables) + write to file
## Features post-processing: scale to range [0,1] using min-max normalization
remainingFeatures <- 1:ncol(target_features_normalized)
remainingFeatures <- setdiff(remainingFeatures, c(2:568))
target_features_normalized[,remainingFeatures] <- 
    normalize(target_features_normalized[,remainingFeatures], method = 'range')

## write to file
write_csv(subset(target_features_normalized, !is.na((target_IDs))), 
          paste(path, 'train/targetFeaturesNormalized.csv', sep = ''))
write_csv(subset(target_features_normalized, is.na((target_IDs))), 
          paste(path, 'test/targetFeaturesNormalized.csv', sep = ''))
