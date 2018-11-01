library(tidyverse)
library(BBmisc)

train_compound_features <- read_csv("data/train/compoundFeatures.csv")
test_compound_features <- read_csv("data/test/compoundFeatures.csv")
compound_features <- full_join(train_compound_features, test_compound_features)
compound_features_normalized <- normalize(compound_features, method = "range")
write_csv(subset(compound_features_normalized, !is.na((compound_IDs))), 
          "data/train/compoundFeaturesNormalized.csv")
write_csv(subset(compound_features_normalized, is.na((compound_IDs))), 
          "data/test/compoundFeaturesNormalized.csv")

train_target_features <- read_csv("data/train/targetFeatures.csv")
test_target_features <- read_csv("data/test/targetFeatures.csv")
target_features <- full_join(train_compound_features, test_target_features)
target_features_normalized <- normalize(target_features, method = "range")
write_csv(subset(target_features_normalized, !is.na((target_IDs))), 
          "data/train/targetFeaturesNormalized.csv")
write_csv(subset(target_features_normalized, is.na((target_IDs))), 
          "data/test/targetFeaturesNormalized.csv")
