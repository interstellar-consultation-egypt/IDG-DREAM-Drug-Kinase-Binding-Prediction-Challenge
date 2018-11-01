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
