

## load necessary packages
library(tidyverse)
source('ml.R')


## load data
train_dta <- read_csv('data/train/train.csv')
test_dta  <- read_csv('data/test/test.csv')
submission <- read_csv('data/round_1_template.csv')
numCompoundFeatures <- 165


## TRAINING SET, TRAINING SET LABELS & TEST SET
trainingSetLabels <- train_dta[,3] %>% as.matrix()
train_dta <- train_dta[,-(1:3)]

        
## INDICES OF COMPOUND/TARGET FEATURES
compFeatIndx <- 1:numCompoundFeatures
targFeatIndx <- (numCompoundFeatures+1):ncol(train_dta)


## PREDICTIONS
predictions <- ensembler(train_dta,
                         trainingSetLabels,
                         test_dta,
                         compFeatIndx,
                         targFeatIndx)
predictions <- 
    predictions %>% 
    as_data_frame() %>% 
    rename(`pKd_[M]_pred` = pkd)


## SUBMISSION FILE!
submission %>% 
    cbind(predictions) %>% 
    write_csv('data/round_1_submission_1.csv')

