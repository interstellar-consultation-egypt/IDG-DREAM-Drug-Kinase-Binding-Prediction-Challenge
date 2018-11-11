

## load necessary packages
library(tidyverse)
source('ml.R')


## load data
path <- 'data/dataset1/'
train_dta  <- read_csv(paste(path, 'train/train.csv', sep = ''))
test_dta   <- read_csv(paste(path, 'test/test.csv', sep = ''))
submission <- read_csv(paste(path, 'round_1_template.csv', sep = ''))
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


## RESULTS
## RMSE     PEARSON     SPEARMAN    CI          F1          AVG_AUC
## 1.1816	0.2012	    0.1898	    0.5566	    0.1194	    0.579

