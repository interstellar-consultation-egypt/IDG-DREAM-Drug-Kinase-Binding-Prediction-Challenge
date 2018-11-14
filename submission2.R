

## load necessary packages
library(tidyverse)
source('ml.R')


## load data
path <- 'data/dataset1/'
train_dta  <- read_csv(paste(path, 'train/train.csv', sep = ''))
test_dta   <- read_csv(paste(path, 'test/test.csv', sep = ''))
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
                         targFeatIndx,
                         predFunParams = list(numLearners = 20,
                                              bag = F,
                                              r = 0.7,
                                              baseLearner = 'rls_kron_graphreg',
                                              lambda = 2^-6,
                                              p = 1))
predictions <- 
    predictions %>% 
    as_data_frame() %>% 
    rename(`pKd_[M]_pred` = pkd)


## SUBMISSION FILE!
submission %>% 
    cbind(predictions) %>% 
    write_csv('data/round_1_submission_2.csv')


## RESULTS
## RMSE     PEARSON     SPEARMAN    CI          F1          AVG_AUC
## 1.6785   -0.0396     -0.0544     0.5054      0.1978      0.4813

