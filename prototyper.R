

## load experimental design & machine learning functions
library(tidyverse)
source('cv.R')
source('ml.R')


## FULL CV SETTINGS
cvSettings <- c('S1', 'S2', 'S3', 'S4')
loos <- c(T,F)     ## leave one out?
k <- 5             ## number of folds (used only when 'loo = F')
pools <- c(T,F)    ## use pooling CV?
evalMetrics <- c('rmse', 'pearson', 'spearman', 'f1score')

## ------------------------------------------

## PARTIAL CV SETTINGS    [OVERRIDES FULL CV SETTINGS]
cvSettings <- c('S2', 'S4')
loos <- F
k <- 10
pools <- T

## ------------------------------------------


##########
### TO DO: remaining eval metrics  -->  Concordance index (CI), average AUC
##########


## load data
train_dta <- read_csv('data/train/train.csv')
numCompoundFeatures <- 165


## start the whole thing!
for (cvSetting in cvSettings) {
    for (loo in loos) {
        for (pool in pools) {
            
            looString  <- if (loo)  'Y' else 'N'
            poolString <- if (pool) 'Y' else 'N'
            
            cat('\n====================================================\n\n')
            cat('CV Setting:\t',    cvSetting,  '\n')
            cat('Leave One Out:\t', looString,  '\n')
            cat('Pooling CV:\t',    poolString, '\n')
            cat('Start Time:\t',    as.character(Sys.time()), '\n\n\n')
            
                
            ##############
            ## rls_kron ##
            ##############

            currentMethod <- 'rls_kron'
            predFun <- match.fun(currentMethod)
            cat('Current Method:\t', currentMethod, '\n\n')

            ## print header of evaluation results table
            cat('lambda')
            for (em in evalMetrics)
                cat('\t', substr(em, start = 1, stop = 6))
            cat('\n')

            ## grid search over parameters: lambda
            results <- expand.grid(lambda=(2^(-5:5)))
            results <-
                results %>%
                cbind(matrix(0, nrow(results), length(evalMetrics)))
            colnames(results)[3:ncol(results)] <- evalMetrics
            predFunParams <- list()
            for (i in 1:nrow(results)) {
                predFunParams$lambda <- results$lambda[i]

                ## get evaluation results
                if (loo) {
                    evalResults <- loocv(train_dta,
                                         numCompoundFeatures,
                                         cvSetting,
                                         pool,
                                         evalMetrics,
                                         predFun,
                                         predFunParams)

                } else {
                    evalResults <- kfoldcv(train_dta,
                                           numCompoundFeatures,
                                           cvSetting,
                                           k,
                                           pool,
                                           evalMetrics,
                                           predFun,
                                           predFunParams)
                }

                ## print parameter values + their evaluation results
                cat(round(predFunParams$lambda,3))
                for (j in 1:length(evalMetrics)) {
                    cat('\t', round(evalResults[j], 3))
                }

                cat('\t', as.character(Sys.time()), '\n')
            }

            cat('\n')
            
                
            # #############
            # ## rls_avg ##
            # #############
            # 
            # currentMethod <- 'rls_kron'
            # predFun <- match.fun(currentMethod)
            # cat('Current Method:\t', currentMethod, '\n\n')
            # 
            # ## print header of evaluation results table
            # cat('lambda')
            # for (em in evalMetrics)
            #     cat('\t', substr(em, start = 1, stop = 6))
            # cat('\n')
            # 
            # ## grid search over parameters: lambda
            # results <- expand.grid(lambda=(2^(-5:5)))
            # results <- 
            #     results %>% 
            #     cbind(matrix(0, nrow(results), length(evalMetrics)))
            # colnames(results)[3:ncol(results)] <- evalMetrics
            # predFunParams <- list()
            # for (i in 1:nrow(results)) {
            #     predFunParams$lambda <- results$lambda[i]
            #     
            #     ## get evaluation results
            #     if (loo) {
            #         evalResults <- loocv(train_dta,
            #                              numCompoundFeatures,
            #                              cvSetting,
            #                              pool,
            #                              evalMetrics,
            #                              predFun,
            #                              predFunParams)
            #         
            #     } else {
            #         evalResults <- kfoldcv(train_dta,
            #                                numCompoundFeatures,
            #                                cvSetting,
            #                                k,
            #                                pool,
            #                                evalMetrics,
            #                                predFun,
            #                                predFunParams)
            #     }
            #     
            #     ## print parameter values + their evaluation results
            #     cat(round(predFunParams$lambda,3))
            #     for (j in 1:length(evalMetrics)) {
            #         cat('\t', round(evalResults[j], 3))
            #     }
            #     
            #     cat('\t', as.character(Sys.time()), '\n')
            # }
            # 
            # cat('\n')
            
                
            ##############
            ## EnsemKRR ##
            ##############
            
            currentMethod <- 'ensembler'
            predFun <- match.fun(currentMethod)
            cat('Current Method:\t', currentMethod, '\n\n')
            
            ## print header of evaluation results table
            cat('lambda\tnmLrnr\tbag\tr\tbsLrnr')
            for (em in evalMetrics)
                cat('\t', substr(em, start = 1, stop = 6))
            cat('\n')
            
            ## grid search over parameters: lambda
            results <- expand.grid(lambda=0.125,
                                   numLearners = c(20,50,100,500),
                                   bag = F,
                                   r = 0.7,    #seq(0.1,0.9,0.1),
                                   baseLearner = 'rls_kron',
                                   stringsAsFactors = F)
            results <- 
                results %>% 
                cbind(matrix(0, nrow(results), length(evalMetrics)))
            predFunParams <- list()
            for (i in 1:nrow(results)) {
                predFunParams$lambda      <- results$lambda[i]
                predFunParams$numLearners <- results$numLearners[i]
                predFunParams$bag         <- results$bag[i]
                predFunParams$r           <- results$r[i]
                predFunParams$baseLearner <- results$baseLearner[i]
                
                ## get evaluation results
                if (loo) {
                    evalResults <- loocv(train_dta,
                                         numCompoundFeatures,
                                         cvSetting,
                                         pool,
                                         evalMetrics,
                                         predFun,
                                         predFunParams)
                    
                } else {
                    evalResults <- kfoldcv(train_dta,
                                           numCompoundFeatures,
                                           cvSetting,
                                           k,
                                           pool,
                                           evalMetrics,
                                           predFun,
                                           predFunParams)
                }
                
                ## print parameter values + their evaluation results
                cat(round(predFunParams$lambda,3), '\t',
                    predFunParams$numLearners, '\t',
                    predFunParams$bag, '\t',
                    predFunParams$r, '\t',
                    substr(predFunParams$baseLearner, start = 1, stop = 6))
                for (j in 1:length(evalMetrics)) {
                    cat('\t', round(evalResults[j], 3))
                }
                
                cat('\t', as.character(Sys.time()), '\n')
            }
            
            cat('\n====================================================\n\n')
            
        }
    }
}
