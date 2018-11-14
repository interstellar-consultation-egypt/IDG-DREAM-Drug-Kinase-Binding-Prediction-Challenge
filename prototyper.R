

## load experimental design & machine learning functions
library(tidyverse)
source('cv.R')
source('ml.R')


## ------------------------------------------

## FULL CV SETTINGS
cvSettings <- c('S1', 'S2', 'S3', 'S4')
loos <- c(F,T)     ## leave one out?
k <- 10            ## number of folds (used only when 'loo = F')
pools <- c(F,T)    ## use pooling CV?
evalMetrics <- c('rmse', 'pearson', 'spearman', 'f1score', 'ci', 'avg_auc')

## ---------------

## PARTIAL CV SETTINGS    [OVERRIDES FULL CV SETTINGS]
evalMetrics <- c('rmse', 'pearson', 'spearman', 'f1score')
loos <- F          ## k-fold cv 
pools <- F         ## averaging cv

##########
### TO DO:   cv.R --> fix eval. metrics: Concordance index (CI), average AUC
##########

## ------------------------------------------


## load data
path <- 'data/dataset1/'
train_dta <- read_csv(paste(path, 'train/train.csv', sep = ''))
numCompoundFeatures <- 165


## ------------------------------------------


##############
## rls_kron ##
##############

currentMethod <- 'rls_kron'
predFun <- match.fun(currentMethod)
for (loo in loos) {
    for (cvSetting in cvSettings) {
        for (pool in pools) {
            
            looString  <- if (loo)  'Y' else 'N'
            poolString <- if (pool) 'Y' else 'N'
            
            cat('\n====================================================\n\n')
            cat('CV Setting:\t',    cvSetting,  '\n')
            cat('Leave One Out:\t', looString)
            if (!loo)
                cat('\t[k=', k, ']', sep = '')
            cat('\n')
            cat('Pooling CV:\t',    poolString, '\n')
            cat('Start Time:\t',    as.character(Sys.time()), '\n\n\n')
            
            ## print header of evaluation results table
            cat('Current Method:\t', currentMethod, '\n')
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

            cat('\n====================================================\n\n')
            
        }
    }
}

cat('\n################################################################')
cat('\n################################################################')
cat('\n################################################################')
cat('\n################################################################\n\n')


#############
## rls_avg ##
#############

currentMethod <- 'rls_avg'
predFun <- match.fun(currentMethod)
for (loo in loos) {
    for (cvSetting in cvSettings) {
        for (pool in pools) {
            
            looString  <- if (loo)  'Y' else 'N'
            poolString <- if (pool) 'Y' else 'N'
            
            cat('\n====================================================\n\n')
            cat('CV Setting:\t',    cvSetting,  '\n')
            cat('Leave One Out:\t', looString)
            if (!loo)
                cat('\t[k=', k, ']', sep = '')
            cat('\n')
            cat('Pooling CV:\t',    poolString, '\n')
            cat('Start Time:\t',    as.character(Sys.time()), '\n\n\n')

            ## print header of evaluation results table
            cat('Current Method:\t', currentMethod, '\n')
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

            cat('\n====================================================\n\n')
            
        }
    }
}

cat('\n################################################################')
cat('\n################################################################')
cat('\n################################################################')
cat('\n################################################################\n\n')


##############
## EnsemKRR ##
##############

currentMethod <- 'ensembler'
predFun <- match.fun(currentMethod)
for (loo in loos) {
    for (cvSetting in cvSettings) {
        for (pool in pools) {
            
            looString  <- if (loo)  'Y' else 'N'
            poolString <- if (pool) 'Y' else 'N'
            
            cat('\n====================================================\n\n')
            cat('CV Setting:\t',    cvSetting,  '\n')
            cat('Leave One Out:\t', looString)
            if (!loo)
                cat('\t[k=', k, ']', sep = '')
            cat('\n')
            cat('Pooling CV:\t',    poolString, '\n')
            cat('Start Time:\t',    as.character(Sys.time()), '\n\n\n')

            ## print header of evaluation results table
            cat('Current Method:\t', currentMethod, '\n')
            cat('lambda\tnmLrnr\tbag\tr\tbsLrnr')
            for (em in evalMetrics)
                cat('\t', substr(em, start = 1, stop = 6))
            cat('\n')

            ## grid search over parameters: lambda
            results <- expand.grid(baseLearner = c('rls_avg'),
                                   bag = F,
                                   numLearners = 20,
                                   r = seq(0.1,0.9,0.1),
                                   lambda=(2^(-5:-2)),
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

cat('\n################################################################')
cat('\n################################################################')
cat('\n################################################################')
cat('\n################################################################\n\n')


##################
## EnsemRLSkron ##
##################

currentMethod <- 'ensembler'
predFun <- match.fun(currentMethod)
for (loo in loos) {
    for (cvSetting in cvSettings) {
        for (pool in pools) {
            
            looString  <- if (loo)  'Y' else 'N'
            poolString <- if (pool) 'Y' else 'N'
            
            cat('\n====================================================\n\n')
            cat('CV Setting:\t',    cvSetting,  '\n')
            cat('Leave One Out:\t', looString)
            if (!loo)
                cat('\t[k=', k, ']', sep = '')
            cat('\n')
            cat('Pooling CV:\t',    poolString, '\n')
            cat('Start Time:\t',    as.character(Sys.time()), '\n\n\n')

            ## print header of evaluation results table
            cat('Current Method:\t', currentMethod, '\n')
            cat('lambda\tnmLrnr\tbag\tr\tbsLrnr')
            for (em in evalMetrics)
                cat('\t', substr(em, start = 1, stop = 6))
            cat('\n')

            ## grid search over parameters: lambda
            results <- expand.grid(baseLearner = c('rls_kron'),
                                   bag = F,
                                   numLearners = 20,
                                   r = seq(0.1,0.9,0.1),
                                   lambda=(2^(-5:-2)),
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

cat('\n################################################################')
cat('\n################################################################')
cat('\n################################################################')
cat('\n################################################################\n\n')


#######################
## rls_kron_graphreg ##
#######################

currentMethod <- 'rls_kron_graphreg'
predFun <- match.fun(currentMethod)
for (loo in loos) {
    for (cvSetting in cvSettings) {
        for (pool in pools) {
            
            looString  <- if (loo)  'Y' else 'N'
            poolString <- if (pool) 'Y' else 'N'
            
            cat('\n====================================================\n\n')
            cat('CV Setting:\t',    cvSetting,  '\n')
            cat('Leave One Out:\t', looString)
            if (!loo)
                cat('\t[k=', k, ']', sep = '')
            cat('\n')
            cat('Pooling CV:\t',    poolString, '\n')
            cat('Start Time:\t',    as.character(Sys.time()), '\n\n\n')

            ## print header of evaluation results table
            cat('Current Method:\t', currentMethod, '\n')
            cat('lambda\tp')
            for (em in evalMetrics)
                cat('\t', substr(em, start = 1, stop = 6))
            cat('\n')

            ## grid search over parameters: lambda, p
            results <- expand.grid(lambda=(2^(-5:5)), p=1:10)
            results <-
                results %>%
                cbind(matrix(0, nrow(results), length(evalMetrics)))
            colnames(results)[3:ncol(results)] <- evalMetrics
            predFunParams <- list()
            for (i in 1:nrow(results)) {
                predFunParams$lambda <- results$lambda[i]
                predFunParams$p <- results$p[i]

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
                cat(round(predFunParams$lambda,3), '\t', predFunParams$p)
                for (j in 1:length(evalMetrics)) {
                    cat('\t', round(evalResults[j], 3))
                }

                cat('\t', as.character(Sys.time()), '\n')
            }

            cat('\n====================================================\n\n')
            
        }
    }
}

cat('\n################################################################')
cat('\n################################################################')
cat('\n################################################################')
cat('\n################################################################\n\n')

