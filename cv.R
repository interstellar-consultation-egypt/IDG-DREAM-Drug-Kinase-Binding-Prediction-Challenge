#===============================================================================

## About this file:
## This file is for performing tasks pertaining experimental design on pairwise
## input methods (particularly DTI prediction). It is assumed that the data
## given to this file (i.e. for training prediction models and testing them) is
## as follows:
##     > Column 1: compound ID
##     > Column 2: target ID
##     > Column 3: response variable
##     > Remaining Columns: compound features followed by target features

#===============================================================================

library(tidyverse)

#===============================================================================

evaluate <- function(predictions, labels, evalMetrics) {
    
    ## type casting to prevent problems later
    predictions <- as.double(predictions)
    labels <- as.double(labels)
    
    ## to hold evaluation results
    evalResults <- vector(mode = 'double', length(evalMetrics))
    
    ## evaluate!
    i <- 0
    for (evalMetric in evalMetrics) {
        i <- i + 1
        
        if (evalMetric == 'rmse')        ## root mean square error
            evalResults[i] <- Metrics::rmse(labels, predictions)
        
        if (evalMetric == 'pearson')     ## pearson correlation coefficient
            evalResults[i] <- cor(labels, predictions)
        
        if (evalMetric == 'spearman')    ## spearman correlation coefficient
            evalResults[i] <- cor(labels, predictions, method = 'spearman')
        
        if (evalMetric == 'f1score')     ## F1 score
            evalResults[i] <- Metrics::f1((labels > 7) %>% as.numeric(), 
                                          (predictions > 7) %>% as.numeric())
        
        if (evalMetric == 'ci')          ## concordance index
            evalResults[i] <- 0
        
        if (evalMetric == 'avg_auc')     ## average AUC
            evalResults[i] <- 0
    
        
        ##########
        ### TO DO: fix eval. metrics: Concordance index (CI), average AUC
        ##########
        
    }
    
    ## return evaluation results
    evalResults
}

#===============================================================================

# Helper function: determines the folds to be used in the k-fold CV experiment
get_folds <- function(allData, cvSetting, k) {
    
    ## distinct compounds
    compound_IDs <- as.vector(unlist(unique(allData[,1])))
    
    ## distinct targets
    target_IDs <- as.vector(unlist(unique(allData[,2])))
  
    if (cvSetting == 'S2') {
        len <- length(compound_IDs)    ## len <- numCompounds
    } else if (cvSetting == 'S3') {
        len <- length(target_IDs)      ## len <- numTargets
    } else {
        len <- nrow(allData)           ## len <- numInstances
    }
  
    rand_ind <- sample(len)    ## randomize order
  
    ## number of folds (i.e. 'k') cannot be greater than 'len'
    k <- if (k > len) len else k
    
    
    ##########
    ### TO DO: display warning (k>len  -->  k=len)
    ##########
    
    
    ## get the different folds of the k-fold experiment
    folds = list()
    for (i in 1:k) {
        if (cvSetting == 'S2') {
            ## if it is leave-drug-out CV, then testCompounds = i
            testCompounds <- if (k == len) i else
                ## otherwise, testCompounds = ...
                rand_ind[(floor((i-1)*len/k)+1) : floor(i*len/k)]
            
            ## retrieve current fold
            folds[[i]] <-
                allData %>% 
                filter(compound_id %in% compound_IDs[testCompounds])
        
        } else if (cvSetting == 'S3') {
            ## if it is leave-target-out CV, then testTargets = i
            testTargets <- if (k == len) i else
                ## otherwise, testTargets = ... 
                rand_ind[(floor((i-1)*len/k)+1) : floor(i*len/k)]
            
            ## retrieve current fold
            folds[[i]] <-
                allData %>% 
                filter(target_id %in% target_IDs[testTargets])
        
        } else {
            ## if it is leave-one-out CV OR leave-both-drug-&-target-out CV,
            ## then testTargets = i
            testInstances <- if (k == len) i else 
                ## otherwise, testTargets = ... 
                rand_ind[(floor((i-1)*len/k)+1) : floor(i*len/k)]
            
            ## retrieve current fold
            folds[[i]] <-
                allData[testInstances,]
        }
    }
    
    ## return folds info
    folds
}

#===============================================================================

kfoldcv <- function(allData, 
                    numCompoundFeatures,
                    cvSetting = 'S1', 
                    k = 5, 
                    poolingCV = T,
                    evalMetrics = c('rmse'),
                    predFun = match.fun('rls_kron'),
                    predFunParams = list(),
                    dimReduction = 'none') {
    
    ## for reproducibility
    set.seed(12345)
    
    ## the k-fold experiment --------------------
    folds <- get_folds(allData, cvSetting, k)    ## get the k folds
    k = length(folds)                            ## in case 'k' was modified
    predictions <- list()                        ## to hold predictions
    for (i in 1:k) {
        ## get test set
        testSet <- folds[[i]]
    
        ## get training set
        if (cvSetting == 'S4') {    ## if it is the special case of S4...
            ## exclude any records containing compounds or targets 
            ## that exist in the test set
            trainingSet <- 
                allData %>% 
                anti_join(testSet, by = 'compound_id') %>% 
                anti_join(testSet, by = 'target_id')
      
        } else {    ## otherwise...
            ## exclude test instances from training set
            trainingSet <- 
                allData %>% 
                anti_join(testSet, by = colnames(testSet))
                ## the 'by = colnames(...)' suppress unwanted messages/warnings
        }
        
        ## TRAINING SET, TRAINING SET LABELS & TEST SET
        trainingSetLabels <- trainingSet[,3] %>% as.matrix()
        trainingSet <- trainingSet[,-(1:3)]
        testSet <- testSet[,-(1:3)]
        
        ## INDICES OF COMPOUND/TARGET FEATURES
        compFeatIndx <- 1:numCompoundFeatures
        targFeatIndx <- (numCompoundFeatures+1):ncol(trainingSet)
        
        ## perform training + prediction for current fold i
        predictions[[i]] <- 
            predFun(trainingSet,
                    trainingSetLabels,
                    testSet,
                    compFeatIndx,
                    targFeatIndx,
                    predFunParams,
                    dimReduction)
    }
    
    
    # ## remove unneeded variables
    # rm(trainingSet, trainingSetLabels, testSet, compFeatIndx, targFeatIndx)
    
    
    ## compute evaluation metrics --------------------
    if (poolingCV) {    ## if pooling CV is to be used...
        ## POOL the predictions...
        predictions <- unlist(predictions)
        labels <- c()
        for (i in 1:k) {
            labels <- c(labels, folds[[i]][,3])
        }
        labels <- unlist(labels)
    
        ## ...then get evaluation results
        evalResults <- evaluate(predictions, 
                                labels, 
                                evalMetrics)
    
    } else {    ## average CV
        ## to hold evaluation results
        evalResults <- matrix(0, 
                              nrow = length(folds), 
                              ncol = length(evalMetrics))
    
        ## Get evaluation results for each fold...
        for (i in 1:k) {
            labels <- folds[[i]][,3] %>% as.matrix()
            evalResults[i,] <- evaluate(predictions[[i]], 
                                        labels, 
                                        evalMetrics)
        }
    
        ## ...then AVERAGE them
        # evalResultsSD <- apply(evalResults, 2, sd)
        evalResults <- apply(evalResults, 2, mean)
    }
    
    ## return evaluation results
    evalResults
}

#===============================================================================

loocv <- function(allData,
                  numCompoundFeatures,
                  cvSetting = 'S1', 
                  poolingCV = T,
                  evalMetrics = c('rmse'),
                  predFun = match.fun('rls_kron'),
                  predFunParams = list(),
                  dimReduction = 'none') {
    
    ## determine value of 'k' based on 'cvSetting'
    if (cvSetting == 'S2') {               ## leave-drug-out
        k = nrow(unique(allData[,1]))
    } else if (cvSetting == 'S3') {        ## leave-target-out
        k = nrow(unique(allData[,2]))
    } else {    ## traditional leave-one-out OR leave-both-drug-&-target-out
        k = nrow(allData)
    }
    
    ## get evaluation results
    evalResults <- kfoldcv(allData, 
                           numCompoundFeatures,
                           cvSetting, 
                           k, 
                           poolingCV, 
                           evalMetrics, 
                           predFun, 
                           predFunParams,
                           dimReduction)
    
    ## return evaluation results
    evalResults
}

#===============================================================================