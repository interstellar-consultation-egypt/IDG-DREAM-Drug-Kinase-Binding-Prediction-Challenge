library(dplyr)

#===============================================================================

evaluate <- function(predictions, labels, evalMetrics) {
    evalResults <- vector(mode = 'double', length(evalMetrics))
    i <- 0
    for (evalMetric in evalMetrics) {
        i <- i + 1
        if (evalMetric == 'rmse')
            evalResults[i] <- Metrics::rmse(labels, predictions)
        if (evalMetric == 'pearson')
            evalResults[i] <- cor(labels, predictions)
        if (evalMetric == 'spearman')
            evalResults[i] <- cor(labels, predictions, method = 'spearman')
        if (evalMetric == 'f1score')
            evalResults[i] <- Metrics::f1(labels, predictions)
    
        ### TO DO: remaining metrics  -->  Concordance index (CI), AUC
    }
}

#===============================================================================

get_folds <- function(allData, cvSetting, k) {
    compound_IDs <- unique(allData$compound_id)    ## distinct compounds
    target_IDs <- unique(allData$target_id)        ## distinct targets
  
    if (cvSetting == 'S2') {
        len <- length(compound_IDs)    ## len <- numDrugs
    } else if (cvSetting == 'S3') {
        len <- length(target_IDs)      ## len <- numTargets
    } else {
        len <- nrow(allData)           ## len <- numInstances
    }
  
    set.seed(12345)            ## for reproducibility
    rand_ind <- sample(len)    ## randomize order
  
    ## number of folds (i.e. 'k') cannot be greater than 'len'
    k <- if (k > len) len else k
    
    ### TO DO: display warning (k>len  -->  k=len)
  
    ## get the different folds of the k-fold experiment
    folds = list()
    for (i in 1:k) {
        if (cvSetting == 'S2') {
            ## if it is leave-drug-out CV, then testCompounds = i
            testCompounds <- if (k == len) i else
                ## otherwise, testCompounds = ...
                rand_ind[(floor((i-1)*len/k)+1) : floor(i*len/k)]
            testCompounds <- compound_IDs(testCompounds)
            folds[[i]] <-
                allData %>% 
                filter(compound_id %in% testCompounds)
        
        } else if (cvSetting == 'S3') {
            ## if it is leave-target-out CV, then testTargets = i
            testTargets <- if (k == len) i else
                ## otherwise, testTargets = ... 
                rand_ind[(floor((i-1)*len/k)+1) : floor(i*len/k)]
            testTargets <- target_IDs(testTargets)
            folds[[i]] <-
                allData %>% 
                filter(target_id %in% testTargets)
        
        } else {
            testInstances <- if (k == len) i else 
                rand_ind[(floor((i-1)*len/k)+1) : floor(i*len/k)]
            folds[[i]] <-
                allData[testInstances,]
        }
    }
}

#===============================================================================

kfoldcv <- function(allData, 
                    cvSetting = 'S1', 
                    k = 5, 
                    poolingCV = T,
                    evalMetrics = c('rmse'),
                    predFun,
                    predFunParams = list()) {
  
  
    ## the k-fold experiment --------------------
    folds <- get_folds(allData, cvSetting, k)    ## get the k folds
    k = length(folds)                            ## in case 'k' was modified
    predictions <- list()
    for (i in 1:k) {
        ## get test set
        testSet <- folds[[i]]
    
        ## get training set
        if (cvSetting == 'S4') {    ## if it is the special case of S4...
            ## exclude records containing compounds/targets existing in test set
            trainingSet <- 
                allData %>% 
                anti_join(testSet, by = 'compound_id') %>% 
                anti_join(testSet, by = 'target_id')
      
        } else {    ## otherwise...
            ## exclude test instances from training set
            trainingSet <- 
                allData %>% 
                anti_join(testSet)
        }
    
        ## perform training + prediction for current fold i
        predictions[[i]] <- predFun(trainingSet, 
                                    testSet, 
                                    predFunParams)
    }
    
    
    ## compute evaluation metrics --------------------
    if (poolingCV) {    ## if pooling CV is to be used...
        ## POOL the predictions...
        predictions <- unlist(predictions)
        labels <- c()
        for (i in 1:k) {
            labels <- c(labels, folds[[i]]$pkd)
        }
    
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
            evalResults[i,] <- evaluate(predictions, 
                                        labels, 
                                        evalMetrics)
        }
    
        ## ...then AVERAGE them
        evalResults <- apply(evalResults, 2, mean)
    }
  
    evalResults
}

#===============================================================================

loocv <- function(allData,
                  cvSetting = 'S1', 
                  poolingCV = T,
                  evalMetrics = c('rmse'),
                  predFun,
                  predFunParams = list()) {
    ## allData = dataframe containing all the data
    ## cvSetting = type of CV being considered (i.e. one of 'S1','S2','S3','S4')
    ## poolingCV = pool results together before evaluation
    ## evalMetrics = evaluation metrics to compute
    ## predFun = prediction model function to be called (training + prediction)
    ## predFunParams = any parameters the prediction model may need
  
  
    ## determine function to use based on 'cvSetting'
    if (cvSetting == 'S2') {    ## leav-drug-out
        k = length(unique(allData$compound_id))
    } else if (cvSetting == 'S3') {    ## leav-target-out
        k = length(unique(allData$target_id))
    } else {    ## traditional leave-one-out OR leave-both-drug-&-target-out
        k = nrow(allData)
    }
  
  
    ## get evaluation results
    evalResults <- kfoldcv(allData, 
                           cvSetting, 
                           k, 
                           poolingCV, 
                           evalMetrics, 
                           predFun, 
                           predFunParams)
  
  
    evalResults
}

#===============================================================================