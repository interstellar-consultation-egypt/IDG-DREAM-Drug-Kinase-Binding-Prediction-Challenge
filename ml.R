#===============================================================================

## Assumptions
## 'trainingSet' and 'testSet' have the same number of columns
## The columns consist of compound features followed by target features

#===============================================================================


rls_avg <- function(trainingSet,
                    trainingSetLabels,
                    testSet,
                    compFeatIndx,
                    targFeatIndx,
                    predFunParams = list(lambda = 0.125),
                    dimReduction = 'none') {
    
    ## PARAMETER VALUES
    lambda <- predFunParams$lambda
    
    ## DIMENSIONALITY REDUCTION
    if (dimReduction != 'none') {
        dim_red(trainingSet,
                testSet,
                compFeatIndx,
                targFeatIndx,
                dimReduction)
    }
    
    ## KERNELS
    Kc_train <- rbf_kernel(as.matrix(trainingSet[,compFeatIndx]))
    Kt_train <- rbf_kernel(as.matrix(trainingSet[,targFeatIndx]))
    Kc_test  <- rbf_kernel(as.matrix(testSet[,compFeatIndx]),
                           as.matrix(trainingSet[,compFeatIndx]))
    Kt_test  <- rbf_kernel(as.matrix(testSet[,targFeatIndx]),
                           as.matrix(trainingSet[,targFeatIndx]))
    
    ## PREDICTIONS
          y <- trainingSetLabels
    lambdaI <- lambda * diag(ncol(Kc_train))    ## lambda * identity matrix
     yhat_c <- Kc_test %*% (solve(Kc_train + lambdaI) %*% y)
     yhat_t <- Kt_test %*% (solve(Kt_train + lambdaI) %*% y)
       yhat <- (yhat_c + yhat_t) / 2
    
    ## return predicted labels
    yhat
}


#===============================================================================


rls_kron <- function(trainingSet,
                     trainingSetLabels,
                     testSet,
                     compFeatIndx,
                     targFeatIndx,
                     predFunParams = list(lambda = 0.125),
                     dimReduction = 'none') {
    
    ## PARAMETER VALUES
    lambda <- predFunParams$lambda
    
    ## DIMENSIONALITY REDUCTION
    if (dimReduction != 'none') {
        dim_red(trainingSet,
                testSet,
                compFeatIndx,
                targFeatIndx,
                dimReduction)
    }
    
    ## KERNELS
    Kc_train <- rbf_kernel(as.matrix(trainingSet[,compFeatIndx]))
    Kt_train <- rbf_kernel(as.matrix(trainingSet[,targFeatIndx]))
    Kc_test  <- rbf_kernel(as.matrix(testSet[,compFeatIndx]),
                           as.matrix(trainingSet[,compFeatIndx]))
    Kt_test  <- rbf_kernel(as.matrix(testSet[,targFeatIndx]),
                           as.matrix(trainingSet[,targFeatIndx]))
    Ktrain   <- Kc_train * Kt_train
    Ktest    <- Kc_test * Kt_test
    
    ## PREDICTIONS
          y <- trainingSetLabels
    lambdaI <- lambda * diag(ncol(Kc_train))    ## lambda * identity matrix
       yhat <- Ktest %*% (solve(Ktrain + lambdaI) %*% y)
    
    ## return predicted labels
    yhat
}


#===============================================================================


ensembler <- function(trainingSet,
                      trainingSetLabels,
                      testSet,
                      compFeatIndx,
                      targFeatIndx,
                      predFunParams = list(numLearners = 20,
                                           bag = F,
                                           r = 0.7,
                                           baseLearner = 'rls_kron',
                                           lambda = 0.125),
                      dimReduction = 'none') {
    
    ## PARAMETER VALUES
    numLearners <- predFunParams$numLearners
            bag <- predFunParams$bag
              r <- predFunParams$r
    baseLearner <- match.fun(predFunParams$baseLearner)
         lambda <- predFunParams$lambda
    
    ## ENSEMBLE TIME!
    yhat <- 0
    for (i in 1:numLearners) {
        trainingSet_i <- trainingSet
        trainingSetLabels_i <- trainingSetLabels
        testSet_i <- testSet
        
        ## BAGGING
        if (bag) {
            baggedInstances <- sample(nrow(trainingSet), replace = T)
            trainingSet_i <- trainingSet[baggedInstances,]
            trainingSetLabels_i <- trainingSetLabels[baggedInstances]
        }
        
        ## FEATURE SUBSPACING
        numCompoundFeatures <- length(compFeatIndx)
        numTargetFeatures <- length(targFeatIndx)
        if (r < 1) {
            ## COMPOUND FEATURES
            compoundFeatures <- 
                sample(numCompoundFeatures, floor(numCompoundFeatures*r))
            numCompoundFeatures <- length(compoundFeatures)
            
            ## TARGET FEATURES
            targetFeatures <- 
                sample(numTargetFeatures, floor(numTargetFeatures*r))
            numTargetFeatures <- length(targetFeatures)
            
            ## SET OF SUBSPACED FEATURES
            selectedFeatures <- 
                c(compoundFeatures, (targetFeatures + numCompoundFeatures))
            
            ## MODIFIED TRAINING AND TEST SETS
            trainingSet_i <- trainingSet_i[, selectedFeatures]
            testSet_i <- testSet_i[, selectedFeatures]
        }
        
        ## PREDICTIONS
        yhat <- yhat + baseLearner(trainingSet,
                                   trainingSetLabels,
                                   testSet,
                                   1:numCompoundFeatures,
                                   (numCompoundFeatures+1):ncol(trainingSet),
                                   predFunParams,
                                   dimReduction)
    }
    
    ## return predicted labels
    yhat <- yhat / numLearners
    yhat
}


#===============================================================================
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#===============================================================================


rbf_kernel <- function(X, Y, sigma) {
    ## if argument 'Y' was not passed to the function...
    if (missing(Y))
        Y <- X
    
    ## SIGMA
    if (missing(sigma))     ## if 'sigma' was not passed to the function...
        sigma <- ncol(X)    ## sigma = number of features
    
    ## compute RBF kernel
    rbfKernel <- rdetools::rbfkernel(X, sigma, Y)
    
    ## return RBF kernel
    rbfKernel
}


#===============================================================================


dim_red <- function(trainingSet,
                    testSet,
                    compFeatIndx,
                    targFeatIndx,
                    dimReduction) {
    
    ## type this in the console and hit Enter:
    ## dimRed::dimRedMethodList()
    ##
    ## The output is a list of dimensionality reduction that are available
}


#===============================================================================


# check_missing <- function(trainingSet, 
#                           trainingSetLabels,
#                           testSet, 
#                           numCompoundFeat) {
#  
#     ## if any important arguments missing...   
#     if (missing(trainingSet) || 
#         missing(trainingSetLabels) || 
#         missing(testSet) || 
#         missing(numCompoundFeat))
#         ## STOP!!!
#         stop('function \'ml.R::rls()\' did not receive argument(s): ',
#              if (missing(trainingSet)) '\'trainingSet\' ' else '',
#              if (missing(trainingSetLabels)) '\'trainingSetLabels\' ' else '',
#              if (missing(testSet)) '\'testSet\' ' else '',
#              if (missing(numCompoundFeat)) '\'numCompoundFeat\' ' else '')
# }


#===============================================================================