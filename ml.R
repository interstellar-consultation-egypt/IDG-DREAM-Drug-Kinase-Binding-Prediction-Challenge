library(dimRed)
#===============================================================================

## Assumptions
## 'train' and 'test' have the same number of columns
## The columns consist of compound features followed by target features

#===============================================================================


rls_avg <- function(train,
                    trainLabels,
                    test,
                    compFeatIndx,
                    targFeatIndx,
                    predFunParams = list(lambda = 0.125),
                    dimReduction = 'none') {

    ## separate into compound and target features
    trainCompoundFeatures <- train[,compFeatIndx]
    trainTargetFeatures   <- train[,targFeatIndx]
    testCompoundFeatures <- test[,compFeatIndx]
    testTargetFeatures   <- test[,targFeatIndx]
    
    ## DIMENSIONALITY REDUCTION
    if (dimReduction != 'none') {
        trainCompoundFeatures <- dim_red(trainCompoundFeatures, dimReduction)
        trainTargetFeatures   <- dim_red(trainTargetFeatures,   dimReduction)
        testCompoundFeatures  <- dim_red(testCompoundFeatures,  dimReduction)
        testTargetFeatures    <- dim_red(testTargetFeatures,    dimReduction)
    }
    
    ## PARAMETER VALUES
    lambda <- predFunParams$lambda
    
    ## KERNELS
    Kc_train <- rbf_kernel(as.matrix(trainCompoundFeatures))
    Kt_train <- rbf_kernel(as.matrix(trainTargetFeatures))
    Kc_test  <- rbf_kernel(as.matrix(testCompoundFeatures),
                           as.matrix(trainCompoundFeatures))
    Kt_test  <- rbf_kernel(as.matrix(testTargetFeatures),
                           as.matrix(trainTargetFeatures))
    
    ## PREDICTIONS
          y <- trainLabels
    lambdaI <- lambda * diag(ncol(Kc_train))    ## lambda * identity matrix
     yhat_c <- Kc_test %*% (solve(Kc_train + lambdaI) %*% y)
     yhat_t <- Kt_test %*% (solve(Kt_train + lambdaI) %*% y)
       yhat <- (yhat_c + yhat_t) / 2
    
    ## return predicted labels
    yhat
}


#===============================================================================


rls_kron <- function(train,
                     trainLabels,
                     test,
                     compFeatIndx,
                     targFeatIndx,
                     predFunParams = list(lambda = 0.125),
                     dimReduction = 'none') {

    ## separate into compound and target features
    trainCompoundFeatures <- train[,compFeatIndx]
    trainTargetFeatures   <- train[,targFeatIndx]
    testCompoundFeatures <- test[,compFeatIndx]
    testTargetFeatures   <- test[,targFeatIndx]
    
    ## DIMENSIONALITY REDUCTION
    if (dimReduction != 'none') {
        trainCompoundFeatures <- dim_red(trainCompoundFeatures, dimReduction)
        trainTargetFeatures   <- dim_red(trainTargetFeatures,   dimReduction)
        testCompoundFeatures  <- dim_red(testCompoundFeatures,  dimReduction)
        testTargetFeatures    <- dim_red(testTargetFeatures,    dimReduction)
    }
    
    ## PARAMETER VALUES
    lambda <- predFunParams$lambda
    
    ## KERNELS
    Kc_train <- rbf_kernel(as.matrix(trainCompoundFeatures))
    Kt_train <- rbf_kernel(as.matrix(trainTargetFeatures))
    Kc_test  <- rbf_kernel(as.matrix(testCompoundFeatures),
                           as.matrix(trainCompoundFeatures))
    Kt_test  <- rbf_kernel(as.matrix(testTargetFeatures),
                           as.matrix(trainTargetFeatures))
    Ktrain   <- Kc_train * Kt_train
    Ktest    <- Kc_test * Kt_test
    
    ## PREDICTIONS
          y <- trainLabels
    lambdaI <- lambda * diag(ncol(Kc_train))    ## lambda * identity matrix
       yhat <- Ktest %*% (solve(Ktrain + lambdaI) %*% y)
    
    ## return predicted labels
    yhat
}


#===============================================================================


rls_kron_graphreg <- function(train,
                              trainLabels,
                              test,
                              compFeatIndx,
                              targFeatIndx,
                              predFunParams = list(lambda = 0.125, p = 5),
                              dimReduction = 'none') {

    ## separate into compound and target features
    trainCompoundFeatures <- train[,compFeatIndx]
    trainTargetFeatures   <- train[,targFeatIndx]
    testCompoundFeatures <- test[,compFeatIndx]
    testTargetFeatures   <- test[,targFeatIndx]
    
    ## DIMENSIONALITY REDUCTION
    if (dimReduction != 'none') {
        trainCompoundFeatures <- dim_red(trainCompoundFeatures, dimReduction)
        trainTargetFeatures   <- dim_red(trainTargetFeatures,   dimReduction)
        testCompoundFeatures  <- dim_red(testCompoundFeatures,  dimReduction)
        testTargetFeatures    <- dim_red(testTargetFeatures,    dimReduction)
    }
    
    ## PARAMETER VALUES
    lambda <- predFunParams$lambda
    p <- predFunParams$p
    
    ## KERNELS
    Kc_train <- rbf_kernel(as.matrix(trainCompoundFeatures))
    Kt_train <- rbf_kernel(as.matrix(trainTargetFeatures))
    Kc_test  <- rbf_kernel(as.matrix(testCompoundFeatures),
                           as.matrix(trainCompoundFeatures))
    Kt_test  <- rbf_kernel(as.matrix(testTargetFeatures),
                           as.matrix(trainTargetFeatures))
    Ktrain   <- Kc_train * Kt_train
    Ktest    <- Kc_test * Kt_test
    
    ## SPARSIFICATION (keep top 5 values in each row)
    if (p > 0) {
        Ktrain <- sparsifier(Ktrain, p)
        Ktest <- sparsifier(Ktest, p)
    }
    
    ## PREDICTIONS
          y <- trainLabels
    lambdaI <- lambda * diag(ncol(Kc_train))    ## lambda * identity matrix
       yhat <- Ktest %*% (solve(Ktrain + lambdaI) %*% y)
    
    ## return predicted labels
    yhat
}


#===============================================================================


ensembler <- function(train,
                      trainLabels,
                      test,
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
        train_i <- train
        trainLabels_i <- trainLabels
        test_i <- test
        
        ## BAGGING
        if (bag) {
            baggedInstances <- sample(nrow(train), replace = T)
            train_i <- train[baggedInstances,]
            trainLabels_i <- trainLabels[baggedInstances]
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
            train_i <- train_i[, selectedFeatures]
            test_i <- test_i[, selectedFeatures]
        }
        
        ## PREDICTIONS
        yhat <- yhat + baseLearner(train,
                                   trainLabels,
                                   test,
                                   1:numCompoundFeatures,
                                   (numCompoundFeatures+1):ncol(train),
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


dim_red <- function(dta, dimReduction) {
    
    ## For inspiration, type this in the console and hit Enter:
    ## 
    ##      dimRed::dimRedMethodList()
    ##
    ## The output is a list of dimensionality reduction techniques to consider
    
    
    ## DIMENSIONALITY REDUCTION
    if (dimReduction == 'pca') {
      dta <- prcomp(dta, retx = T, center = T)$x[,1:10]
        
    } else if (dimReduction == 'kpca') {
      dta <- kpca(dta, kernel = "rbfdot",
             kpar=list(sigma=0.2),features = 10 )
        
    } else if (dimReduction == 'isomap') {
      dta <- embed(dat, "Isomap", mute = NULL, knn = 10)
        
    } else if (dimReduction == 'lapeig') {
      leim <- LaplacianEigenmaps()
      emb <- leim@fun(dat, leim@stdpars, ndim = 10)
      dta <- emb@data@data
      
    } else if (dimReduction == 'mds') {
      mds <- MDS()
      emb <- mds@fun(dat, mds@stdpars, ndim = 10)
      
    } else if (dimReduction == 'mds') {
        ##...
        
    }
    
    
    ## return dimensionality-reduced data
    dta
}


#===============================================================================


sparsifier <- function(mtrx, p) {
    ## sparsify
    for (i in 1:nrow(mtrx)) {
        row_i <- mtrx[i,]
        sorted <- sort(row_i, decreasing = T, index.return = T)
        toBeKeptIndices <- sorted$ix[1:p]
        toBeZeroedIndices <- setdiff(1:ncol(mtrx), toBeKeptIndices)
        mtrx[i,toBeZeroedIndices] <- 0
    }
    
    ## return sparsified similarity matrix
    mtrx
}


#===============================================================================


# check_missing <- function(train, 
#                           trainLabels,
#                           test, 
#                           numCompoundFeat) {
#  
#     ## if any important arguments missing...   
#     if (missing(train) || 
#         missing(trainLabels) || 
#         missing(test) || 
#         missing(numCompoundFeat))
#         ## STOP!!!
#         stop('function \'ml.R::rls()\' did not receive argument(s): ',
#              if (missing(train)) '\'train\' ' else '',
#              if (missing(trainLabels)) '\'trainLabels\' ' else '',
#              if (missing(test)) '\'test\' ' else '',
#              if (missing(numCompoundFeat)) '\'numCompoundFeat\' ' else '')
# }


#===============================================================================