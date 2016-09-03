
################################################################################## Auxiliar functions

## Calculate the square root of a matrix
## Input:
## A: a matrix
## Output: the square root of matrix A

sqrtMatrix <- function(A) {
    V <- eigen(A)$vectors
    Dtemp <- qr.solve(V) %*% A %*% V
    D <- matrix(0, ncol(A), ncol(A))
    diag(D) <- sqrt(diag(Dtemp))
    return(V%*%D%*%qr.solve(V))
}

## Block matrix determinant

detMatrix <- function(Sigma) {

    A <- Sigma[[1]]
    B <- Sigma[[2]]
    C <- Sigma[[3]]
    D <- Sigma[[4]]

    res <- det(A) * det(D - C %*% qr.solve(A) %*% B)
    return(res)
}

################################################################################## Performs the Partial Canonical Correlation Analysis (PCCA)
## Input:
## X: a matrix containing the time series in the columns (predictors)
## Y: a matrix containing the time series in the columns (response variables)
## groups: an array identifying which group each gene in X belongs
## Output:
## B: an array containing the canonical correlation coefficients for each
## predictor
## boot: boolean (if FALSE, likelihood ratio test is performed)
## pvalue: matrix containing the p-values
## DF: degrees of freedom (of the chi-square test)
## est: statistic used in the test

gGranger <- function(X, Y, groups, boot=FALSE) {

    B <- matrix(0, max(groups), 1)
    pvalue <- matrix(0, max(groups),1)


    for (gr in 1:max(groups)) {
        in_ <- which(groups==gr)
        out_ <- which(groups != gr)

        SigmaXX <- cov(as.matrix(X[,in_]))
        SigmaYY <- cov(as.matrix(Y))
        SigmaXY <- cov(as.matrix(X[,in_]),as.matrix(Y))
        SigmaYX <- t(SigmaXY)

        if (length(out_ > 0)) {
            SigmaXZ <- cov(as.matrix(X[,in_]),
                               as.matrix(X[,out_]))
            SigmaZX <- t(SigmaXZ)
            SigmaZZ <- cov(as.matrix(X[,out_]))
            SigmaZY <- cov(as.matrix(X[,out_]),as.matrix(Y))
            SigmaYZ <- t(SigmaZY)

            invSigmaZZ <- qr.solve(SigmaZZ)

            SigmaXX.Z <- SigmaXX - SigmaXZ %*% invSigmaZZ %*% SigmaZX
            SigmaXY.Z <- SigmaXY - SigmaXZ %*% invSigmaZZ %*% SigmaZY
            SigmaYX.Z <- SigmaYX - SigmaYZ %*% invSigmaZZ %*% SigmaZX
            SigmaYY.Z <- SigmaYY - SigmaYZ %*% invSigmaZZ %*% SigmaZY
        }
        else {
            SigmaXX.Z <- SigmaXX
            SigmaXY.Z <- SigmaXY
            SigmaYX.Z <- SigmaYX
            SigmaYY.Z <- SigmaYY
        }

        Sigma <- list()
        Sigma[[1]] <- SigmaXX.Z
        Sigma[[2]] <- SigmaXY.Z
        Sigma[[3]] <- SigmaYX.Z
        Sigma[[4]] <- SigmaYY.Z

        B[gr,1] <- sqrt(eigen(qr.solve(sqrtMatrix(SigmaXX.Z)) %*%
                   SigmaXY.Z %*% qr.solve(SigmaYY.Z) %*% SigmaYX.Z %*%
                   qr.solve(sqrtMatrix(SigmaXX.Z)))$values[1])

        if (boot == FALSE) {
            ## Likelihood ratio test
            est <- (nrow(X)-length(out_)-1-0.5*
                   (length(in_)+ncol(Y)+1))*
                   log((det(SigmaXX.Z)*det(SigmaYY.Z))/detMatrix(Sigma))

            ## Degrees of freedom (p X q)
            DF <- length(in_) * ncol(Y)
            pvalue[gr,1] <- 1 - pchisq(est, DF)
        }
    }

    res <- list()
    res$B <- B
    if (boot == FALSE) {
        res$pvalue <- pvalue 
    }
    return(res)
}

################################################################################## Identify the existence of Granger causality from a set of time-series to 
## another set of time-series.
## Block bootstrap
## Input:
## Z: contains the time series in the columns
## groups: an array containing the labels of the groups
## bootMax: number of bootstrap samples (default=1000)
## tamBlock: length of the block (default=10)
## Output:
## pvalue: matrix containing the p-values
## beta: matrix containing the coefficients obtained using Partial CCA.

gGranger.like <- function(Z, groups) {

    numSample <- nrow(Z)
    numGenesX <- length(groups)

    pvalue <- matrix(0, max(groups), max(groups))
    beta <- matrix(0, max(groups), max(groups))

    ## Normalization for mean of zero and variance of one
    for (i in 1:numGenesX) {
        Z[,i] <- (Z[,i]-mean(Z[,i]))/sd(Z[,i])
    }

    for (gr in 1:max(groups)) {

        numGenesY <- length(which(groups==gr))

        corT <- gGranger(as.matrix(Z[1:(numSample-1),]), as.matrix(Z[2:numSample,which(groups==gr)]), groups, boot=FALSE)

        for (i in 1:max(groups)) {
            pvalue[i,gr] <- corT$pvalue[i,1]
            beta[i,gr] <- corT$B[i,1]
        }
    }

    res <- list()
    res$beta <- beta
    res$pvalue <- pvalue

    return(res)
}



################################################################################## Identify the existence of Granger causality from a set of time-series to 
## another set of time-series.
## Block bootstrap
## Input:
## Z: contains the gene expressions in the columns
## groups: an array containing the labels of the groups
## bootMax: number of bootstrap samples (default=1000)
## tamBlock: length of the block (default=10)
## Output:
## pvalue: a matrix containing the p-values
## beta: matrix containing the coefficients obtained using Partial CCA.

gGranger.boot <- function(Z, groups, tamBlock=10, bootMax=1000) {

    numSample <- nrow(Z)
    numGenesX <- length(groups)

    pvalue <- matrix(0, max(groups), max(groups))
    beta <- matrix(0, max(groups), max(groups))

    ## Normalization for mean of zero and variance of one
    for (i in 1:numGenesX) {
        Z[,i] <- (Z[,i]-mean(Z[,i]))/sd(Z[,i])
    }

    for (gr in 1:max(groups)) {

        numGenesY <- length(which(groups==gr))

        corT <- gGranger(Z[1:(numSample-1),], Z[2:numSample,which(groups==gr)], groups, boot=TRUE)$B
        distr <- matrix(0, max(groups), bootMax)

        X <- Z
        Y <- matrix(Z[,which(groups==gr)], numSample, length(which(groups==gr)))

        ## Block bootstrap
        for (boot in 1:bootMax) {

            blksX <- sample(seq(1:(numSample-tamBlock+1)),((numSample/tamBlock)+1), replace=TRUE)
            blksY <- sample(seq(1:(numSample-tamBlock+1)),((numSample/tamBlock)+1), replace=TRUE)

            X_ <- matrix(0, ((numSample/tamBlock)+1)*tamBlock, numGenesX)
            Y_ <- matrix(0, ((numSample/tamBlock)+1)*tamBlock, numGenesY)

            for (i in 1:((numSample/tamBlock)+1)) {
                X_[((i-1)*tamBlock+1):(i*tamBlock),] <- 
                                             X[blksX[i]:(blksX[i]+tamBlock-1),]
                Y_[((i-1)*tamBlock+1):(i*tamBlock),] <- 
                                             Y[blksY[i]:(blksY[i]+tamBlock-1),]
            }

            Xtemp <- matrix(0, numSample, numGenesX)
            Ytemp <- matrix(0, numSample, numGenesY)

            Xtemp <- as.matrix(X_[1:numSample,])
            Ytemp <- as.matrix(Y_[1:numSample,])

     
            distr[,boot] <- abs(gGranger(matrix(Xtemp[1:(numSample-1),],
                               (numSample-1), ncol(Xtemp)),
                                matrix(Ytemp[2:numSample,], (numSample-1),
                                ncol(Ytemp)), groups, boot=TRUE)$B)
        }

        for (i in 1:max(groups)) {
            pvalue[i,gr] <- 1-rank(c(abs(corT[i,1]),distr[i,]),
                            ties="first")[1]/(bootMax+1)
            beta[i,gr] <- corT[i,1]
        }
    }

    res <- list()
    res$beta <- beta
    res$pvalue <- pvalue

    return(res)
}

###############################################################################

