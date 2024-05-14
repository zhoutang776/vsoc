rm(list = ls())
library(SuperLearner)

registerDoParallel(cores=parallel::detectCores())
set.seed(1234, kind = "L'Ecuyer-CMRG")

source("vsoc.R")

# =========================== Parameter Definition =============================
SL.gam.customize <- function(...){
    SuperLearner::SL.gam(..., deg.gam = 5)
}
SL.xgboost.customize <- function(...) {
    SL.xgboost(
        ...,
        shrinkage = 0.3,
        max_depth = 10,
        ntrees = 100,
        minobspernode = 5,
        interaction_constraints = list(c(0, 1))
    )
}

learners <- c("SL.mean", "SL.loess", "SL.xgboost", "SL.gam.customize", "SL.earth")
learners <- c("SL.mean", "SL.loess", "SL.gam.customize", "SL.earth", "")
learners <- c("SL.mean", "SL.polymars", "SL.gam.customize", "SL.earth")
learners <- c("SL.xgboost.customize", "SL.gam.customize", "SL.earth")

sample.sizes <- c(500, 1000, 2000, 3000, 4000, 5000)
ranking.algorithms <- list("linear", "gam", "earth")
ranking.algorithms <- list("gam")

alpha <- 0.05
V <- 10
seed <- 123
iter <- 500

# =========================== Function Definition ==============================
generate.samples <- function(n){
    # generate sample using Y = (X1+0.5) * (X2+1) + sqrt(max(X2+5, 0)) + 5*sqrt((X3-0.2)^2+1) + epsilon
    # and epsilon ~ normal(0, f(X4, X5))
    # In this case, it contains: 
    #   1) null variable
    #   2) heteroskedasticity
    #   3) interaction term
    #   4) non-linear term 
    #   5) correlation between variables
    
    p <- 5
    mu <- rep(0, p)  # mean vector
    rho <- 0.4 #  control independent or not
    Sigma <- matrix(rho, ncol = p, nrow = p)  # covariance matrix
    diag(Sigma) <- 1
    X <- mvtnorm::rmvnorm(n, mean = mu, sigma = Sigma)
    
    mean.function <- function(X) (X[, 1] + 0.5) * (X[, 2] + 1) +
        sqrt(pmax(X[, 2] + 5, 0)) +
        5 * sqrt((X[, 3] - 0.2)^2 + 1)
    
    std.function <- function(X) sqrt( 1 + abs(X[, 4]) + abs(X[, 5]) )
    
    epsilon <- rnorm(n, 0, std.function(X))
    y <- mean.function(X) + epsilon
    X <- as.data.frame(X)
    y <- as.matrix(y, ncol = 1)
    return(list(X = X, y = y, epsilon = epsilon, mu = mu, Sigma = Sigma))
}


conditional_normal <- function(index_1, index_2, X, mu, Sigma,
                               keepdim = TRUE){
    # assume X ~ multivariate normal distribution
    # return the conditional mean/covariance matrix given X_(index_2)
    
    # expect to receive
    # X a n-by-(p+q) matrix
    # mu a (p+q)-length vector
    # Sigma a (p+q)-by-(p+q) matrix
    X <- as.matrix(X)
    X_1 <- as.matrix(X[, index_1, drop = FALSE])
    X_2 <- as.matrix(X[, index_2, drop = FALSE])
    index <- c(index_1, index_2)
    mu <- as.vector(mu[index])
    Sigma <- as.matrix(Sigma[index, index])
    
    n <- nrow(X_1)
    p <- ncol(X_1); q <- ncol(X_2)
    # sanity check!
    valid <- TRUE
    if(nrow(X_1) != nrow(X_2)) valid <- FALSE
    if(p + q != length(mu)) valid <- FALSE
    if(!isSymmetric.matrix(Sigma)) valid <- FALSE
    if(!valid) {
        stop("wrong parameters")
    }
    
    mu_1 <- mu[1:p]
    mu_2 <- mu[(p + 1):(p + q)]
    Sigma_11 <- matrix(Sigma[1:p, 1:p], nrow = p, ncol = p)
    Sigma_12 <- matrix(Sigma[1:p, (p + 1):(p + q)], nrow = p, ncol = q)
    Sigma_21 <- matrix(t(Sigma_12), nrow = q, ncol = p)
    Sigma_22 <- matrix(Sigma[(p + 1):(p + q), (p + 1):(p + q)], nrow = q, ncol = q)
    
    mean <- t(mu_1 + Sigma_12 %*% solve(Sigma_22) %*% (t(X_2) - mu_2))
    cov <- Sigma_11 - Sigma_12 %*% solve(Sigma_22) %*% Sigma_21
    
    if(keepdim) {
        mu <- matrix(NA, ncol = (p + q), nrow = n)
        mu[, index_1] <- mean
        Sigma <- matrix(NA, ncol = (p + q), nrow = (p + q))
        Sigma[index_1, index_1] <- cov
        mean <- mu
        cov <- Sigma
    }
    return(list(mean = mean, cov = cov))
}

# =========================== Population Parameters ==============================
# linear rank
linear.ranking <- function(X, y){
    lars.fit <- lars::lars(as.matrix(X), y, type = "lar")
    rank <- unlist(lars.fit$actions)
    stopifnot(length(rank) == ncol(X))
    
    glmnet.fit <- glmnet::cv.glmnet(as.matrix(X), y)
    lambda.min.index <- glmnet.fit$index[1]
    lambda.1se.index <- glmnet.fit$index[2]
    select <- c(glmnet.fit$nzero[lambda.min.index], glmnet.fit$nzero[lambda.1se.index])
    
    return(list(rank=as.vector(rank), select=as.numeric(select)))
}
get.linear.truth <- function(n = 1e5, m = 1e4, seed = 123){
    return(
        list(optimal.rank=list(c(1, 2, 3, 4, 5), c(1, 2, 3, 5, 4)),
             optimal.psi=c(0.1385120, 0.2255010, 0.6614436, 0.6981669, 0.7446607),
             optimal.auc=2.026698,
             optimal.sel=0.2204812)
    )
    
    # data <- generate.samples(5e6)
    # rank <- linear.ranking(data$X, data$y)
    
    rank <- c(1, 2, 3, 4, 5)
    data <- generate.samples(n)
    X <- data$X
    y <- data$y
    epsilon <- data$epsilon
    mu <- data$mu
    Sigma <- data$Sigma
    p <- ncol(X)
    
    
    var_y <- mean((y - mean(y))^2)
    mse_full <- mean((epsilon - mean(epsilon))^2)
    
    # first one
    temp <- conditional_normal(rank[1], rank[2:p], X, mu, Sigma)
    mean <- temp$mean; cov <- temp$cov
    Ey_given_X2345 <- (mean[, rank[1]] + 0.5) * (X[, 2] + 1) +
        sqrt(pmax(X[, 2] + 5, 0)) +
        5 * sqrt((X[, 3] - 0.2)^2 + 1)
    mse_redu_1 <- mean((y - Ey_given_X2345)^2)
    
    # second one
    temp <- conditional_normal(rank[1:2], rank[3:p], X, mu, Sigma)
    mean <- temp$mean; cov <- temp$cov
    Ey_given_X345 <- matrix(0, nrow = n)
    for(i in 1:n) {
        data.temp <- matrix(NA, nrow = m, ncol = p)
        index <- rank[1:2]
        temp <- mvtnorm::rmvnorm(m, mean = mean[i, index], sigma = as.matrix(cov[index, index]))
        data.temp[, index] <- temp
        
        Ey_given_X345[i] <- mean((data.temp[, 1] + 0.5) * (data.temp[, 2] + 1) +
                                     sqrt(pmax(data.temp[, 2] + 5, 0))) +
            5 * sqrt((X[i, 3] - 0.2)^2 + 1)
        
    }
    mse_redu_12 <- mean((y - Ey_given_X345)^2)
    
    # third one
    temp <- conditional_normal(rank[1:3], rank[4:p], X, mu, Sigma)
    mean <- temp$mean; cov <- temp$cov
    Ey_given_X45 <- matrix(0, nrow = n)
    for(i in 1:n) {
        data.temp <- matrix(NA, nrow = m, ncol = p)
        index <- rank[1:3]
        temp <- mvtnorm::rmvnorm(m, mean = mean[i, index], sigma = as.matrix(cov[index, index]))
        data.temp[, index] <- temp
        
        Ey_given_X45[i] <- mean(
            (data.temp[, 1] + 0.5) * (data.temp[, 2] + 1) +
                sqrt(pmax(data.temp[, 2] + 5, 0)) +
                5 * sqrt((data.temp[, 3] - 0.2)^2 + 1)
        )
    }
    mse_redu_123 <- mean((y - Ey_given_X45)^2)
    
    # forth one
    temp <- conditional_normal(rank[1:4], rank[p], X, mu, Sigma)
    mean <- temp$mean; cov <- temp$cov
    Ey_given_X5 <- matrix(0, nrow = n)
    for(i in 1:n) {
        data.temp <- matrix(NA, nrow = m, ncol = p)
        index <- rank[1:4]
        temp <- mvtnorm::rmvnorm(m, mean = mean[i, index], sigma = as.matrix(cov[index, index]))
        data.temp[, index] <- temp
        
        Ey_given_X5[i] <- mean(
            (data.temp[, 1] + 0.5) +
                (data.temp[, 2] + 1) +
                sqrt(pmax(data.temp[, 2] + 5, 0)) +
                5 * sqrt((data.temp[, 3] - 0.2)^2 + 1)
        )
    }
    mse_redu_1234 <- mean((y - Ey_given_X5)^2)
    
    # fifth one
    mse_redu_12345 <- var_y
    
    
    mse_redu <- c(mse_redu_1, mse_redu_12, mse_redu_123, mse_redu_1234, mse_redu_12345)
    
    theta <- mse_redu - mse_full
    
    psi <- theta / var_y
    
    auc <- sum(psi) - (psi[1] + psi[p]) / 2
    return(
        list(optimal.rank = list(c(1, 2, 3, 4, 5), c(1, 2, 3, 5, 4)),
             optimal.mse_full = mse_full,
             optimal.mse_redu = mse_redu,
             optimal.theta = theta,
             optimal.psi = psi,
             optimal.auc = auc)
    )
}

# gam rank
gam.ranking <- function(X, y){
    model <- mgcv::bam(y ~ s(V1) + s(V2) + s(V3) + s(V4) + s(V5),
                       data = data.frame(y = y, X))
    temp <- summary(model)
    rank <- sort(temp$chi.sq, decreasing = TRUE, index.return = TRUE)$ix
    select <- sum(temp$s.table[, 4] < 0.05)
    return(list(rank=rank, select=select))
}
get.gam.truth <- function(n = 1e5, m = 1e4, seed = 123){
    return(
        list(optimal.rank = list(c(3, 1, 2, 4, 5), c(3, 1, 2, 5, 4)),
             optimal.psi = c(0.3865451, 0.5497132, 0.6613689, 0.6980622, 0.7446607),
             optimal.auc = 2.474747,
             optimal.sel = 0.1396124)
    )
    
    # data <- generate.samples(1e5, seed=seed)
    # rank <- gam.ranking(data$X, data$y)
    
    rank <- c(3, 1, 2, 4, 5)
    data <- generate.samples(n, seed = seed)
    X <- data$X
    y <- data$y
    epsilon <- data$epsilon
    mu <- data$mu
    Sigma <- data$Sigma
    p <- ncol(X)
    
    # without any covariates
    var_y <- var(y)
    predictiveness_full <- 1 - var(epsilon) / var_y
    
    # first one
    temp <- conditional_normal(rank[1], rank[2:p], X, mu, Sigma)
    mean <- temp$mean; cov <- temp$cov
    Ey_given_X1245 <- matrix(0, nrow = n)
    for(i in 1:n) {
        data.temp <- matrix(NA, nrow = m, ncol = p)
        index <- rank[1]
        temp <- mvtnorm::rmvnorm(m, mean = mean[i, index], sigma = as.matrix(cov[index, index]))
        data.temp[, index] <- temp
        
        Ey_given_X1245[i] <- mean(
            (X[i, 1] + 0.5) * (X[i, 2] + 1) +
                sqrt(pmax(X[i, 2] + 5, 0)) +
                5 * sqrt((data.temp[, 3] - 0.2)^2 + 1)
        )
    }
    predictiveness_reduce_3 <- 1 - mean((y - Ey_given_X1245)^2) / var_y
    
    # second one
    temp <- conditional_normal(rank[1:2], rank[3:p], X, mu, Sigma)
    mean <- temp$mean; cov <- temp$cov
    Ey_given_X245 <- matrix(0, nrow = n)
    for(i in 1:n) {
        data.temp <- matrix(NA, nrow = m, ncol = p)
        index <- rank[1:2]
        temp <- mvtnorm::rmvnorm(m, mean = mean[i, index], sigma = as.matrix(cov[index, index]))
        data.temp[, index] <- temp
        
        Ey_given_X245[i] <- mean(
            (data.temp[, 1] + 0.5) * (X[i, 2] + 1) +
                sqrt(pmax(X[i, 2] + 5, 0)) +
                5 * sqrt((data.temp[, 3] - 0.2)^2 + 1)
        )
    }
    predictiveness_reduce_31 <- 1 - mean((y - Ey_given_X245)^2) / var_y
    
    # third one
    temp <- conditional_normal(rank[1:3], rank[4:p], X, mu, Sigma)
    mean <- temp$mean; cov <- temp$cov
    Ey_given_X45 <- matrix(0, nrow = n)
    for(i in 1:n) {
        data.temp <- matrix(NA, nrow = m, ncol = p)
        index <- rank[1:3]
        temp <- mvtnorm::rmvnorm(m, mean = mean[i, index], sigma = as.matrix(cov[index, index]))
        data.temp[, index] <- temp
        
        Ey_given_X45[i] <- mean(
            (data.temp[, 1] + 0.5) * (data.temp[, 2] + 1) +
                sqrt(pmax(data.temp[, 2] + 5, 0)) +
                5 * sqrt((data.temp[, 3] - 0.2)^2 + 1)
        )
    }
    predictiveness_reduce_312 <- 1 - mean((y - Ey_given_X45)^2) / var_y
    
    # forth one
    temp <- conditional_normal(rank[1:4], rank[p], X, mu, Sigma)
    mean <- temp$mean; cov <- temp$cov
    Ey_given_X5 <- matrix(0, nrow = n)
    for(i in 1:n) {
        data.temp <- matrix(NA, nrow = m, ncol = p)
        index <- rank[1:4]
        temp <- mvtnorm::rmvnorm(m, mean = mean[i, index], sigma = as.matrix(cov[index, index]))
        data.temp[, index] <- temp
        
        Ey_given_X5[i] <- mean(
            (data.temp[, 1] + 0.5) * (data.temp[, 2] + 1) +
                sqrt(pmax(data.temp[, 2] + 5, 0)) +
                5 * sqrt((data.temp[, 3] - 0.2)^2 + 1)
        )
    }
    predictiveness_reduce_3124 <- 1 - mean((y - Ey_given_X5)^2) / var_y
    
    # fifth one
    predictiveness_reduce_31245 <- 1 - var_y / var_y
    
    psi <- vector()
    # psi for the X1
    psi[1] <- predictiveness_full - predictiveness_reduce_3
    # psi for the X1, X2
    psi[2] <- predictiveness_full - predictiveness_reduce_31
    # psi for the X1, X2, X3
    psi[3] <- predictiveness_full - predictiveness_reduce_312
    # psi for the X1, X2, X3, X4
    psi[4] <- predictiveness_full - predictiveness_reduce_3124
    # psi for the X1, X2, X3, X4, X5
    psi[5] <- predictiveness_full - predictiveness_reduce_31245
    
    auc <- sum(psi) - (psi[1] + psi[p]) / 2
    return(
        list(optimal.rank = list(c(3, 1, 2, 4, 5), c(3, 1, 2, 5, 4)),
             optimal.psi = psi,
             optimal.auc = auc)
    )
}

# earth rank
earth.ranking <- function(X, y){
    model <- earth::earth(y ~ ., data = data.frame(y = y, X))
    str_rank <- row.names(earth::evimp(model, trim = FALSE))
    rank <- as.numeric(gsub("V([0-9]+).*$", "\\1", str_rank))
    select <- sum(!endsWith(str_rank, "unused"))
    return(list(rank=rank, select=select))
}
get.earth.truth <- function(n = 1e4, m = 1e4, seed = 123){
    return(
        list(optimal.rank = list(c(3, 1, 2, 4, 5), c(3, 1, 2, 5, 4)),
             optimal.psi = c(0.3865451, 0.5497132, 0.6613689, 0.6980622, 0.7446607),
             optimal.auc = 2.474747,
             optimal.sel = 0.2204563)
    )
    # data <- generate.samples(1e6, seed=seed)
    # rank <- earth.ranking(data$X, data$y)
}



sim <- function(n = 1000, alpha = 0.05, V = 5, bootstrap = FALSE, B = 500, 
                rank.algorithm = function(X, y){ seq_len(ncol(X)) },
                learners = c("SL.mean", "SL.lm")){
    data <- generate.samples(n = n)
    X <- as.data.frame(scale(data$X))
    y <- scale(data$y)
    
    rank.select <- rank.algorithm(X, y)
    result <- vsoc(X = X, y = y, alpha = alpha, V = V, bootstrap = F, B = NULL,
                   rank = rank.select$rank, learners = learners)
    
    
    if(is.null(result)) {
        return(list(X = X, y = y))
    } else {
        psi_est <- result$psi_est; psi_var <- result$psi_var;
        auc_est <- result$auc_est; auc_var <- result$auc_var;
        
        est <- c(psi_est, auc_est)
        var <- c(psi_var, auc_var)
    }
    return(list(est = est, var = var, rank = result$rank, select = rank.select$select))
}


# ================================= Running ====================================
# reset the random seed for each ranking algorithm
set.seed(1234, kind = "L'Ecuyer-CMRG")
linear_res <- sim(n=n, alpha = alpha, V = V, bootstrap = F, B = NULL, 
               rank.algorithm = linear.ranking,
               learners = learners) 

set.seed(1234, kind = "L'Ecuyer-CMRG")
gam_res <- sim(n=n, alpha = alpha, V = V, bootstrap = F, B = NULL, 
                  rank.algorithm = gam.ranking,
                  learners = learners) 

set.seed(1234, kind = "L'Ecuyer-CMRG")
earth_res <- sim(n=n, alpha = alpha, V = V, bootstrap = F, B = NULL, 
                  rank.algorithm = earth.ranking,
                  learners = learners) 














