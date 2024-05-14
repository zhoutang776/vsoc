rm(list = ls())
if(Sys.info()[['sysname']] == "Darwin") {
    setwd("/Users/zhoutang/Onedrive/Research/VIM/code/vsoc")
}else {
    setwd("/home/zhoutang/vsoc")
}

library(foreach)
library(doParallel)
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


average.res <- NULL
for(j in seq_along(ranking.algorithms)) {
    # reset the random seed for each ranking algorithm
    set.seed(1234, kind = "L'Ecuyer-CMRG")
    
    rank.algorithm <- switch(ranking.algorithms[[j]],
                             "linear" = linear.ranking, "gam" = gam.ranking, "earth" = earth.ranking
    )
    get.truth.function <- switch(ranking.algorithms[[j]],
                                 "linear" = get.linear.truth, "gam" = get.gam.truth, "earth" = get.earth.truth
    )
    optimal.param <- get.truth.function()

    optimal.psi <- optimal.param$optimal.psi
    optimal.auc <- optimal.param$optimal.auc
    optimal.rank <- optimal.param$optimal.rank
    optimal.sel <- optimal.param$optimal.sel


    for(i in seq_along(sample.sizes)) {
        n <- sample.sizes[i]
        
        print(paste("ranking =", ranking.algorithms[[j]], ", n = ", n))
        
        start <- proc.time()[3]
        results <- foreach(1:iter, .inorder=FALSE, .combine="rbind") %dopar% {
            sim(n=n, alpha = alpha, V = V, bootstrap = F, B = NULL, 
                rank.algorithm = rank.algorithm,
                learners = learners)    
        }
        end <- proc.time()[3]
        print(sprintf("sample size n=%4d , elapse %4.2f seconds", n, (end-start)))
        
        save(results, file = sprintf("%s-%d-%d.RData", ranking.algorithms[[j]], n, V))

        est <- matrix(unlist(results[, 1]), nrow=6)
        std <- sqrt(matrix(unlist(results[, 2]), nrow=6)/n)
        
        bias <- est - c(optimal.psi, optimal.auc)
        sqrt.n.bias <- sqrt(n) * apply(bias, 1, mean)
        
        n.variance <- n * apply(std^2, 1, mean)
        
        lb <- est - 1.96 * std
        ub <- est + 1.96 * std
        
        cover <- (lb <= c(optimal.psi, optimal.auc)) & (c(optimal.psi, optimal.auc) <= ub)
        coverage <- apply(cover, 1, mean)
    
        rank_cover <- sapply(results[, 3], function(rank){
            cover_flag <- FALSE
            for(i in seq_along(optimal.rank)) {
                if(all(rank == optimal.rank[[i]])) {
                    cover_flag <- TRUE
                    break
                }
            }
            return(cover_flag)
        })
        rank_coverage <- mean(rank_cover)
        res <- data.frame(sqrt.n.bias = sqrt.n.bias,
                     n.variance = n.variance,
                     coverage = coverage,
                     rank_coverage = rank_coverage, 
                     n = n,
                     algo = ranking.algorithms[[j]])
        
        average.res <- rbind(average.res, res)
    }
}




average.res$algo[average.res$algo == "earth"] <- "MARS"
average.res$algo[average.res$algo == "linear"] <- "LASSO"
average.res$algo[average.res$algo == "gam"] <- "GAM"


# saveRDS(average.res, "average.res_April_27_cv.rds")
# saveRDS(average.res, "average.res_April_14_no_cv.rds")

# cv_average.res <- readRDS("average.res_April_27_cv.rds")
# nocv_average.res <- readRDS("average.res_April_27_no_cv.rds")

# cv_average.res <- readRDS("twisted_average.res_April_10.rds")
# nocv_average.res <- readRDS("average.res_April_14_no_cv.rds")

# ================================= VROC plot ====================================

cv_average.res["cv"] <- "Yes"
nocv_average.res["cv"] <- "No"
average.res <- rbind(nocv_average.res, cv_average.res)
average.res$cv <- factor(average.res$cv, levels=c("Yes", "No"))


library(ggplot2)
library(reshape2)
require(gridExtra)
library(ggpubr)
library(latex2exp)



fig1 <- ggplot(data=average.res[seq(6, nrow(average.res), 6), ]) +
    geom_line(aes(x=n, y=(sqrt.n.bias), colour=algo, linetype=cv)) +
    geom_point(aes(x=n, y=(sqrt.n.bias), colour=algo, shape=cv), size=2) +
    scale_shape_manual(name="Cross fitting", values=c(16, 5)) + 
    scale_linetype_manual(name="Cross fitting", values=c("solid", "dashed")) + 
    geom_hline(yintercept = 0, linetype="dashed") +
    coord_cartesian(ylim=c(-4, 4)) +
    scale_x_continuous(
        breaks=c(500, 1000, 2000, 3000, 4000, 5000), 
        label=c("500", "1K", "2K", "3K", "4K", "5K")
    ) +
    ylab(TeX("$n^{1/2} \\times$  bias of estimate")) + xlab("") +
    # labs(color='Ranking algorithm', linetype='Ranking algorithm') + 
    theme(legend.position = 'none')



fig2 <- ggplot(data=average.res[seq(6, nrow(average.res), 6), ]) +
    geom_line(aes(x=n, y=sqrt(n.variance), colour=algo, linetype=cv)) +
    geom_point(aes(x=n, y=sqrt(n.variance), colour=algo, shape=cv), size=2) +
    
    scale_shape_manual(name="Cross fitting", values=c(16, 5)) + 
    scale_linetype_manual(name="Cross fitting", values=c("solid", "dashed")) + 
    
    geom_hline(yintercept = 0, linetype="dashed") +
    coord_cartesian(ylim=c(0, 4)) +
    # scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(
        breaks=c(500, 1000, 2000, 3000, 4000, 5000), 
        label=c("500", "1K", "2K", "3K", "4K", "5K")
    ) +
    ylab(TeX("$n^{1/2} \\times$ estimated SD")) + xlab("Sample size") +
    # labs(color='Ranking algorithm', linetype='Ranking algorithm') + 
    theme(legend.position = 'none')

fig3 <- ggplot(data=average.res[seq(6, nrow(average.res), 6), ]) +
    geom_line(aes(x=n, y=coverage, colour=algo, linetype=cv)) +
    geom_point(aes(x=n, y=coverage, colour=algo, shape=cv), size=2) +
    
    scale_shape_manual(name="Cross fitting", values=c(16, 5)) + 
    scale_linetype_manual(name="Cross fitting", values=c("solid", "dashed")) + 
    scale_color_manual(name="Ranking algorithm", values=c("#F8766D", "#00BA38", "#619CFF")) + 

    geom_hline(yintercept = 0.95, linetype="dashed") +
    coord_cartesian(ylim=c(0.50, 1)) +
    scale_y_continuous(breaks=c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0),
                       label=c("0.50", "0.60", "0.70", "0.80", "0.90", "0.95", "1.00")
    ) +
    scale_x_continuous(
        breaks=c(500, 1000, 2000, 3000, 4000, 5000), 
        label=c("500", "1K", "2K", "3K", "4K", "5K")
    ) +
    ylab("Coverage of 95% CI") + xlab("") +
    # labs(color='Ranking algorithm', linetype='Ranking algorithm') +
    labs(color='Ranking algorithm') +
    theme(legend.key.width = unit(3, "line"), panel.spacing = unit(0.5, "lines"), legend.position="bottom",
          # strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12), 
          # text = element_text(size = 10)
    ) +
    guides(color=guide_legend(override.aes = list(size=1.8)), shape=guide_legend(override.aes = list(size=2.5)))

ggarrange(fig1, fig2, fig3, ncol=3, common.legend = TRUE, legend="bottom", widths = c(1, 1, 1.03),
          legend.grob=get_legend(fig3))
    



# AUVROC <- rep(FALSE, nrow(average.res))
# AUVROC[seq(6, 108, 6)] <- TRUE
# cond <- AUVROC & (average.res$algo == "GAM" | average.res$algo == "EARTH") & (average.res$n >= 3000)
# average.res[cond, ] 
# 
# average.res[cond, "coverage"] <- average.res[cond, "coverage"] +0.001




# -------------------------------- VROC plot -------------------------------- 

set.seed(1234, kind = "L'Ecuyer-CMRG")
linear.result <- vsoc(X = X, y = y, alpha = alpha, V = V, bootstrap = bootstrap, B = B,
                      rank.algorithm = linear.ranking, learners = learners)
set.seed(1234, kind = "L'Ecuyer-CMRG")
earth.result <- vsoc(X = X, y = y, alpha = alpha, V = V, bootstrap = bootstrap, B = B,
                     rank.algorithm = earth.ranking, learners = learners)

n <- 1000
true_LASSO <- c(0.1385120, 0.2255010, 0.6614436, 0.6981669, 0.7446607, 2.026698)
true_GAM <- c(0.3865451, 0.5497132, 0.6613689, 0.6980622, 0.7446607, 2.474747)
true_EARTH <- c(0.3865451, 0.5497132, 0.6613689, 0.6980622, 0.7446607, 2.474747)

LASSO <- c(0.1508122, 0.2292867, 0.6828857, 0.7088305, 0.7540599, 2.0674663)
GAM <- c(0.3981514, 0.5772118, 0.6828857, 0.7088305, 0.7540599, 2.5450337)
EARTH <- c(0.3981514, 0.5772118, 0.6828857, 0.7088305, 0.7540599, 2.5450337)

var_LASSO <- c(0.2675945, 0.4358583, 0.5411378, 0.4599158, 0.3499316, 3.2856025)
var_GAM <- c(0.7652492, 0.8720060, 0.5411378, 0.3690459, 0.3499316, 5.1474752)
var_EARTH <- c(0.7652492, 0.8720060, 0.5411378, 0.3690459, 0.3499316, 5.1474752)

lb_EARTH <- EARTH - 1.96 * sqrt(var_EARTH/n)
ub_EARTH <- EARTH + 1.96 * sqrt(var_EARTH/n)

lb_LASSO <- LASSO - 1.96 * sqrt(var_LASSO/n)
ub_LASSO <- LASSO + 1.96 * sqrt(var_LASSO/n)

lb_GAM <- GAM - 1.96 * sqrt(var_GAM/n)
ub_GAM <- GAM + 1.96 * sqrt(var_GAM/n)


uniform_lb_EARTH <- EARTH - 2.569214 * sqrt(var_EARTH/n)
uniform_ub_EARTH <- EARTH + 2.569214 * sqrt(var_EARTH/n)

uniform_lb_LASSO <- LASSO - 2.569214 * sqrt(var_LASSO/n)
uniform_ub_LASSO <- LASSO + 2.569214 * sqrt(var_LASSO/n)

uniform_lb_GAM <- GAM - 2.569214 * sqrt(var_GAM/n)
uniform_ub_GAM <- GAM + 2.569214 * sqrt(var_GAM/n)


lasso <- data.frame(true=true_LASSO[-6], psi=LASSO[-6], lb=lb_LASSO[-6], ub=ub_LASSO[-6], 
                    uniform_lb=uniform_lb_LASSO[-6], uniform_ub=uniform_ub_LASSO[-6], 
                    nvar=1:5, algo="LASSO")

gam <- data.frame(true=true_GAM[-6], psi=GAM[-6], lb=lb_GAM[-6], ub=ub_GAM[-6], 
                  uniform_lb=uniform_lb_GAM[-6], uniform_ub=uniform_ub_GAM[-6], 
                  nvar=1:5, algo="GAM")

earth <- data.frame(true=true_EARTH[-6], psi=EARTH[-6], lb=lb_EARTH[-6], ub=ub_EARTH[-6], 
                    uniform_lb=uniform_lb_EARTH[-6], uniform_ub=uniform_ub_EARTH[-6], 
                    nvar=1:5, algo="MARS")

vsoc <- data.frame(rbind(lasso, gam, earth))


ggplot(data=vsoc) +
    geom_line(aes(x=nvar, y=psi, colour=algo, linetype=algo)) +
    geom_line(aes(x=nvar, y=uniform_lb, colour=algo, linetype=algo)) +
    geom_line(aes(x=nvar, y=uniform_ub, colour=algo, linetype=algo)) +
    # scale_linetype_manual(values=c("dashed", "dotdash", "twodash")) + 
    geom_point(aes(x=nvar, y=true, shape=algo), size=1.5) +
    geom_point(aes(x=nvar, y=psi, colour=algo, shape=algo)) +
    geom_errorbar(mapping=aes(x=nvar, ymin=lb, ymax=ub, colour=algo, linetype=NULL), width=0.1) +
    coord_cartesian(ylim=c(0, 1)) + scale_y_continuous(limits=c(0, 1), expand = c(0, 0)) + 
    ylab("Estimated R-squared predictiveness measures") + xlab("Number of selected covariates") +
    labs(color='Ranking algorithm', linetype='Ranking algorithm', shape='Ranking algorithm') +
    theme_bw() + 
    theme(legend.key.width = unit(3, "line"), 
          # panel.spacing = unit(1, "lines"), 
          # legend.position="bottom",
          # strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12),
          # text = element_text(size = 10)
          ) 
  

    # + guides(color=guide_legend(nrow=2))





# -------------------------------- selection -------------------------------- 

# estimate plot

set.seed(12345, kind = "L'Ecuyer-CMRG")
n <- 1000
data <- generate.samples(n = n)
X <- scale(data$X)
y <- scale(data$y)


# GAM selection
model <- mgcv::bam(y ~ s(V1) + s(V2) + s(V3) + s(V4) + s(V5),
                   data = data.frame(y = y, X))
temp <- summary(model)
sort(temp$chi.sq, decreasing = TRUE, index.return = TRUE)$ix
sum(temp$s.table[, 4] < 0.05)
# LASSO selection
glmnet::cv.glmnet(x=as.matrix(X), y=y)
# MARS selection
model <- earth::earth(y ~ ., data = data.frame(y = y, X))
model


n <- 1000
true_LASSO <- c(0.1385120, 0.2255010, 0.6614436, 0.6981669, 0.7446607)
true_GAM <- c(0.3865451, 0.5497132, 0.6613689, 0.6980622, 0.7446607)
true_EARTH <- c(0.3865451, 0.5497132, 0.6613689, 0.6980622, 0.7446607)

LASSO <- c(0.1508122, 0.2292867, 0.6828857, 0.7088305, 0.7540599)
GAM <- c(0.3981514, 0.5772118, 0.6828857, 0.7088305, 0.7540599)
EARTH <- c(0.3981514, 0.5772118, 0.6828857, 0.7088305, 0.7540599)

var_LASSO <- c(0.2675945, 0.4358583, 0.5411378, 0.4599158, 0.3499316)
var_GAM <- c(0.7652492, 0.8720060, 0.5411378, 0.3690459, 0.3499316)
var_EARTH <- c(0.7652492, 0.8720060, 0.5411378, 0.3690459, 0.3499316)

lb_EARTH <- EARTH - 1.96 * sqrt(var_EARTH/n)
ub_EARTH <- EARTH + 1.96 * sqrt(var_EARTH/n)

lb_LASSO <- LASSO - 1.96 * sqrt(var_LASSO/n)
ub_LASSO <- LASSO + 1.96 * sqrt(var_LASSO/n)

lb_GAM <- GAM - 1.96 * sqrt(var_GAM/n)
ub_GAM <- GAM + 1.96 * sqrt(var_GAM/n)


uniform_lb_EARTH <- EARTH - 2.569214 * sqrt(var_EARTH/n)
uniform_ub_EARTH <- EARTH + 2.569214 * sqrt(var_EARTH/n)

uniform_lb_LASSO <- LASSO - 2.569214 * sqrt(var_LASSO/n)
uniform_ub_LASSO <- LASSO + 2.569214 * sqrt(var_LASSO/n)

uniform_lb_GAM <- GAM - 2.569214 * sqrt(var_GAM/n)
uniform_ub_GAM <- GAM + 2.569214 * sqrt(var_GAM/n)


lasso <- data.frame(true=true_LASSO, psi=LASSO, lb=lb_LASSO, ub=ub_LASSO, 
                    uniform_lb=uniform_lb_LASSO, uniform_ub=uniform_ub_LASSO, 
                    nvar=1:5, algo="LASSO")

gam <- data.frame(true=true_GAM, psi=GAM, lb=lb_GAM, ub=ub_GAM, 
                  uniform_lb=uniform_lb_GAM, uniform_ub=uniform_ub_GAM, 
                  nvar=1:5, algo="GAM")

earth <- data.frame(true=true_EARTH, psi=EARTH, lb=lb_EARTH, ub=ub_EARTH, 
                    uniform_lb=uniform_lb_EARTH, uniform_ub=uniform_ub_EARTH, 
                    nvar=1:5, algo="MARS")



lasso.min <- lasso[3, ]
lasso.1se <- lasso[1, ]
original <- c(NA,0,NA,NA,NA,NA,0,NA)

lasso.min <- rbind(lasso.min, original)
lasso.1se <- rbind(lasso.1se, original)

# lasso.min$algo <- "LASSO w/ min"
# lasso.1se$algo <- "LASSO w/ 1se"
lasso.min$algo <- "LASSO"

gam.pvalue <- gam[5, ]
gam.pvalue <- rbind(gam.pvalue, original)
gam.pvalue$algo <- "GAM"

earth.vim <- earth[3, ]
earth.vim <- rbind(earth.vim, original)
earth.vim$algo <- "MARS"

selected.vsoc <- rbind(lasso.min, gam.pvalue, earth.vim)

ggplot(data=selected.vsoc) +
    geom_line(aes(x=nvar, y=psi, colour=algo, linetype=algo)) +
    # scale_color_manual(name="Selection algorithm", values=c("#F8766D", "#00BA38", "#619CFF"),
    #                    labels = unname(
    #                        TeX(c("GAM", "LASSO", "MARS"))
    #                        )) +
    geom_point(aes(x=nvar, y=psi, colour=algo)) +
    geom_point(aes(x=nvar, y=true)) +
    geom_errorbar(mapping=aes(x=nvar, ymin=lb, ymax=ub, colour=algo, linetype=NULL), width=0.1) +
    coord_cartesian(ylim=c(0, 1)) +
    ylab("Estimated R-squared predictiveness measures") + xlab("Number of selected covariates") +
    labs(color='Selection algorithm', linetype='Selection algorithm') +
    theme(legend.key.width = unit(2, "line"), panel.spacing = unit(0.5, "lines"), legend.position="bottom",
    )











