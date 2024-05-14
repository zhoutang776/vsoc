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
        ntrees = 100,
        # minobspernode = 5,
    )
}

learners <- c("SL.mean", "SL.loess", "SL.xgboost", "SL.gam.customize", "SL.earth")
learners <- c("SL.mean", "SL.loess", "SL.gam.customize", "SL.earth", "")
learners <- c("SL.mean", "SL.polymars", "SL.gam.customize", "SL.earth")
learners <- c("SL.xgboost.customize", "SL.gam.customize", "SL.earth")

sample.sizes <- c(500, 1000, 2000, 3000, 4000, 5000)
ranking.algorithms <- list("linear", "gam", "earth")

sample.sizes <- c(10000, 20000)
ranking.algorithms <- list("gam")

alpha <- 0.05
V <- 5
seed <- 123
iter <- 200
bootstrap <- T
B <- 500

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
             optimal.auc=2.026698)
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
             optimal.auc = 2.474747)
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
             optimal.auc = 2.474747)
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

    result <- vsoc(X = X, y = y, alpha = alpha, V = V, bootstrap = bootstrap, B = B,
                   rank.algorithm = rank.algorithm, learners = learners)

    if(is.null(result)) {
        return(list(X = X, y = y))
    } else {
        psi_est <- result$psi_est; psi_var <- result$psi_var;
        auc_est <- result$auc_est; auc_var <- result$auc_var;
        auc_se <- sqrt(auc_var)
        
        est <- c(psi_est, auc_est)
        var <- c(psi_var, auc_var)
        cis <- NULL
        
        if(bootstrap){
            boot_auc_est <- result$bootstrap_auc_est
            boot_auc_var <- result$bootstrap_auc_var
            boot_auc_se <- sqrt(boot_auc_var)
            
            cis <- get_boot_ci(bs_est=boot_auc_est, bs_est0=auc_est, bs_se=auc_se, est=auc_est, se=auc_se)
            cis$wald_ci <- result$auc_ci
        }
    }
    return(list(est = est, var = var, rank = result$rank, cis=cis, result=result))
}


# ================================= Running ====================================


average.res <- NULL
final.boot_coverage <- NULL
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
    
    
    for(i in seq_along(sample.sizes)) {
        n <- sample.sizes[i]
        
        start <- proc.time()[3]
        results <- foreach(1:iter, .inorder=FALSE, .combine="rbind") %dopar% {
            sim(n=n, alpha = alpha, V = V, bootstrap = bootstrap, B = B, 
                rank.algorithm = rank.algorithm,
                learners = learners)    
        }
        end <- proc.time()[3]
        print(sprintf("sample size n=%4d , elapse %4.2f seconds", n, (end-start)))
        
        save(results, file = sprintf("boot-%s-%d.RData", ranking.algorithms[[j]], n))
        
        est <- matrix(unlist(results[, 1]), nrow=6)
        std <- sqrt(matrix(unlist(results[, 2]), nrow=6)/n)
        
        bias <- est - c(optimal.psi, optimal.auc)
        sqrt.n.bias <- sqrt(n) * apply(bias, 1, mean)
        
        n.variance <- n * apply(std^2, 1, mean)
        
        lb <- est - 1.96 * std
        ub <- est + 1.96 * std
        
        cover <- (lb <= c(optimal.psi, optimal.auc)) & (c(optimal.psi, optimal.auc) <= ub)
        coverage <- apply(cover, 1, mean)
        
        boot_cover <- sapply(results[, 4], function(cis) {
            wald_cover <- (cis$wald_ci[1] <= optimal.auc) & (optimal.auc <= cis$wald_ci[2])
            perc_cover <- (cis$perc_ci[1] <= optimal.auc) & (optimal.auc <= cis$perc_ci[2])
            stud_cover <- (cis$stud_ci[1] <= optimal.auc) & (optimal.auc <= cis$stud_ci[2])
            efron_cover <- (cis$efron_ci[1] <= optimal.auc) & (optimal.auc <= cis$efron_ci[2])
            boot_wald_cover <- (cis$boot_wald_ci[1] <= optimal.auc) & (optimal.auc <= cis$boot_wald_ci[2])
            return(
                c(wald_cover, perc_cover, stud_cover, efron_cover, boot_wald_cover)
            )
        })
        boot_coverage <- matrix(apply(boot_cover, 1, mean), nrow=1)
        boot_coverage <- data.frame(
            "wald"=boot_coverage[1],
            "percentile"=boot_coverage[2],
            "percentile-t"=boot_coverage[3],
            "efron"=boot_coverage[4],
            "boot_wald"=boot_coverage[5],
            "size"=n,
            "alog" = ranking.algorithms[[j]]
        )
        print(boot_coverage)
        final.boot_coverage <- rbind(final.boot_coverage, boot_coverage)
        
        rank_cover <- sapply(results[, 3], function(rank){
            cover_flag <- FALSE
            for(i in seq_along(optimal.rank)) {
                if(all(rank[1:5] == optimal.rank[[i]])) {
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



# save(average.res, final.boot_coverage, file="bootstrap_average.res_April_27.RData")


# ================================= Plot ====================================

# load(file="bootstrap_average.res_April_27.RData")

library(ggplot2)
library(reshape2)
require(gridExtra)
library(ggpubr)
library(latex2exp)

colnames(final.boot_coverage) <- c("Wald", "Percentile", "Percentile-t", "Efron", "Bootstrap Wald", "size", "algo")
final.boot_coverage <- final.boot_coverage[, c(1, 5, 2, 3, 4, 6, 7)]

sub_boot_coverage <- final.boot_coverage[final.boot_coverage$algo == "gam", ]
fig1 <- ggplot(data=melt(sub_boot_coverage, id=c("size", "algo"), value.name = "coverage")) +
    geom_line(aes(x=size, y=coverage, colour=variable, linetype=variable)) +
    geom_point(aes(x=size, y=coverage, colour=variable)) +
    geom_hline(yintercept = 0.95, linetype="dashed") +
    coord_cartesian(ylim=c(0.85, 1)) +
    scale_y_continuous(breaks=c(0.85, 0.9, 0.95, 1.0),
                       label=c("0.85", "0.90", "0.95", "1.00")
    ) +
    scale_x_continuous(
        breaks=c(500, 1000, 2000, 3000, 4000, 5000), 
        label=c("500", "1K", "2K", "3K", "4K", "5K")
    ) +
    ylab("Coverage of 95% CI") + xlab("") +
    labs(color='Ranking algorithm', linetype='Ranking algorithm') +
    theme(legend.key.width = unit(2, "line"), panel.spacing = unit(0.5, "lines"), legend.position="bottom",
          # strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12), 
          # text = element_text(size = 10)
    ) + 
    guides(color=guide_legend(nrow=1))

sub_boot_coverage <- final.boot_coverage[final.boot_coverage$algo == "linear", ]
fig2 <- ggplot(data=melt(sub_boot_coverage, id=c("size", "algo"), value.name = "coverage")) +
    geom_line(aes(x=size, y=coverage, colour=variable, linetype=variable)) +
    geom_point(aes(x=size, y=coverage, colour=variable)) +
    geom_hline(yintercept = 0.95, linetype="dashed") +
    coord_cartesian(ylim=c(0.85, 1)) +
    scale_y_continuous(breaks=c(0.85, 0.9, 0.95, 1.0),
                       label=c("0.85", "0.90", "0.95", "1.00")
    ) +
    scale_x_continuous(
        breaks=c(500, 1000, 2000, 3000, 4000, 5000), 
        label=c("500", "1K", "2K", "3K", "4K", "5K")
    ) +
    ylab("") + xlab("Sample size") +
    labs(color='Ranking algorithm', linetype='Ranking algorithm') +
    theme(legend.key.width = unit(2, "line"), panel.spacing = unit(0.5, "lines"), legend.position="bottom",
          # strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12), 
          # text = element_text(size = 10)
    ) + 
    guides(color=guide_legend(nrow=1))

sub_boot_coverage <- final.boot_coverage[final.boot_coverage$algo == "earth", ]
fig3 <- ggplot(data=melt(sub_boot_coverage, id=c("size", "algo"), value.name = "coverage")) +
    geom_line(aes(x=size, y=coverage, colour=variable, linetype=variable)) +
    geom_point(aes(x=size, y=coverage, colour=variable)) +
    geom_hline(yintercept = 0.95, linetype="dashed") +
    coord_cartesian(ylim=c(0.85, 1)) +
    scale_y_continuous(breaks=c(0.85, 0.9, 0.95, 1.0),
                       label=c("0.85", "0.90", "0.95", "1.00")
    ) +
    scale_x_continuous(
        breaks=c(500, 1000, 2000, 3000, 4000, 5000), 
        label=c("500", "1K", "2K", "3K", "4K", "5K")
    ) +
    ylab("") + xlab("") +
    labs(color='Ranking algorithm', linetype='Ranking algorithm') +
    theme(legend.key.width = unit(2, "line"), panel.spacing = unit(0.5, "lines"), legend.position="bottom",
          # strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12), 
          # text = element_text(size = 10)
    ) + 
    guides(color=guide_legend(nrow=1))

ggarrange(fig1, fig2, fig3, ncol=3, common.legend = TRUE, legend="bottom", widths = c(1, 1, 1),
          legend.grob=get_legend(fig3))
    








