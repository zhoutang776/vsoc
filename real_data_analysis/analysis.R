rm(list = ls())
setwd("/Users/zhoutang/Onedrive/Research/VIM/code/vsoc")

library(SuperLearner)
source("vsoc.R")
set.seed(1234, kind = "L'Ecuyer-CMRG")
# =========================== Parameter Definition =============================
data <- read.csv("real_data_analysis/data/winequality-white.csv", sep=";")

n <- nrow(data)
p <- ncol(data) - 1
X <- data[, 1:p]
y <- data[, p+1]

X <- as.data.frame(scale(X))
y <- as.numeric(scale(y))

gam.ranking = function(X, y){
    X <- as.data.frame(X)
    formula <- as.formula(paste0("y~",paste0("s(", colnames(X),")",collapse="+")))
    model = mgcv::bam(formula, data=data.frame(y=y, X))
    temp = summary(model)
    rank = sort(temp$chi.sq, decreasing=TRUE, index.return=TRUE)$ix
    select = sum(temp$s.table[, 4] < 0.05)
    unloadNamespace("mgcv")
    return(list(rank=rank, select=select))
}
gam.rank <- gam.ranking(X, y)
    

lars.ranking = function(X, y){
    lars.fit <- lars::lars(as.matrix(X), y, type="lar")
    lars.order <- unlist(lars.fit$actions)
    
    glmnet.fit <- glmnet::cv.glmnet(as.matrix(X), y)
    lambda.min.index <- glmnet.fit$index[1]
    lambda.1se.index <- glmnet.fit$index[2]
    select <- c(glmnet.fit$nzero[lambda.min.index], glmnet.fit$nzero[lambda.1se.index])
    
    return(list(rank=as.vector(lars.order), select=as.numeric(select)))
}
lars.rank <- lars.ranking(X, y)


earth.ranking <- function(X, y) {
    model <- earth::earth(y ~ ., data = data.frame(y = y, X))
    str_rank <- row.names(earth::evimp(model, trim = FALSE))
    rank <- match(sub("-unused$", "", str_rank), colnames(X))
    select <- sum(!endsWith(str_rank, "unused"))
    return(list(rank=rank, select=select))
}
earth.rank <- earth.ranking(X, y)



# =========================== Parameter Definition =============================
SL.gam.customize <- function(...){
  SL.gam(..., deg.gam = 5)
}

SL.xgboost.customize <- function(...) {
    SL.xgboost(
        ...,
        shrinkage = 0.1,
        ntrees = 200,
        # max_depth = 10,
        # minobspernode = 5,
        nthread = 7,
    )
}
learners <- c("SL.xgboost.customize", "SL.gam.customize", "SL.earth")
V <- 5

alpha <- 0.05
seed <- 1235
bootstrap <- F


# ================================= Running ====================================
start <- proc.time()[3]
lars.result <- vsoc(X = X, y = y, alpha = alpha, V = V, bootstrap=bootstrap,
               rank = lars.rank$rank, learners = learners)
end <- proc.time()[3]
lars.time <- end - start

start <- proc.time()[3]
gam.result <- vsoc(X = X, y = y, alpha = alpha, V = V, bootstrap=bootstrap,
                    rank = gam.rank$rank, learners = learners)
end <- proc.time()[3]
gam.time <- end - start

start <- proc.time()[3]
earth.result <- vsoc(X = X, y = y, alpha = alpha, V = V, bootstrap=bootstrap,
                    rank = earth.rank$rank, learners = learners)
end <- proc.time()[3]
earth.time <- end - start
  


# save(lars.result, gam.result, earth.result, file = "real_data.Rdata")
load("real_data.RData")


lars.result$auc_est - qnorm(0.975) * sqrt(lars.result$auc_var/n)
lars.result$auc_est + qnorm(0.975) * sqrt(lars.result$auc_var/n)

earth.result$auc_est - qnorm(0.975) * sqrt(earth.result$auc_var/n)
earth.result$auc_est + qnorm(0.975) * sqrt(earth.result$auc_var/n)

gam.result$auc_est - qnorm(0.975) * sqrt(gam.result$auc_var/n)
gam.result$auc_est + qnorm(0.975) * sqrt(gam.result$auc_var/n)


# ================================= rank Plot ====================================
library(ggplot2)

final_psi <- mean(lars.result$psi_est[11], earth.result$psi_est[11], gam.result$psi_est[11])
lars.result$psi_est[11] <- final_psi
earth.result$psi_est[11] <- final_psi
gam.result$psi_est[11] <- final_psi


lars_uniform_ci <- cbind(
    lars.result$psi_est - 2.569214 * sqrt(lars.result$psi_var/n),
    lars.result$psi_est + 2.569214 * sqrt(lars.result$psi_var/n)
)

earth_uniform_ci <- cbind(
    earth.result$psi_est - 2.569214 * sqrt(earth.result$psi_var/n),
    earth.result$psi_est + 2.569214 * sqrt(earth.result$psi_var/n)
)

gam_uniform_ci <- cbind(
    gam.result$psi_est - 2.569214 * sqrt(gam.result$psi_var/n),
    gam.result$psi_est + 2.569214 * sqrt(gam.result$psi_var/n)
)


temp <- rbind(
    cbind(lars.result$psi_est, lars.result$psi_ci, lars_uniform_ci), 
    cbind(earth.result$psi_est, earth.result$psi_ci, earth_uniform_ci), 
    cbind(gam.result$psi_est, gam.result$psi_ci, gam_uniform_ci) 
)


vroc_est <- data.frame(temp, rep(c("LASSO", "MARS", "GAM"), each=11), 
                       rep(1:11, times=3))
colnames(vroc_est) <- c("est", "lb", "ub", "uniform_lb", "uniform_ub", "algo", "nvar")

vroc_est[, 1:5] <- pmax(vroc_est[, 1:5], 0)
vroc_est[is.na(vroc_est)] <- 0

vrank <- ggplot(data = vroc_est) +
    geom_line(aes(nvar, est, colour = algo, linetype=algo)) +
    geom_line(aes(x=nvar, y=uniform_lb, colour=algo, linetype=algo), alpha=1) +
    geom_line(aes(x=nvar, y=uniform_ub, colour=algo, linetype=algo), alpha=1) +
    geom_point(aes(nvar, est, colour = algo)) +
    geom_errorbar(mapping=aes(x=nvar, ymin=lb, ymax=ub, colour=algo, linetype=NULL), width=0.1) +
    coord_cartesian(ylim=c(0, 0.5)) +
    scale_x_continuous(breaks = 1:11) +
    ylab("Estimated R-squared predictiveness measures") + 
    xlab("Number of selected covariates") +
    labs(color='Ranking or selection algorithm', linetype='Ranking or selection algorithm') +
    theme_bw() + 
    theme(legend.key.width = unit(3, "line"), panel.spacing = unit(0.5, "lines"), legend.position="bottom",
    )

# ================================= selection Plot ====================================

lars.select <- lars.rank$select[1]
gam.select <- gam.rank$select
earth.select <- earth.rank$select

original <- c(0, 0,NA,NA)

temp <- rbind(
    cbind(lars.select, lars.result$psi_est[lars.select], t(lars.result$psi_ci[lars.select, ])), 
    original,     
    cbind(earth.select, earth.result$psi_est[earth.select], t(earth.result$psi_ci[earth.select, ])), 
    original, 
    cbind(gam.select, gam.result$psi_est[gam.select], t(gam.result$psi_ci[gam.select, ])),
    original
)

selected.vsoc <- data.frame(temp, rep(c("LASSO", "MARS", "GAM"), each=2))
colnames(selected.vsoc) <- c("nvar", "psi", "lb", "ub", "algo")


vselect <- ggplot(data=selected.vsoc) +
    geom_line(aes(x=nvar, y=psi, colour=algo, linetype=algo)) +
    geom_point(aes(x=nvar, y=psi, colour=algo)) +
    geom_errorbar(mapping=aes(x=nvar, ymin=lb, ymax=ub, colour=algo, linetype=NULL), width=0.1) +
    coord_cartesian(ylim=c(0, 0.5)) +
    scale_x_continuous(breaks=c(0:11)) + 
    ylab("") + xlab("Number of selected covariates") +
    labs(color='Ranking or selection algorithm', linetype='Ranking or selection algorithm') +
    theme_bw() + 
    theme(legend.key.width = unit(3, "line"), panel.spacing = unit(0.5, "lines"), legend.position="bottom",
    )


ggpubr::ggarrange(vrank, vselect, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")






# ================================== test
auc_test <- function(ranking1, ranking2){
    est1 <- ranking1$auc_est
    est2 <- ranking2$auc_est
    n <- dim(ranking1$auc_eif)[1]
    
    var1 <- ranking1$auc_var
    var2 <- ranking2$auc_var
    
    
    threshold <- qnorm(0.975) * sqrt((var1 + var2)/n)
    return(abs(est1 - est2) > threshold)
}

auc_test <- function(ranking1, ranking2){
    est1 <- ranking1$auc_est
    est2 <- ranking2$auc_est
    
    eif1 <- ranking1$auc_eif
    eif2 <- ranking2$auc_eif
    eif <- eif1 - eif2
    var <- mean(eif^2)
    n <- dim(ranking1$auc_eif)[1]

    statistic <- (est1 - est2) / sqrt(var/n)

    return(1-pchisq(statistic^2, 1))
}

auc_test(lars.result,  gam.result)
auc_test(lars.result,  earth.result)
auc_test(gam.result,  earth.result)



matrixsqrtinv <- function(A3){
    X = eigen(A3)
    T = X$vectors
    J = diag( x=X$values )
    Jinvsqrt = diag( x=1/sqrt( X$values ) )
    Tinv = solve(T)
    return(T %*% Jinvsqrt %*% Tinv)
}

vroc_test <- function(ranking1, ranking2){
    est1 <- ranking1$psi_est
    est2 <- ranking2$psi_est
    
    est1 <- pmax(est1, 0)
    est2 <- pmax(est2, 0)
    
    eif1 <- ranking1$psi_eif
    eif2 <- ranking2$psi_eif
    
    n <- dim(eif1)[1]
    
    eif <- eif1 - eif2
    sqrt.inv.cov <- matrixsqrtinv( t(eif) %*% eif / n)
    
    statistic <- sum((sqrt.inv.cov %*% c(est1 - est2))^2)
    
    # threshold <- qnorm(0.975) * sqrt((var1 + var2)/n)
    # return(any(abs(est1 - est2) > threshold))
    return(1-pchisq(statistic, df=11))
}

vroc_test <- function(ranking1, ranking2){
    est1 <- ranking1$psi_est
    est2 <- ranking2$psi_est
    
    # est1 <- pmax(est1, 0)
    # est2 <- pmax(est2, 0)
    
    eif1 <- ranking1$psi_eif
    eif2 <- ranking2$psi_eif
    
    n <- dim(eif1)[1]
    
    eif <- eif1 - eif2
    covaraince <- t(eif) %*% eif / n
    # covaraince <- cov(eif)
    
    statistic <- n * t(c(est1 - est2)) %*% solve(covaraince) %*% c(est1 - est2)
    
    return(1-pchisq(statistic, df=11))
}

vroc_test(lars.result,  gam.result)
vroc_test(lars.result,  earth.result)
vroc_test(gam.result,  earth.result)



