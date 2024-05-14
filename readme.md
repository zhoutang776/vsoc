# VSOC: Nonparametric Assessment of Variable Selection and Ranking Algorithms

## Introduction
Selecting from or ranking a set of candidates variables in terms of their capacity
for predicting an outcome of interest is an important task in many scientific fields.
A variety of methods for variable selection and ranking have been proposed in the
literature. In practice, it can be challenging to know which method is most appropriate for a given dataset. We propose methods of comparing variable selection and ranking algorithms. We first introduce measures of the quality of variable selection and ranking algorithms. We then define estimators of our proposed
measures, and establish asymptotic results for our estimators. 

More detail may be found in our papers on [Nonparametric Assessment of Variable Selection and Ranking Algorithms](https://arxiv.org/abs/2308.11593)


## Example
```
rm(list = ls())
library(SuperLearner)
set.seed(1234, kind = "L'Ecuyer-CMRG")
source("vsoc.R")

# =========================== Function Definition ==============================
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

# linear rank
lasso.ranking <- function(X, y){
    lars.fit <- lars::lars(as.matrix(X), y, type = "lar")
    rank <- unlist(lars.fit$actions)
    stopifnot(length(rank) == ncol(X))
    
    glmnet.fit <- glmnet::cv.glmnet(as.matrix(X), y)
    lambda.min.index <- glmnet.fit$index[1]
    lambda.1se.index <- glmnet.fit$index[2]
    select <- c(glmnet.fit$nzero[lambda.min.index], glmnet.fit$nzero[lambda.1se.index])
    
    return(list(rank=as.vector(rank), select=as.numeric(select)))
}

# =========================== Run VSOC ==============================
alpha <- 0.05
V <- 10
n <- 500
learners <- c("SL.xgboost.customize", "SL.gam.customize", "SL.earth")

data <- generate.samples(n = n)
X <- as.data.frame(scale(data$X))
y <- scale(data$y)

rank.select <- lasso.ranking(X, y)
result <- vsoc(X = X, y = y, alpha = alpha, V = V, bootstrap = F, B = NULL,
               rank = rank.select$rank, learners = learners)
```



