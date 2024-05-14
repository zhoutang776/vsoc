vsoc <- function(X, y, rank = NULL, rank.algorithm = NULL,
                 V = 5, learners = c("SL.glmnet", "SL.gam", "SL.mean"),
                 alpha = 0.05, bootstrap = FALSE, B = 500){
    n <- nrow(X); p <- ncol(X)
    if(is.null(dim(y))) {
        y <- as.matrix(y)
    }

    # stopifnot(n %% V == 0)
    stopifnot(!is.null(rank) || !is.null(rank.algorithm))
    
    select <- NULL
    if(is.null(rank)) {
        rank.select <- rank.algorithm(X, y)
        rank <- rank.select$rank
        select <- rank.select$select
    }
    # this list stores all of the fitted learners. each learner is fitted on the full dataset instead of
    # bootstrap dataset
    all_y_pred <- list()
    # compute the full regression function firstly because this is the same for all of data and bootstrap data
    indx <- 1:p
    str_indx <- paste(as.character(indx), collapse = ",")
    all_y_pred[[str_indx]] <- vsoc_regression(X = X, y = y, V = V, learners = learners)$y_pred
    y_pred_full <- matrix(all_y_pred[[str_indx]], ncol=1)
    y_pred_redu <- matrix(0, nrow = n, ncol = p)
    
    for(i in 1:p) {
        # calculate the variable importance of first i variables, which means we need to exclude the first i variables
        if(i != p) {
            indx <- sort(rank[(i + 1):p])
            str_indx <- paste(as.character(indx), collapse = ",")
        } else {
            indx <- NULL
            str_indx <- "NULL"
        }
        all_y_pred[[str_indx]] <- vsoc_regression(
            X = X[, indx, drop = FALSE], y = y, learners = learners, V = V
        )$y_pred
        y_pred_redu[, i] <- all_y_pred[[str_indx]]
    }
    
    temp <- vsoc_est_eif(y_pred_full = y_pred_full, y_pred_redu = y_pred_redu, y = y, cv_fold = cv_fold)
    
    mse_full_est <- temp$mse_full_est
    mse_full_eif <- temp$mse_full_eif
    mse_full_var <- apply(mse_full_eif^2, 2, mean)
    mse_full_ci <- mse_full_est + sqrt(mse_full_var / n) * qnorm(c(alpha / 2, 1 - alpha / 2))
    
    
    mse_redu_est <- temp$mse_redu_est
    mse_redu_eif <- temp$mse_redu_eif
    mse_redu_var <- apply(mse_redu_eif^2, 2, mean)
    mse_redu_ci <- mse_redu_est + sqrt(mse_redu_var / n) %o% qnorm(c(alpha / 2, 1 - alpha / 2))
    
    theta_est <- temp$theta_est
    theta_eif <- temp$theta_eif
    theta_var <- apply(theta_eif^2, 2, mean)
    theta_ci <- theta_est + sqrt(theta_var / n) %o% qnorm(c(alpha / 2, 1 - alpha / 2))
    
    psi_est <- temp$psi_est
    psi_eif <- temp$psi_eif
    # compute the asymptotic covariance matrix
    psi_cov <- (t(psi_eif) %*% psi_eif) / n
    psi_var <- diag(psi_cov)
    psi_ci <- psi_est + sqrt(psi_var / n) %o% qnorm(c(alpha / 2, 1 - alpha / 2))
    
    auc_est <- temp$auc_est
    auc_eif <- temp$auc_eif
    auc_cov <- (t(auc_eif) %*% auc_eif) / n
    auc_var <- diag(auc_cov)
    # calculate the confidence interval
    auc_ci <- auc_est + sqrt(auc_var / n) * qnorm(c(alpha / 2, 1 - alpha / 2))
    
    # psi_est <- isoreg(psi_est)$yf
    # if est < 0, set to zero and print warning
    
    bootstrap_auc_est <- NULL
    bootstrap_auc_var <- NULL
    bootstrap_psi_est <- NULL
    bootstrap_psi_var <- NULL
    bootstrap_rank_list <- list()
    bootstrap_select_list <- list()
    if(bootstrap) {
        # means calculating the se from bootstrap.
        bootstrap_auc_est <- matrix(NA, nrow = B, ncol = 1)
        bootstrap_auc_var <- matrix(NA, nrow = B, ncol = 1)
        
        bootstrap_psi_est <- matrix(NA, nrow = B, ncol = p)
        bootstrap_psi_var <- matrix(NA, nrow = B, ncol = p)
        
        for(b in 1:B) {
            bootstrap_index <- sample(n, n, replace = TRUE)

            bootstrap_y <- y[bootstrap_index,]
            bootstrap_X <- X[bootstrap_index,]
            bootstrap_rank.select <- rank.algorithm(X = bootstrap_X, y = bootstrap_y)
            bootstrap_rank <- bootstrap_rank.select$rank
            bootstrap_select <- bootstrap_rank.select$select
            
            bootstrap_rank_list[[b]] <- bootstrap_rank
            bootstrap_select_list[[b]] <- bootstrap_select
                
            # compute full predictiveness
            bootstrap_y_pred_full <- matrix(0, nrow = n, ncol = 1)
            str_indx <- paste(as.character(1:p), collapse = ",")
            bootstrap_y_pred_full[, 1] <- all_y_pred[[str_indx]][bootstrap_index]

            # compute reduced predictiveness recursively
            bootstrap_y_pred_redu <- matrix(0, nrow = n, ncol = p)
            for(i in 1:p) {
                if(i == p) {
                    indx <- NULL
                    str_indx <- "NULL"
                } else {
                    indx <- sort(bootstrap_rank[(i + 1):p])
                    str_indx <- paste(as.character(indx), collapse = ",")
                }
                if(!(str_indx %in% names(all_y_pred))) {
                    # means we need to fit a learner before we calculate the predict!
                    all_y_pred[[str_indx]] <- vsoc_regression(
                        X = X[, indx, drop = FALSE], y = y, V = V, learners = learners
                    )$y_pred
                }
                bootstrap_y_pred_redu[, i] <- all_y_pred[[str_indx]][bootstrap_index]
            }

            bootstrap_temp <- vsoc_est_eif(y_pred_full = bootstrap_y_pred_full,
                                           y_pred_redu = bootstrap_y_pred_redu,
                                           y = bootstrap_y)
            
            bootstrap_auc_est[b,] <- bootstrap_temp$auc_est
            bootstrap_auc_eif <- bootstrap_temp$auc_eif
            bootstrap_auc_var[b, ] <- apply(bootstrap_auc_eif^2, 2, mean)
            bootstrap_psi_est[b, ] <- bootstrap_temp$psi_est
            bootstrap_psi_var[b, ] <- apply(bootstrap_temp$psi_eif^2, 2, mean)
            
        }
    } 

    return(
        list(
            psi_est = psi_est, psi_var = psi_var, psi_ci = psi_ci, psi_cov = psi_cov, psi_eif=psi_eif,
            auc_est = auc_est, auc_var = auc_var, auc_ci = auc_ci, auc_eif=auc_eif,
            theta_est = theta_est, theta_var = theta_var, theta_ci = theta_ci, theta_eif=theta_eif,
            mse_full_est = mse_full_est, mse_full_ci = mse_full_ci, mse_full_var = mse_full_var, mse_full_eif=mse_full_eif,
            mse_redu_est = mse_redu_est, mse_redu_ci = mse_redu_ci, mse_redu_var = mse_redu_var, mse_redu_eif=mse_redu_eif,
            rank = rank, select = select, 
            bootstrap_rank_list = bootstrap_rank_list, bootstrap_select_list = bootstrap_select_list, 
            bootstrap_auc_est = bootstrap_auc_est, bootstrap_auc_var = bootstrap_auc_var,
            bootstrap_psi_est = bootstrap_psi_est, bootstrap_psi_var = bootstrap_psi_var
        )
    )

}


vsoc_regression <- function(X, y, learners = c("SL.glmnet", "SL.gam", "SL.mean"), V = 5, cv_fold = NULL){
    n <- nrow(X)
    p <- ncol(X)
    if(is.null(dim(y))) {
        y <- as.matrix(y)
    }
    # set up the cross-fitting
    if(is.null(cv_fold)) cv_fold <- .make_folds(y, V = V)
    else{ stopifnot(length(cv_fold) == n)}

    y_pred <- matrix(0, nrow = n, ncol = 1)
    for(v in 1:V) {
        # fit super learner on the full covariates
        if(V == 1) {
            # means do not use cross validation to estimate the nuisance parameter
            train_index <- rep(TRUE, n)
            test_index <- rep(TRUE, n)
        } else {
            train_index <- cv_fold != v
            test_index <- cv_fold == v
        }
        if(p == 0) {
            y_pred[test_index,] <- mean(y[train_index,])
            next
        }

        if(length(learners) == 1) {
            # If we have multiple learners, SuperLearner will use cross validation to determine the
            # coefficient of each learner. It will be computationally efficient if we call the learner
            # directly when there only have one learner.

            learner <- eval(parse(text = learners))
            model <- learner(Y = y[train_index],
                             X = X[train_index, , drop = FALSE],
                             newX = X[test_index, , drop = FALSE],
                             family = stats::gaussian(),
                             obsWeights = rep(1, sum(train_index))
            )
            y_pred[test_index,] <- model$pred
        } else {
            # use multiple learners to fit the response variable
            model <- SuperLearner::SuperLearner(
                Y = y[train_index, ],
                X = X[train_index, , drop = FALSE],
                newX = X[test_index, , drop = FALSE],
                family = stats::gaussian(),
                SL.library = learners
            )
            # print(model$coef)
            # get predictions on the validation fold
            y_pred[test_index,] <- model$SL.predict
        }
    }
    return(list(y_pred = y_pred, cv_fold = cv_fold))
}


vsoc_est_eif <- function(y_pred_full, y_pred_redu, y, ...){
    # get two kind of estimator and their eif
    # unscale version theta with its eif
    # e.g. theta = int (mu_full - mu_redu)^2 dP = mse_redu - mse_full,
    # where mse_redu = int (y - mu_redu)^2 dP, mse_full = int (y - mu_full)^2 dP
    #
    # scale version psi with its eif
    # e.g. psi = theta / sigma^2
    
    p <- ncol(y_pred_redu)
    n <- nrow(y_pred_redu)
    
    vec <- rep(1, p)
    vec[1] <- 0.5
    vec[p] <- 0.5
    mn_y <- mean(y)
    var_y <- mean((y - mn_y)^2)
    sigma_eif <- (y - mn_y)^2 - var_y
    
    mse_full_est <- mean((y - y_pred_full)^2)
    mse_full_eif <- (y - y_pred_full)^2 - mse_full_est
    
    mse_redu_est <- rep(0, p)
    mse_redu_eif <- matrix(0, nrow=n, ncol=p)
    
    theta_est <- rep(0, p)
    theta_eif <- matrix(0, nrow=n, ncol=p)
    
    psi_est <- rep(0, p)
    psi_eif <- matrix(0, nrow=n, ncol=p)
    
    for(i in 1:p){
        mse_redu_est[i] <- mean((y - y_pred_redu[, i])^2)
        mse_redu_eif[, i] <- (y - y_pred_redu[, i])^2 - mse_redu_est[i]
        
        theta_est[i] <- mse_redu_est[i] - mse_full_est
        theta_eif[, i] <- mse_redu_eif[, i] - mse_full_eif
        
        psi_est[i] <- theta_est[i] / var_y
        psi_eif[, i] <- theta_eif[, i] / var_y + sigma_eif * (-theta_est[i] / var_y^2)
    }
    auc_est <- sum(psi_est * vec)
    auc_eif <- psi_eif %*% vec
    
    return(
        list(auc_est = auc_est, auc_eif = auc_eif,
             psi_est = psi_est, psi_eif = psi_eif,
             theta_est = theta_est, theta_eif = theta_eif,
             mse_full_est = mse_full_est, mse_full_eif = mse_full_eif,
             mse_redu_est = mse_redu_est, mse_redu_eif = mse_redu_eif)
    )
}


vsoc_ci <- function(n, est, var = NULL, alpha = 0.05, bootstrap = FALSE){
    if(nrow(est) == n){
        return(t(apply(est, 2, function(x) quantile(x, c(alpha / 2, 1 - alpha / 2)))))
    } else {
        return(est + sqrt(var / n) %o% qnorm(c(alpha / 2, 1 - alpha / 2)))
    }
}

.make_folds <- function(y, V = 2, seed = NULL){
    n <- length(y)
    folds <- sample(rep(1:V, length.out = n))
    return(folds)
}

get_wald_ci <- function(alpha=0.05, bs_est=NULL, bs_est0=NULL, bs_se=NULL, est=NULL, se=NULL){
    if(is.null(bs_est) || is.null(bs_est0)){
        stop("wrong arguments in Wald method")
    }
    if(is.null(est)) est <- bs_est0
    
    bs_bias <- 0
    wald_ci <- est - bs_bias + qnorm(c(alpha/2, 1-alpha/2)) * sd(bs_est) 
    
    return(as.numeric(wald_ci))
}

get_efron_ci <- function(alpha=0.05, bs_est=NULL, bs_est0=NULL, bs_se=NULL, est=NULL, se=NULL){
    if(is.null(bs_est) || is.null(bs_est0)){
        stop("wrong arguments in Efron's method")
    }
    
    efron_ci <- quantile(bs_est, c(alpha/2, 1-alpha/2), type = 6)
    
    return(as.numeric(efron_ci))
}

get_perc_ci <- function(alpha=0.05, bs_est=NULL, bs_est0=NULL, bs_se=NULL, est=NULL, se=NULL){
    if(is.null(bs_est) || is.null(bs_est0)){
        stop("wrong arguments in percentile method")
    }
    if(is.null(est)) est <- bs_est0
    
    # basic
    t_stats <- (bs_est - bs_est0)
    lb <- est - quantile(t_stats, c(1-alpha/2), type = 6)
    ub <- est - quantile(t_stats, c(alpha/2), type = 6) 
    perc_ci <- cbind(lb, ub)  
    
    return(as.numeric(perc_ci))
}

# only this will use the second estimator in the boot
get_stud_ci <- function(alpha=0.05, bs_est=NULL, bs_est0=NULL, bs_se=NULL, est=NULL, se=NULL){
    if(is.null(bs_est) || is.null(bs_est0) || is.null(se) || is.null(bs_se)){
        stop("wrong arguments in percentile-t method")   
    }
    if(is.null(est)) est <- bs_est0
    
    
    t_stats <- (bs_est - bs_est0) / bs_se
    lb <-  est - quantile(t_stats, 1-alpha/2, type=6) * se
    ub <-  est - quantile(t_stats, alpha/2, type=6) * se
    stud_ci <- cbind(lb, ub) 
    
    return(as.numeric(stud_ci))
}


get_boot_ci <- function(alpha=0.05, bs_est=NULL, bs_est0=NULL, bs_se=NULL, est=NULL, se=NULL){
    # boot_out$t0: est, variance
    
    # Normal
    wald_ci <- get_wald_ci(alpha=alpha, bs_est=bs_est, bs_est0=bs_est0, bs_se=bs_se, est=est, se=se) 
    # Basic
    perc_ci <- get_perc_ci(alpha=alpha, bs_est=bs_est, bs_est0=bs_est0, bs_se=bs_se, est=est, se=se) 
    # Studentized
    stud_ci <- get_stud_ci(alpha=alpha, bs_est=bs_est, bs_est0=bs_est0, bs_se=bs_se, est=est, se=se)
    # Percentile
    efron_ci <-get_efron_ci(alpha=alpha, bs_est=bs_est, bs_est0=bs_est0, bs_se=bs_se, est=est, se=se) 
    
    return(
        list(boot_wald_ci=wald_ci, perc_ci=perc_ci, stud_ci=stud_ci, efron_ci=efron_ci)
    )
}

