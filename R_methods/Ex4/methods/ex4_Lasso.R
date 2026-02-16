#' Lasso Method for Importance Sampling Variable Selection - Example 4
#' 
#' Uses L1-regularized logistic regression for classification-based IS (D=50).
#' Independent sampling structure (iid normal).
#'
#' @param n_reps Number of independent replications (default: 25)
#' @param ncores Number of parallel cores (default: 25)
#' @param seed_base Base seed for reproducibility (default: 127)
#' @param ROOT Root directory for sourcing shared functions
#' @return List of results for each replication

library(doParallel)

run_lasso <- function(n_reps = 25, 
                      ncores = 25, 
                      seed_base = 127,
                      ROOT = getwd()) {
  
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl), add = TRUE)
  clusterExport(cl, "ROOT", envir = environment())
  
  POE_est <- foreach(j = 1:n_reps, 
                     .combine = c, 
                     .errorhandling = 'remove',
                     .packages = c("pracma", "mvtnorm", "cubature", "lhs", 
                                   "loo", "GauPro", "DiceKriging", "glmnet")) %dopar% {
    
    source(file.path(ROOT, "AMISE_bivariate.R"))
    source(file.path(ROOT, "functions_ex4.R"))
    source(file.path(ROOT, "Additive_S.R"))
    source(file.path(ROOT, "cv.R"))
    
    set.seed(seed_base + j)
    print(paste0('experiment:', j))
    
    D <- 50
    Nm <- 200
    Nm_init <- 200
    l_target <- 18.42
    l <- l_target
    n_iter <- 5
    
    sample_X <- randomLHS(Nm_init, D) * 10 - 5
    sample_Y <- get_Y(sample_X)
    sample_Z <- (sample_Y > l)
    sample_X_all <- sample_X
    sample_Z_all <- sample_Z
    
    est_vec <- c()
    index <- c()
    sample_Y_all <- sample_Y
    sample_S_all <- c()
    
    for (simu in 1:n_iter) {
      fit <- cv.glmnet(sample_X_all, sample_Z_all, alpha = 1, family = "binomial")
      
      coef_best <- coef(fit, s = "lambda.min")
      index <- c(index, list(coef_best[-1, ]))
      
      # Rejection sampling (iid normal)
      x <- c()
      k_all <- c()
      i <- 0
      
      while ((nrow(x) < Nm) || is.null(x)) {
        i <- i + 1
        x_new <- rmvnorm(1, rep(0, D))
        u <- runif(1, min = 0, max = dmvnorm(x_new))
        k <- predict(fit, newx = x_new, type = "response")
        
        if (u <= sqrt(k) * dmvnorm(x_new)) {
          x <- rbind(x, x_new)
          k_all <- c(k_all, k)
        }
      }
      
      sample_X_new <- x
      sample_Y_new <- get_Y(sample_X_new)
      sample_Z_new <- (sample_Y_new > l)
      sample_S_new <- k_all
      sample_S_all <- c(sample_S_all, sample_S_new)
      sample_H_new <- sqrt(sample_S_new)
      est_temp_new <- get_Est(sample_Y_new, sample_H_new, l_target)
      sample_Y_all <- c(sample_Y_all, sample_Y_new)
      est_vec <- c(est_vec, est_temp_new)
      
      sample_X_all <- rbind(sample_X_all, sample_X_new)
      sample_Z_all <- c(sample_Z_all, sample_Z_new)
    }
    
    # list(est_vec = est_vec, sample_X_all = sample_X_all, sample_Y_all = sample_Y_all,
    #     index = index, sample_S_all = sample_S_all)
    mean(est_vec)
  }
  
  return(POE_est)
}
