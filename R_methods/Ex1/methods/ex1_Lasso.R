#' Lasso-logit Method for Importance Sampling Variable Selection
#' 
#' Uses L1-regularized logistic regression to identify important variables
#' and construct the importance sampling density.
#'
#' param n_reps Number of independent replications (default: 25)
#' param ncores Number of parallel cores (default: 25)
#' param seed_base Base seed for reproducibility (default: 123)
#' param ROOT Root directory for sourcing shared functions
#' param save_details If TRUE, returns detailed results; otherwise just POE estimates
#' return Vector of POE estimates (or list with details if save_details = TRUE)

library(doParallel)

run_lasso <- function(n_reps = 25, 
                      ncores = 25, 
                      seed_base = 123,
                      ROOT = getwd(),
                      save_details = FALSE) {
  
  # Set up parallel cluster
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl), add = TRUE)  # ensures cleanup even on error
  
  # Export ROOT to workers so they can source files
  clusterExport(cl, "ROOT", envir = environment())
  
  # Run experiments
  POE_est <- foreach(j = 1:n_reps, 
                     .combine = c, 
                     .errorhandling = 'remove',
                     .packages = c("pracma", "mvtnorm", "cubature", 
                                   "lhs", "loo", "GauPro", 
                                   "DiceKriging", "glmnet")) %dopar% {
    
    # Source shared functions inside each worker
    source(file.path(ROOT, "AMISE_bivariate.R"))
    source(file.path(ROOT, "functions_ex1.R"))
    source(file.path(ROOT, "Additive_S.R"))
    
    set.seed(seed_base + j)
    
    # Parameters
    D <- 4
    Nm <- 200
    Nm_init <- 200
    l_target <- 18.99
    l <- l_target
    n_iter <- 5
    
    # Get initial samples
    sample_X <- randomLHS(Nm_init, D) * 10 - 5
    sample_Y <- get_Y(sample_X)
    sample_Z <- (sample_Y > l)
    sample_X_all <- sample_X
    sample_Z_all <- sample_Z
    sample_Y_all <- sample_Y
    
    est_vec <- c()
    index <- c()
    sample_S_all <- c()
    
    for (simu in 1:n_iter) {
      # Fit Lasso-regularized logistic regression
      fit <- cv.glmnet(sample_X_all, sample_Z_all, alpha = 1, family = "binomial")
      
      # Extract coefficients at optimal lambda
      coef_best <- coef(fit, s = "lambda.min")
      index <- c(index, list(coef_best[-1, ]))
      
      # Rejection sampling from importance density
      x <- c()
      k_all <- c()
      i <- 0
      
      while ((nrow(x) < Nm) || is.null(x)) {
        i <- i + 1
        x_new <- rmvnorm(1, rep(0, D))
        u <- runif(1, min = 0, max = dmvnorm(x_new))
        k <- predict(fit, newx = x_new, type = "response")
        
        if (u <= sqrt(predict(fit, newx = x_new, type = "response")) * dmvnorm(x_new)) {
          x <- rbind(x, x_new)
          k_all <- c(k_all, k)
        }
      }
      
      # Compute estimates with new samples
      sample_X_new <- x
      sample_Y_new <- get_Y(sample_X_new)
      sample_Z_new <- (sample_Y_new > l)
      sample_S_new <- k_all
      sample_S_all <- c(sample_S_all, sample_S_new)
      sample_H_new <- sqrt(sample_S_new)
      est_temp_new <- get_Est(sample_Y_new, sample_H_new, l_target)
      est_vec <- c(est_vec, est_temp_new)
      
      # Update sample bank
      sample_X_all <- rbind(sample_X_all, sample_X_new)
      sample_Y_all <- c(sample_Y_all, sample_Y_new)
      sample_Z_all <- c(sample_Z_all, sample_Z_new)
    }
    
    # Return mean estimate for this replication
    mean(est_vec)
  }
  
  return(POE_est)
}
