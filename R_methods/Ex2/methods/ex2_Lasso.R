#' Lasso-logit Method for Importance Sampling Variable Selection - Example 3
#' 
#' Uses L1-regularized logistic regression to identify important variables
#' and construct the importance sampling density.
#'
#' @param n_reps Number of independent replications (default: 25)
#' @param ncores Number of parallel cores (default: 25)
#' @param seed_base Base seed for reproducibility (default: 1234)
#' @param ROOT Root directory for sourcing shared functions
#' @param save_details If TRUE, returns detailed results
#' @return List of results for each replication

library(doParallel)

run_lasso <- function(n_reps = 25, 
                      ncores = 25, 
                      seed_base = 1235,
                      ROOT = getwd(),
                      save_details = FALSE) {
  
  # Set up parallel cluster
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl), add = TRUE)
  
  # Export ROOT to workers
  clusterExport(cl, "ROOT", envir = environment())
  
  # Run experiments
  POE_est <- foreach(j = 1:n_reps, 
                     .combine = c, 
                     .errorhandling = 'remove',
                     .packages = c("pracma", "mvtnorm", "cubature", "lhs", 
                                   "loo", "GauPro", "DiceKriging", "glmnet")) %dopar% {
    
    # Source shared functions inside each worker
    source(file.path(ROOT, "AMISE_bivariate.R"))
    source(file.path(ROOT, "functions_ex2.R"))
    source(file.path(ROOT, "Additive_S.R"))
    
    set.seed(seed_base + j)
    print(paste0('experiment:', j))
    
    # Parameters for Example 3
    D <- 5
    combo <- t(combn(seq(1:D), 2))
    Nm <- 200
    Nm_init <- 200
    l_target <- 24.2
    l <- l_target
    n_iter <- 5
    
    # Get initial samples
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
      # Fit Lasso-regularized logistic regression
      fit <- cv.glmnet(sample_X_all, sample_Z_all, alpha = 1, family = "binomial")
      
      # Extract coefficients at optimal lambda
      coef_best <- coef(fit, s = "lambda.min")
      
      S <- predict(fit, newx = sample_X_all, type = "response")
      
      index <- c(index, list(coef_best[-1, ]))
      
      # Rejection sampling with correlated structure for Example 3
      x <- c()
      k_all <- c()
      i <- 0
      
      while ((nrow(x) < Nm) || is.null(x)) {
        i <- i + 1
        # Example 3 sampling structure: X1, X4, X5 iid; X2|X1, X3|X1 conditional
        mu_i <- c(0, 0, 0)
        Sigma <- diag(3)
        x_i <- rmvnorm(1, mu_i, Sigma)
        X1 <- x_i[, 1]
        X2 <- rnorm(1, X1, 1)
        X3 <- rnorm(1, X1, 1)
        X4 <- x_i[, 2]
        X5 <- x_i[, 3]
        x_new <- cbind(X1, X2, X3, X4, X5)
        
        # Prior density
        prior_density <- dmvnorm(x_i) * dnorm(X2, X1, 1) * dnorm(X3, X1, 1)
        u <- runif(1, min = 0, max = prior_density)
        k <- predict(fit, newx = x_new, type = "response")
        
        if (u <= sqrt(k) * prior_density) {
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
      
      # Update sample bank
      sample_X_all <- rbind(sample_X_all, sample_X_new)
      sample_Z_all <- c(sample_Z_all, sample_Z_new)
    }
    
    # Return results
    # list(est_vec = est_vec, 
    #      sample_X_all = sample_X_all, 
    #      sample_Y_all = sample_Y_all, 
    #      index = index, 
    #      sample_S_all = sample_S_all)
    mean(est_vec)
  }
  
  return(POE_est)
}
