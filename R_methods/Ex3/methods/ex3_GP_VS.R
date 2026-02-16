#' GP-VS Method for Importance Sampling Variable Selection - Example 3
#' 
#' Uses Gaussian Process with greedy CE ordering and validation-based subset selection (D=10).
#'
#' @param n_reps Number of independent replications (default: 25)
#' @param ncores Number of parallel cores (default: 25)
#' @param seed_base Base seed for reproducibility (default: 123)
#' @param ROOT Root directory for sourcing shared functions
#' @return List of results for each replication

library(doParallel)

run_gp_vs <- function(n_reps = 25, 
                      ncores = 25, 
                      seed_base = 123,
                      ROOT = getwd()) {
  
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl), add = TRUE)
  clusterExport(cl, "ROOT", envir = environment())
  
  POE_est <- foreach(j = 1:n_reps, 
                     .combine = c, 
                     .errorhandling = 'remove',
                     .packages = c("pracma", "mvtnorm", "cubature", "lhs", 
                                   "loo", "DiceKriging", "GauPro")) %dopar% {
    
    source(file.path(ROOT, "AMISE_bivariate.R"))
    source(file.path(ROOT, "functions_ex3.R"))
    source(file.path(ROOT, "Additive_S.R"))
    source(file.path(ROOT, "cv.R"))
    
    set.seed(seed_base + j)
    print(paste0('experiment:', j))
    
    D <- 10
    combo <- t(combn(seq(1:D), 2))
    Nm <- 200
    Nm_init <- 200
    l_target <- 18.74
    l <- l_target
    n_iter <- 5
    
    sample_X <- randomLHS(Nm_init, D) * 10 - 5
    sample_Y <- get_Y(sample_X)
    sample_Z <- (sample_Y > l)
    sample_X_all <- sample_X
    sample_Z_all <- sample_Z
    sample_Y_all <- sample_Y
    
    est_vec <- c()
    index <- c()
    weight_ce_all <- c()
    h_s_pareto_all <- c()
    bandwidth <- c()
    sample_S_all <- c()
    index_all <- c()
    mse_all <- c()
    
    for (simu in 1:n_iter) {
      # Greedy CE ordering with GP
      h_s_all_sample_x <- sampling_ce_order_greedy_GP(Nm, sample_X_all, sample_Y_all, l, 8)
      h_s_all <- h_s_all_sample_x[[1]]
      print(h_s_all)
      sample_x_pareto_all <- h_s_all_sample_x[[2]]
      sample_s_pareto_all <- h_s_all_sample_x[[3]]
      
      h_s_pareto <- h_s_all
      
      # Validation-based subset selection with GP
      h_s_all_select <- S_find_index_validation_GP(sample_X_all, sample_Y_all, h_s_pareto, l, Nm)
      h_s_ce <- h_s_all_select[[1]]
      h_track <- h_s_all_select[[2]]
      
      index_pareto_list <- sapply(h_s_all, function(X) setequal(X, h_s_ce))
      index_pareto <- which(index_pareto_list == 1)
      sample_X_new <- sample_x_pareto_all[[index_pareto]]
      
      index <- c(index, list(h_s_ce))
      index_all <- c(index_all, list(h_track))
      h_s_pareto_all <- c(h_s_pareto_all, list(h_s_pareto))
      
      # Compute estimates
      mu_X_new <- mu(sample_X_new)
      p_all <- 1 - pnorm(l - mu_X_new)
      sample_Y_new <- get_Y(sample_X_new)
      sample_Z_new <- (sample_Y_new > l)
      sample_S_new <- sample_s_pareto_all[[index_pareto]]
      sample_H_new <- sqrt(sample_S_new)
      est_temp_new <- get_Est(sample_Y_new, sample_H_new, l_target)
      
      est_vec <- c(est_vec, est_temp_new)
      mse_all <- c(mse_all, sum((p_all - sample_S_new)^2) / 199)
      
      # Update sample bank
      sample_X_all <- rbind(sample_X_all, sample_X_new)
      sample_Z_all <- c(sample_Z_all, sample_Z_new)
      sample_Y_all <- c(sample_Y_all, sample_Y_new)
      sample_S_all <- c(sample_S_all, sample_S_new)
    }
    
    # list(est_vec = est_vec, sample_X_all = sample_X_all, sample_Y_all = sample_Y_all,
    #      index = index, h_s_pareto = h_s_pareto, index_all = index_all,
    #      sample_S_all = sample_S_all, mse_all = mse_all)
    mean(est_vec)
  }
  
  return(POE_est)
}
