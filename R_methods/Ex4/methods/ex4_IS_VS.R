#' IS-VS Method for Importance Sampling Variable Selection (Proposed Method) - Example 4
#' 
#' Uses greedy CE ordering with Pareto-k filtering and validation-based subset selection.
#' For D=50 (very high-dimensional), uses sampling_ce_order_greedy for efficient subset search.
#' No bandwidth optimization due to 1225 kernel pairs.
#'
#' @param n_reps Number of independent replications (default: 25)
#' @param ncores Number of parallel cores (default: 25)
#' @param seed_base Base seed for reproducibility (default: 127)
#' @param ROOT Root directory for sourcing shared functions
#' @return List of results for each replication

library(doParallel)

run_is_vs <- function(n_reps = 25, 
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
                     .packages = c("pracma", "mvtnorm", "cubature", "lhs", "loo")) %dopar% {
    
    source(file.path(ROOT, "AMISE_bivariate.R"))
    source(file.path(ROOT, "functions_ex4.R"))
    source(file.path(ROOT, "Additive_S.R"))
    source(file.path(ROOT, "cv.R"))
    
    set.seed(seed_base + j)
    print(paste0('experiment:', j))
    
    D <- 50
    combo <- t(combn(seq(1:D), 2))
    n_pairs <- 25 * 49  # C(50,2) = 1225
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
    
    h_init_vec <- c()
    for (i in 1:D) {
      h <- (optimize(function(x) AMISE(sample_X[, i], sample_Z, x, -5, 5), 
                     interval = c(0.1, 1.5)))$minimum
      h_init_vec[i] <- h
    }
    h_init_mat <- matrix(h_init_vec[combo], ncol = 2)
    
    est_vec <- c()
    index <- c()
    weight_ce_all <- c()
    h_s_pareto_all <- c()
    bandwidth <- c()
    sample_S_all <- c()
    sample_Y_all <- sample_Y
    index_all <- c()
    mse_all <- c()
    
    for (simu in 1:n_iter) {
      # No bandwidth optimization for D=50 (too many pairs)
      h_new_mat <- h_init_mat
      
      cost_vec <- c()
      S_all <- c()
      for (c in 1:nrow(combo)) {
        cb <- combo[c, ]
        h <- h_new_mat[c, ]
        S_temp <- S_bi(sample_X_all[, cb], sample_X_all[, cb], sample_Z_all, h[1], h[2])
        S_all <- c(S_all, S_temp)
        cost <- -sum(sample_Z_all * log(S_temp) + (1 - sample_Z_all) * log(1 - S_temp), na.rm = TRUE)
        cost_vec[c] <- cost
      }
      S_all <- matrix(S_all, nrow = length(S_all) / n_pairs)
      weight <- 1 / cost_vec / sum(1 / cost_vec)
      
      # Greedy CE ordering with Pareto-k filtering
      h_s_all_sample_x <- sampling_ce_order_greedy(Nm, sample_X_all, sample_Z_all, 
                                                    h_new_mat, combo, weight, S_all, 8)
      h_s_all <- h_s_all_sample_x[[1]]
      sample_x_pareto_all <- h_s_all_sample_x[[2]]
      
      # Validation-based subset selection
      h_s_pareto <- h_s_all
      h_s_all_select <- S_find_index_validation(sample_X_all, sample_Z_all, h_new_mat, h_s_pareto, weight, Nm)
      h_s_ce <- h_s_all_select[[1]]
      h_track <- h_s_all_select[[2]]
      
      index_pareto_list <- sapply(h_s_all, function(X) setequal(X, h_s_ce))
      index_pareto <- which(index_pareto_list == 1)
      sample_X_new <- sample_x_pareto_all[[index_pareto]]
      
      index <- c(index, list(h_s_ce))
      index_all <- c(index_all, list(h_track))
      h_s_pareto_all <- c(h_s_pareto_all, list(h_s_pareto))
      weight_ce <- weight / sum(weight[h_s_ce])
      weight_ce_all <- c(weight_ce_all, list(weight_ce))
      
      # Compute estimates
      mu_X_new <- mu(sample_X_new)
      p_all <- 1 - pnorm(l - mu_X_new)
      sample_Y_new <- get_Y(sample_X_new)
      sample_Z_new <- (sample_Y_new > l)
      sample_S_new <- S_add2(sample_X_new, sample_X_all, sample_Z_all, h_new_mat, combo, weight_ce, h_s_ce)
      sample_S_all <- c(sample_S_all, sample_S_new)
      sample_H_new <- sqrt(sample_S_new)
      est_temp_new <- get_Est(sample_Y_new, sample_H_new, l_target)
      sample_Y_all <- c(sample_Y_all, sample_Y_new)
      mse_all <- c(mse_all, sum((p_all - sample_S_new)^2) / 199)
      est_vec <- c(est_vec, est_temp_new)
      
      # Update sample bank
      sample_X_all <- rbind(sample_X_all, sample_X_new)
      sample_Z_all <- c(sample_Z_all, sample_Z_new)
      h_init_mat <- h_new_mat
    }
    
    # list(est_vec = est_vec, sample_X_all = sample_X_all, sample_Y_all = sample_Y_all,
         # index = index, h_s_pareto = h_s_pareto, weight_ce_all = weight_ce_all,
         # index_all = index_all, bandwidth = bandwidth, sample_S_all = sample_S_all,
         # mse_all = mse_all)
    mean(est_vec)
  }
  
  return(POE_est)
}
