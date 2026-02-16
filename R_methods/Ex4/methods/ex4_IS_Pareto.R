#' IS-Pareto Method for Importance Sampling Variable Selection - Example 4
#' 
#' Uses greedy CE ordering with Pareto-k diagnostic for subset selection (D=50).
#' Takes the first subset passing Pareto-k threshold (k=1).
#'
#' @param n_reps Number of independent replications (default: 25)
#' @param ncores Number of parallel cores (default: 25)
#' @param seed_base Base seed for reproducibility (default: 127)
#' @param ROOT Root directory for sourcing shared functions
#' @return List of results for each replication

library(doParallel)

run_is_pareto <- function(n_reps = 25, 
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
    n_pairs <- 25 * 49
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
    bandwidth <- c()
    sample_Y_all <- sample_Y
    
    for (simu in 1:n_iter) {
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
      
      # Greedy CE ordering + Pareto-k filter (take first passing subset, k=1)
      h_s_all_sample_x <- sampling_ce_order_greedy(Nm, sample_X_all, sample_Z_all, 
                                                    h_new_mat, combo, weight, S_all, 1)
      h_s_all <- h_s_all_sample_x[[1]]
      sample_x_pareto_all <- h_s_all_sample_x[[2]]
      
      h_s_ce <- h_s_all[[1]]
      sample_X_new <- sample_x_pareto_all[[1]]
      
      weight_ce <- weight / sum(weight[h_s_ce])
      index <- c(index, list(h_s_ce))
      weight_ce_all <- c(weight_ce_all, list(weight_ce))
      
      sample_Y_new <- get_Y(sample_X_new)
      sample_Z_new <- (sample_Y_new > l)
      sample_S_new <- S_add2(sample_X_new, sample_X_all, sample_Z_all, h_new_mat, combo, weight_ce, h_s_ce)
      sample_H_new <- sqrt(sample_S_new)
      est_temp_new <- get_Est(sample_Y_new, sample_H_new, l_target)
      sample_Y_all <- c(sample_Y_all, sample_Y_new)
      est_vec <- c(est_vec, est_temp_new)
      
      sample_X_all <- rbind(sample_X_all, sample_X_new)
      sample_Z_all <- c(sample_Z_all, sample_Z_new)
      h_init_mat <- h_new_mat
    }
    
    # list(est_vec = est_vec, sample_X_all = sample_X_all, sample_Y_all = sample_Y_all,
    #      index = index, weight_ce_all = weight_ce_all, bandwidth = bandwidth)
    mean(est_vec)
  }
  
  return(POE_est)
}
