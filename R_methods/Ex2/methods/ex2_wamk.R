#' WAMK-SIS Method for Importance Sampling - Example 2
#' 
#' Uses Weighted Additive Multi-Kernel approach with all kernel subsets.
#'
#' @param n_reps Number of independent replications (default: 25)
#' @param ncores Number of parallel cores (default: 25)
#' @param seed_base Base seed for reproducibility (default: 1234)
#' @param ROOT Root directory for sourcing shared functions
#' @param save_details If TRUE, returns detailed results
#' @return List of results for each replication

library(doParallel)

run_wamk_sis <- function(n_reps = 25, 
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
                     .packages = c("pracma", "mvtnorm", "cubature", "lhs")) %dopar% {
    
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
    
    # Get starting points of bandwidth
    h_init_vec <- c()
    for (i in 1:D) {
      h <- (optimize(function(x) AMISE(sample_X[, i], sample_Z, x, -5, 5), 
                     interval = c(0.1, 1.5)))$minimum
      h_init_vec[i] <- h
    }
    h_init_mat <- matrix(h_init_vec[combo], ncol = 2)
    
    est_vec <- c()
    index <- c()
    bandwidth <- c()
    weight_s_all <- c()
    sample_Y_all <- sample_Y
    
    for (simu in 1:n_iter) {
      # Choose optimal bandwidth
      h_new_mat <- 0.8 * h_init_mat
      for (c in 1:nrow(combo)) {
        cb <- combo[c, ]
        h_init <- h_init_mat[c, ]
        h <- try(get_h_V2_final(sample_X_all[, cb], sample_Z_all, h_init[1], h_init[2]), silent = TRUE)
        if (!"try-error" %in% class(h)) {
          h_new <- h[nrow(h), ]
          h_new_mat[c, ] <- h_new
        }
      }
      bandwidth <- c(bandwidth, list(h_new_mat))
      
      # Deciding weights
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
      S_all <- matrix(S_all, nrow = length(S_all) / 10)  # 10 = C(5,2)
      weight <- 1 / cost_vec / sum(1 / cost_vec)
      
      # Subset is all kernels
      subset <- 1:10
      
      weight_s <- weight / sum(weight[subset])
      index <- c(index, list(subset))
      weight_s_all <- c(weight_s_all, list(weight_s))
      
      # Get new samples for the next round
      sample_X_new <- sampling_X_new(Nm, sample_X_all, sample_Z_all, h_new_mat, combo, weight_s, subset)
      sample_Y_new <- get_Y(sample_X_new)
      sample_Z_new <- (sample_Y_new > l)
      sample_S_new <- S_add2(sample_X_new, sample_X_all, sample_Z_all, h_new_mat, combo, weight_s, subset)
      sample_H_new <- sqrt(sample_S_new)
      est_temp_new <- get_Est(sample_Y_new, sample_H_new, l_target)
      sample_Y_all <- c(sample_Y_all, sample_Y_new)
      est_vec <- c(est_vec, est_temp_new)
      
      # Update sample bank
      sample_X_all <- rbind(sample_X_all, sample_X_new)
      sample_Z_all <- c(sample_Z_all, sample_Z_new)
      h_init_mat <- h_new_mat
    }
    
    # Return results
    # list(est_vec = est_vec, 
    #      sample_X_all = sample_X_all, 
    #      sample_Y_all = sample_Y_all, 
    #      index = index, 
    #      weight_s_all = weight_s_all, 
    #      bandwidth = bandwidth)
    mean(est_vec)
  }
  
  return(POE_est)
}
