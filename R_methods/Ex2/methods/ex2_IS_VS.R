#' IS-VS Method for Importance Sampling Variable Selection (Proposed Method) - Example 3
#' 
#' Uses Cross-Entropy ordering with Pareto-k diagnostic filtering and 
#' validation-based subset selection for importance sampling.
#'
#' @param n_reps Number of independent replications (default: 25)
#' @param ncores Number of parallel cores (default: 25)
#' @param seed_base Base seed for reproducibility (default: 1234)
#' @param ROOT Root directory for sourcing shared functions
#' @param save_details If TRUE, returns detailed results
#' @return List of results for each replication

library(doParallel)

run_is_vs <- function(n_reps = 25, 
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
                     .packages = c("pracma", "mvtnorm", "cubature", "lhs", "loo")) %dopar% {
    
    # Source shared functions inside each worker
    source(file.path(ROOT, "AMISE_bivariate.R"))
    source(file.path(ROOT, "functions_ex2.R"))
    source(file.path(ROOT, "Additive_S.R"))
    source(file.path(ROOT, "cv.R"))
    
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
    weight_ce_all <- c()
    sample_Y_all <- sample_Y
    h_s_pareto_all <- c()
    index_all <- c()
    
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
      
      # Order subsets by CE
      h_s_all <- S_find_h_s_ce_order(weight, sample_Z_all, S_all)
      
      # Filter subsets by Pareto-k diagnostic
      h_s_pareto <- c()
      sample_x_prepare_all <- c()
      for (i in 1:length(h_s_all)) {
        h_s_ce <- h_s_all[[i]]
        weight_ce <- weight / sum(weight[h_s_ce])
        sample_X_new <- sampling_X_new(Nm, sample_X_all, sample_Z_all, h_new_mat, combo, weight_ce, h_s_ce)
        sample_x_prepare_all <- c(sample_x_prepare_all, list(sample_X_new))
        sample_S_new <- S_add2(sample_X_new, sample_X_all, sample_Z_all, h_new_mat, combo, weight_ce, h_s_ce)
        sample_w <- (1 / sqrt(sample_S_new)) / (sum(1 / sqrt(sample_S_new)))
        diagnostic_pareto <- psis(sample_w, r_eff = 1)
        k_score <- diagnostic_pareto$diagnostics$pareto_k
        print(k_score)
        if (k_score < 0.7) {
          h_s_pareto <- c(h_s_pareto, list(h_s_ce))
        }
        if (length(h_s_pareto) >= 8) {
          break
        }
      }
      
      # Validation-based subset selection
      h_s_all_select <- S_find_index_validation(sample_X_all, sample_Z_all, h_new_mat, h_s_pareto, weight, Nm)
      h_s_ce <- h_s_all_select[[1]]
      h_track <- h_s_all_select[[2]]
      
      index_pareto_list <- sapply(h_s_all, function(x) setequal(x, h_s_ce))
      index_pareto <- which(index_pareto_list == 1)
      sample_X_new <- sample_x_prepare_all[[index_pareto]]
      
      weight_ce <- weight / sum(weight[h_s_ce])
      index <- c(index, list(h_s_ce))
      index_all <- c(index_all, list(h_track))
      h_s_pareto_all <- c(h_s_pareto_all, list(h_s_pareto))
      weight_ce_all <- c(weight_ce_all, list(weight_ce))
      
      # Compute estimates with new samples
      sample_Y_new <- get_Y(sample_X_new)
      sample_Z_new <- (sample_Y_new > l)
      sample_S_new <- S_add2(sample_X_new, sample_X_all, sample_Z_all, h_new_mat, combo, weight_ce, h_s_ce)
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
    #      h_s_pareto_all = h_s_pareto_all, 
    #      weight_ce_all = weight_ce_all, 
    #      index_all = index_all, 
    #      bandwidth = bandwidth)
    mean(est_vec)
  }
  
  return(POE_est)
}
