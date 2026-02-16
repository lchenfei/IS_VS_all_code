#' IS-Pareto Method for Importance Sampling Variable Selection - Example 2
#' 
#' Uses greedy CE ordering with Pareto-k diagnostic for subset selection (D=10).
#'
#' @param n_reps Number of independent replications (default: 25)
#' @param ncores Number of parallel cores (default: 25)
#' @param seed_base Base seed for reproducibility (default: 123)
#' @param ROOT Root directory for sourcing shared functions
#' @return List of results for each replication

library(doParallel)

run_is_pareto <- function(n_reps = 25, 
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
                     .packages = c("pracma", "mvtnorm", "cubature", "lhs", "loo")) %dopar% {
    
    source(file.path(ROOT, "AMISE_bivariate.R"))
    source(file.path(ROOT, "functions_ex3.R"))
    source(file.path(ROOT, "Additive_S.R"))
    
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
    sample_Y_all <- sample_Y
    
    for (simu in 1:n_iter) {
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
      S_all <- matrix(S_all, nrow = length(S_all) / 45)
      weight <- 1 / cost_vec / sum(1 / cost_vec)
      
      # Greedy CE ordering + Pareto-k filter (take first passing subset)
      h_s_all <- S_find_h_s_ce_order_greedy(weight, sample_Z_all, S_all)
      h_s_pareto <- c()
      for (i in 1:length(h_s_all)) {
        h_s_ce <- h_s_all[[i]]
        weight_ce <- weight / sum(weight[h_s_ce])
        sample_X_new <- sampling_X_new(Nm, sample_X_all, sample_Z_all, h_new_mat, combo, weight_ce, h_s_ce)
        sample_S_new <- S_add2(sample_X_new, sample_X_all, sample_Z_all, h_new_mat, combo, weight_ce, h_s_ce)
        sample_w <- (1 / sqrt(sample_S_new)) / (sum(1 / sqrt(sample_S_new)))
        diagnostic_pareto <- psis(sample_w, r_eff = 1)
        k_score <- diagnostic_pareto$diagnostics$pareto_k
        print(k_score)
        if (k_score < 0.7) {
          h_s_pareto <- c(h_s_pareto, list(h_s_ce))
        }
        if (length(h_s_pareto) >= 1) {
          break
        }
      }
      
      index <- c(index, list(h_s_ce))
      h_s_pareto_all <- c(h_s_pareto_all, list(h_s_pareto))
      weight_ce <- weight / sum(weight[h_s_ce])
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
    #     index = index, h_s_pareto = h_s_pareto, weight_ce_all = weight_ce_all,
    #     bandwidth = bandwidth)
    mean(est_vec)
  }
  
  return(POE_est)
}
