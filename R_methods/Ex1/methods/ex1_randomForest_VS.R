#' RF-VS Method for Importance Sampling Variable Selection - Example 1
#' 
#' Uses Random Forest with greedy CE ordering and validation-based subset selection (D=4).
#' Enumerates all 2^D - 1 = 15 variable subsets, orders by CE, applies Pareto-k filtering.
#'
#' @param n_reps Number of independent replications (default: 25)
#' @param ncores Number of parallel cores (default: 25)
#' @param seed_base Base seed for reproducibility (default: 123)
#' @param ROOT Root directory for sourcing shared functions
#' @return List of results for each replication

library(doParallel)

run_rf_vs <- function(n_reps = 25, 
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
                                   "loo", "GauPro", "DiceKriging", "kernlab", "randomForest")) %dopar% {
    
    source(file.path(ROOT, "AMISE_bivariate.R"))
    source(file.path(ROOT, "functions_ex1.R"))
    source(file.path(ROOT, "Additive_S.R"))
    source(file.path(ROOT, "cv.R"))
    
    set.seed(seed_base + j)
    print(paste0('experiment:', j))
    
    D <- 4
    combo <- t(combn(seq(1:D), 2))
    Nm <- 200
    Nm_init <- 200
    l_target <- 18.99
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
    sample_Y_all <- sample_Y
    h_s_pareto_all <- c()
    sample_S_all <- c()
    
    for (simu in 1:n_iter) {
      # Enumerate all 2^D - 1 = 15 variable subsets
      all_combinations <- unlist(lapply(1:D, function(k) combn(1:D, k, simplify = FALSE)), recursive = FALSE)
      
      S_all <- c()
      cost_vec <- c()
      filter <- c()
      
      # Compute CE for each subset
      for (c in 1:length(all_combinations)) {
        cb <- all_combinations[[c]]
        X_train <- sample_X_all[, cb]
        X_train <- as.data.frame(X_train)
        z_factor <- factor(sample_Z_all, levels = c(FALSE, TRUE))
        rf_model <- randomForest(z_factor ~ ., data = X_train, 
                                 importance = TRUE, ntree = 500)
        S_temp <- predict(rf_model, X_train, type = "prob")[, "TRUE"]
        
        if (any(S_temp < 0) || any(S_temp > 1)) {
          cat("Skipping combination", c, "due to invalid probabilities.\n")
          next
        }
        
        S_all <- c(S_all, S_temp)
        cost <- -sum(sample_Z_all * log(S_temp) + (1 - sample_Z_all) * log(1 - S_temp), na.rm = TRUE)
        cost_vec <- c(cost_vec, cost)
        filter <- c(filter, c)
      }
      
      order_s <- order(cost_vec)
      
      # Pareto-k filtering
      h_s_pareto <- c()
      sample_x_prepare_all <- c()
      sample_s_prepare_all <- c()
      
      for (c in order_s) {
        filter_count <- filter[c]
        cb <- all_combinations[[filter_count]]
        sample_new <- sampling_X_new_randomForest(Nm, sample_X_all, sample_Z_all, l, cb)
        sample_X_new <- sample_new[[1]]
        sample_S_new <- sample_new[[2]]
        sample_w <- (1 / sqrt(sample_S_new)) / (sum(1 / sqrt(sample_S_new)))
        diagnostic_pareto <- psis(sample_w, r_eff = 1)
        k_score <- diagnostic_pareto$diagnostics$pareto_k
        print(k_score)
        sample_x_prepare_all <- c(sample_x_prepare_all, list(sample_X_new))
        sample_s_prepare_all <- c(sample_s_prepare_all, list(sample_S_new))
        if (k_score < 0.7) {
          h_s_pareto <- c(h_s_pareto, list(cb))
        }
        if (length(h_s_pareto) >= 8) {
          break
        }
      }
      
      # Validation-based subset selection
      h_s_all_select <- S_find_index_validation_randomForest(sample_X_all, sample_Z_all, h_s_pareto, l, Nm)
      h_s_ce <- h_s_all_select[[1]]
      h_track <- h_s_all_select[[2]]
      
      h_s_all <- all_combinations[order_s]
      
      index_pareto_list <- sapply(h_s_all, function(x) setequal(x, h_s_ce))
      index_pareto <- which(index_pareto_list == 1)
      sample_X_new <- sample_x_prepare_all[[index_pareto]]
      
      index <- c(index, list(h_s_ce))
      h_s_pareto_all <- c(h_s_pareto_all, list(h_s_pareto))
      
      sample_Y_new <- get_Y(sample_X_new)
      sample_Z_new <- (sample_Y_new > l)
      sample_S_new <- sample_s_prepare_all[[index_pareto]]
      sample_S_all <- c(sample_S_all, sample_S_new)
      sample_H_new <- sqrt(sample_S_new)
      est_temp_new <- get_Est(sample_Y_new, sample_H_new, l_target)
      sample_Y_all <- c(sample_Y_all, sample_Y_new)
      est_vec <- c(est_vec, est_temp_new)
      
      sample_X_all <- rbind(sample_X_all, sample_X_new)
      sample_Z_all <- c(sample_Z_all, sample_Z_new)
    }
    
    # list(est_vec = est_vec, sample_X_all = sample_X_all, sample_Y_all = sample_Y_all,
    #     sample_S_all = sample_S_all, index = index, h_s_pareto_all = h_s_pareto_all)
    mean(est_vec)
  }
  
  return(POE_est)
}
