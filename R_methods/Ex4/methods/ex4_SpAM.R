#' SpAM Method for Importance Sampling Variable Selection - Example 4
#' 
#' Uses Sparse Additive Models with CV-based lambda selection (D=50).
#'
#' @param n_reps Number of independent replications (default: 25)
#' @param ncores Number of parallel cores (default: 25)
#' @param seed_base Base seed for reproducibility (default: 127)
#' @param ROOT Root directory for sourcing shared functions
#' @return List of results for each replication

library(doParallel)

run_spam <- function(n_reps = 25, 
                     ncores = 25, 
                     seed_base = 127,
                     ROOT = getwd()) {
  
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl), add = TRUE)
  clusterExport(cl, "ROOT", envir = environment())
  
  POE_est <- foreach(j_simulation = 1:n_reps, 
                     .combine = c, 
                     .errorhandling = 'remove',
                     .packages = c("pracma", "mvtnorm", "cubature", "lhs", 
                                   "loo", "GauPro", "DiceKriging", "SAM")) %dopar% {
    
    source(file.path(ROOT, "AMISE_bivariate.R"))
    source(file.path(ROOT, "functions_ex4.R"))
    source(file.path(ROOT, "Additive_S.R"))
    source(file.path(ROOT, "cv.R"))
    
    set.seed(seed_base + j_simulation)
    print(paste0('experiment:', j_simulation))
    
    D <- 50
    combo <- t(combn(seq(1:D), 2))
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
    sample_Y_all <- sample_Y
    h_s_pareto_all <- c()
    sample_S_all <- c()
    
    for (simu in 1:n_iter) {
      K <- 5
      n <- nrow(sample_X_all)
      folds <- sample(rep(1:K, length.out = n))
      
      set.seed(1238 + j_simulation)
      fit <- samLL(as.matrix(sample_X_all), sample_Z_all)
      
      lambda_master <- fit$lambda
      ce_matrix <- matrix(NA, nrow = K, ncol = length(lambda_master))
      
      X <- as.matrix(sample_X_all)
      y <- as.vector(sample_Z_all)
      n <- nrow(X)
      d <- ncol(X)
      
      Xs <- sweep(X, 2, apply(X, 2, min), "-")
      Xs <- sweep(Xs, 2, apply(X, 2, max) - apply(X, 2, min), "/")
      
      p <- 3
      m <- d * p
      Z <- matrix(0, n, m)
      for (jj in 1:d)
        Z[, ((jj - 1) * p + 1):(jj * p)] <- ns(Xs[, jj], df = p)
      
      Z1 <- cbind(Z, 1)
      n1 <- sum(y == 1)
      g <- -colSums(y * Z1) + (n1 / n) * colSums(Z1)
      g_mat <- matrix(g[1:(p * d)], p, d)
      lambda_max <- max(sqrt(colSums(g_mat^2)))
      
      lambda_grid <- exp(seq(log(1), log(0.1), length = 20)) * lambda_max
      
      for (k in 1:K) {
        train_idx <- which(folds != k)
        test_idx <- which(folds == k)
        
        X_train <- sample_X_all[train_idx, ]
        Z_train <- sample_Z_all[train_idx]
        X_test <- sample_X_all[test_idx, ]
        Z_test <- sample_Z_all[test_idx]
        
        fit_k <- samLL(as.matrix(X_train), Z_train, lambda = lambda_grid)
        pred_probs_k <- predict(fit_k, newdata = as.matrix(X_test), type = "response")$probs
        
        epsilon <- 1e-15
        pred_probs_k <- pmin(pmax(pred_probs_k, epsilon), 1 - epsilon)
        
        for (jj in 1:ncol(pred_probs_k)) {
          ce_matrix[k, jj] <- -sum(Z_test * log(pred_probs_k[, jj]) +
                                     (1 - Z_test) * log(1 - pred_probs_k[, jj]))
        }
      }
      
      ce_mean <- colMeans(ce_matrix)
      best_j <- which.min(ce_mean)
      
      h_s_pareto_all <- c(h_s_pareto_all, fit$func_norm[, best_j])
      
      # Rejection sampling with error handling
      x <- c()
      k_all <- c()
      i <- 0
      
      while ((nrow(x) < Nm) || is.null(x)) {
        i <- i + 1
        x_new <- rmvnorm(1, rep(0, D))
        u <- runif(1, min = 0, max = dmvnorm(x_new))
        
        pred_new_all <- tryCatch({
          predict(fit, newdata = x_new, type = "response")$probs[, best_j]
        }, error = function(e) { NA })
        
        if (is.na(pred_new_all)) next
        k <- pred_new_all
        
        if (u <= sqrt(k) * dmvnorm(x_new)) {
          x <- rbind(x, x_new)
          k_all <- c(k_all, k)
        }
      }
      
      sample_X_new <- x
      sample_Y_new <- get_Y(sample_X_new)
      sample_Z_new <- (sample_Y_new > l)
      sample_S_new <- k_all
      sample_H_new <- sqrt(sample_S_new)
      est_temp_new <- get_Est(sample_Y_new, sample_H_new, l_target)
      sample_Y_all <- c(sample_Y_all, sample_Y_new)
      est_vec <- c(est_vec, est_temp_new)
      
      sample_X_all <- rbind(sample_X_all, sample_X_new)
      sample_Z_all <- c(sample_Z_all, sample_Z_new)
      sample_S_all <- c(sample_S_all, sample_S_new)
    }
    
    # list(est_vec = est_vec, sample_X_all = sample_X_all, sample_Y_all = sample_Y_all,
    #     sample_S_all = sample_S_all, h_s_pareto_all = h_s_pareto_all)
    mean(est_vec)
  }
  
  return(POE_est)
}
