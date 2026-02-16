#' RF-RFE Method for Importance Sampling Variable Selection - Example 2
#' 
#' Uses Random Forest with Recursive Feature Elimination for variable selection (D=10).
#'
#' @param n_reps Number of independent replications (default: 25)
#' @param ncores Number of parallel cores (default: 25)
#' @param seed_base Base seed for reproducibility (default: 123)
#' @param ROOT Root directory for sourcing shared functions
#' @return List of results for each replication

library(doParallel)

run_rf_rfe <- function(n_reps = 25, 
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
                                   "loo", "GauPro", "caret", "randomForest")) %dopar% {
    
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
    sample_Y_all <- sample_Y
    h_s_pareto_all <- c()
    sample_S_all <- c()
    
    for (simu in 1:n_iter) {
      sample_X_all <- as.data.frame(sample_X_all)
      z_factor <- factor(sample_Z_all, levels = c(FALSE, TRUE))
      rf_model <- randomForest(z_factor ~ ., data = sample_X_all, 
                               importance = TRUE, ntree = 500)
      
      CE_summary <- function(data, lev = NULL, model = NULL) {
        epsilon <- 1e-15
        obs <- ifelse(data$obs == TRUE, 1, 0)
        probs <- pmin(pmax(data[, "TRUE"], epsilon), 1 - epsilon)
        ce <- -mean(obs * log(probs) + (1 - obs) * log(1 - probs))
        out <- c(CE = ce)
        names(out) <- "CE"
        return(out)
      }
      
      custom_rfFuncs <- rfFuncs
      custom_rfFuncs$summary <- CE_summary
      custom_rfFuncs$prob <- function(object, x) {
        predict(object, x, type = "prob")
      }
      
      ctrl <- rfeControl(functions = custom_rfFuncs,
                         method = "cv", 
                         allowParallel = FALSE,
                         number = 5)
      
      sizes <- 1:10
      results <- rfe(x = sample_X_all, y = z_factor,
                     sizes = sizes,
                     metric = "CE",
                     rfeControl = ctrl, 
                     maximize = FALSE)
      
      selected_vars <- results$optVariables
      h_s_pareto_all <- c(h_s_pareto_all, list(selected_vars))
      
      X_selected <- sample_X_all[, selected_vars]
      
      final_rf_model <- randomForest(x = X_selected,
                                     y = z_factor,
                                     ntree = 500,
                                     importance = TRUE)
      
      # Rejection sampling (iid normal for Ex2)
      x <- c()
      k_all <- c()
      i <- 0
      
      while ((nrow(x) < Nm) || is.null(x)) {
        i <- i + 1
        x_new <- rmvnorm(1, rep(0, D))
        u <- runif(1, min = 0, max = dmvnorm(x_new))
        x_new <- as.data.frame(x_new)
        x_new_selected <- x_new[, selected_vars]
        k <- predict(final_rf_model, newdata = x_new_selected, type = "prob")[, "TRUE"]
        if (u <= sqrt(k) * dmvnorm(x_new)) {
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
      
      sample_X_all <- rbind(sample_X_all, sample_X_new)
      sample_Z_all <- c(sample_Z_all, sample_Z_new)
    }
    
    # list(est_vec = est_vec, sample_X_all = sample_X_all, sample_Y_all = sample_Y_all,
    #     sample_S_all = sample_S_all, index = index, h_s_pareto_all = h_s_pareto_all)
    mean(est_vec)
  }
  
  return(POE_est)
}
