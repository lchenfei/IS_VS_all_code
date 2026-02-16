#' GP-VS Method for Importance Sampling Variable Selection - Example 1
#' 
#' Uses Gaussian Process with enumeration of all subsets and validation-based selection (D=4).
#' Enumerates all 2^D - 1 = 15 variable subsets, applies Pareto-k filtering.
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
                                   "loo", "GauPro", "DiceKriging")) %dopar% {
    
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
    
    est_vec <- c()
    index <- c()
    sample_Y_all <- sample_Y
    h_s_pareto_all <- c()
    sample_S_all = c()
    mse_all = c()
    
    for (simu in 1:n_iter) {
      # Enumerate all 2^D - 1 = 15 variable subsets
      all_combinations <- unlist(lapply(1:D, function(k) combn(1:D, k, simplify = FALSE)), recursive = FALSE)
      
      S_all = c()
      cost_vec = c()
      filter = c()
      for (c in 1:length(all_combinations)){
        cb = all_combinations[[c]]
        X_train = as.matrix(sample_X_all[, cb])
        gp_model = km(formula = ~1,     # Constant mean
                      design = X_train, # Input variables
                      response = sample_Y_all, # Output variable
                      covtype = "gauss", nugget.estim = TRUE)  # Kernel type
        
        S_fit <- predict(gp_model, X_train, type = "UK",                # Universal Kriging
                         se.compute = TRUE )          # Compute standard error
        
        
        mu_fit = S_fit$mean
        sd_fit = S_fit$sd
        
        S_temp = 1 - pnorm((l - mu_fit)/sd_fit)
        
        S_all = c(S_all, S_temp)
        cost = -sum(sample_Z_all*log(S_temp)+(1-sample_Z_all)*log(1-S_temp), na.rm = T)
        cost_vec = c(cost_vec, cost)
        filter = c(filter, c)
        
      }
      
      order_s = order(cost_vec)
      
      h_s_pareto = c()
      sample_x_prepare_all = c()
      sample_s_prepare_all = c()
      
      for (c in order_s){
        cb = all_combinations[[c]]
        X_train = sample_X_all[, cb]
        
        sample_new = sampling_X_new_GP(Nm,sample_X_all,sample_Y_all, l, cb)
        
        sample_X_new = sample_new[[1]]
        sample_S_new = sample_new[[2]]
        sample_w = (1/sqrt(sample_S_new))/(sum(1/sqrt(sample_S_new)))
        diagnostic_pareto = psis(sample_w, r_eff = 1)
        k_score = diagnostic_pareto$diagnostics$pareto_k
        
        print(k_score)
        
        sample_x_prepare_all = c(sample_x_prepare_all, list(sample_X_new))
        sample_s_prepare_all = c(sample_s_prepare_all, list(sample_S_new))
        if(k_score < 0.7){
          h_s_pareto = c(h_s_pareto, list(cb))
        }
        if(length(h_s_pareto) >= 8){
          break    
          
        }
      }
      
      
      
      h_s_all_select = S_find_index_validation_GP(sample_X_all, sample_Y_all, h_s_pareto, l, Nm)
      h_s_ce =  h_s_all_select[[1]]
      h_track = h_s_all_select[[2]]
      
      h_s_all = all_combinations[order_s]
      
      index_pareto_list = sapply(h_s_all, function(x) setequal(x, h_s_ce))
      index_pareto = which(index_pareto_list == 1)
      sample_X_new = sample_x_prepare_all[[index_pareto]]
      
      index = c(index, list(h_s_ce))
      # index_all = c(index_all, list(h_track))
      h_s_pareto_all = c(h_s_pareto_all, list(h_s_pareto))
      
      mu_X_new = mu(sample_X_new)
      sample_Y_new = get_Y(sample_X_new)
      sample_Z_new = (sample_Y_new>l)
      sample_S_new = sample_s_prepare_all[[index_pareto]]
      sample_S_all = c(sample_S_all, sample_S_new)
      p_all = 1 - pnorm(l - mu_X_new)
      mse_all = c(mse_all, sum((p_all - mu_X_new) ** 2))
      sample_H_new = sqrt(sample_S_new)
      est_temp_new = get_Est(sample_Y_new,sample_H_new,l_target)
      sample_Y_all = c(sample_Y_all, sample_Y_new)
      est_vec = c(est_vec,est_temp_new)
      # update sample bank
      sample_X_all = rbind(sample_X_all,sample_X_new)
      sample_Z_all = c(sample_Z_all,sample_Z_new)
    }
    
    # list(est_vec = est_vec, sample_X_all = sample_X_all, sample_Y_all = sample_Y_all,
    #      index = index, h_s_pareto_all = h_s_pareto_all)
    mean(est_vec)
  }
  
  return(POE_est)
}
