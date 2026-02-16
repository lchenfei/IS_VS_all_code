#' OptiTreeStrata Method - Example 1
#' 
#' Wrapper for OptiTreeStrata stratified sampling method (D=4).
#' Based on case='2' in the original OptiTreeStrata code.
#'
#' @param n_reps Number of independent replications (default: 25)
#' @param ncores Number of parallel cores (default: 25)
#' @param seed_base Base seed for reproducibility (default: 123)
#' @param n_total Total computational budget (default: 1200)
#' @param n_pilot Pilot sample size (default: 200)
#' @param m Number of batches (default: 5)
#' @param ROOT Root directory for sourcing shared functions
#' @return List of results for each replication

library(doParallel)

run_optitree <- function(n_reps = 25, 
                         ncores = 25, 
                         seed_base = 123,
                         n_total = 1200,
                         n_pilot = 200,
                         m = 5,
                         ROOT = getwd()) {
  
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl), add = TRUE)
  clusterExport(cl, "ROOT", envir = environment())
  
  POE_est <- foreach(j = 1:n_reps, 
                     .combine = c, 
                     .errorhandling = 'remove',
                     .packages = c("truncnorm", "splines", "rootSolve", "rmutil", 
                                   "pracma", "mvtnorm", "cubature", "caret", 
                                   "stats", "lhs", "LHD", "R.utils", "crch")) %dopar% {
    
    # Set global variables required by OptiTreeStrata
    Interest <<- 'Failure_prob'
    case <<- '2'  # Ex1 corresponds to case='2' in OptiTreeStrata
    
    # Source OptiTreeStrata functions
    source(file.path(ROOT, "AMISE_univariate_OptiTree.R"))
    source(file.path(ROOT, "function_OptiTreeStrata_independent.R"))
    
    set.seed(seed_base + j)
    print(paste0('OptiTreeStrata Ex1, experiment:', j))
    
    # Parameters for Example 1 (case='2')
    p <<- 4
    l <<- 18.99
    max_x <<- 5
    min_x <<- -5
    
    I <- 2       # Binary tree
    l_RR <- 0.05 # threshold of reduction rate
    nj_min <- 10 # minimal allocation within stratum
    
    n_simulation <- n_total - n_pilot
    
    # Draw pilot samples using OA-LHS
    Pilot_sample <- Draw_Pilot_sample_OA_LHS(n_pilot, p, min_x, max_x)
    Pilot_sample$Z <- Make_Z(Pilot_sample$Y, l)
    # Pilot_sample <- Pilot_sample[1:n_pilot, ]
    
    # Handle case when no exceeding samples in pilot
    Num_taking_budget <- 0
    while ((mean(Pilot_sample$Z) == 0) & (Num_taking_budget < m)) {
      Additional_Pilot <- Draw_Pilot_sample(n_pilot, p, min_x, max_x)
      Additional_Pilot$Z <- Make_Z(Additional_Pilot$Y, l)
      Pilot_sample <- rbind(Pilot_sample, Additional_Pilot)
      Num_taking_budget <- Num_taking_budget + 1
    }
    
    # Find optimal bandwidth
    h_hat <<- find_bandwidth(Pilot_sample)
    
    # Fit proposed model
    Initial_stratification <- Tree_stratification(Pilot_sample, n_simulation, nj_min, l, I, m, l_RR)
    Simulation_result <- Simulation_MSS(Pilot_sample, Initial_stratification, l, n_simulation, 
                                        nj_min, m, I, l_RR, j, m - Num_taking_budget)
    
    Ehat_MSS <- Simulation_result$Ehat_MSS
    
    Ehat_MSS
  }
  
  return(POE_est)
}
