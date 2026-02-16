#' =============================================================================
#' IJDS Paper - Main Runner for Result Table
#' =============================================================================
#' 
#' Unified script to run the 8 core methods across all 4 examples and generate
#' the main result table for the paper.
#' 
#' Methods: IS-VS (proposed), IS-CE, IS-Pareto, WAMK-SIS, Lasso, SpAM, RF-RFE, OptiTreeStrat
#' 
#' Usage:
#'   source("main.R")
#'   
#'   # Run single example with all methods
#'   results <- run_example(1, n_reps = 25, ncores = 25)
#'   
#'   # Run all examples
#'   all_results <- run_all_examples(n_reps = 25, ncores = 25)
#'   
#'   # Generate result table
#'   table <- generate_result_table(all_results)
#' =============================================================================
#' 
#' ============================================
#  Check for required packages
#  ============================================
required_packages <- c("doParallel", "truncnorm", "splines", "rootSolve", "rmutil", 
                       "pracma", "mvtnorm", "cubature", "caret", "stats", "lhs", 
                       "LHD", "R.utils", "crch", "loo", "GauPro", "DiceKriging", "kernlab", 
                       "SAM", "glmnet", "randomForest")

missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if (length(missing_packages) > 0) {
  message("The following packages are missing: ", paste(missing_packages, collapse = ", "))
  message("Please install them using:")
  message('install.packages(c("', paste(missing_packages, collapse = '", "'), '"))')
  stop("Missing packages. Please install and re-run.")
}

library(doParallel)

# =============================================================================
# Configuration
# =============================================================================

ROOT <- getwd()
RESULTS_DIR <- file.path(ROOT, "results")

if (!dir.exists(RESULTS_DIR)) {
  dir.create(RESULTS_DIR, recursive = TRUE)
}

# Core methods for the paper (8 methods)
CORE_METHODS <- c("IS-VS", "IS-CE", "IS-Pareto", "WAMK-SIS", "Lasso", "SpAM", "RF-RFE", "OptiTree")

# Example parameters
EXAMPLE_PARAMS <- list(
  Ex1 = list(D = 4,  l_target = 18.99, n_pairs = 6,    seed_base = 123,  functions = "functions_ex1.R"),
  Ex2 = list(D = 5,  l_target = 24.2,  n_pairs = 10,   seed_base = 1235, functions = "functions_ex2.R"),
  Ex3 = list(D = 10, l_target = 18.74, n_pairs = 45,   seed_base = 123,  functions = "functions_ex3.R"),
  Ex4 = list(D = 50, l_target = 18.42, n_pairs = 1225, seed_base = 127,  functions = "functions_ex4.R")
)

# Method availability by example 
METHOD_AVAILABILITY <- list(
  Ex1 = c("IS-VS", "IS-CE", "IS-Pareto", "WAMK-SIS", "Lasso", "SpAM", "RF-RFE", "OptiTree", "RF_VS", "GP_VS"),
  Ex2 = c("IS-VS", "IS-CE", "IS-Pareto", "WAMK-SIS", "Lasso", "SpAM", "RF-RFE", "OptiTree", "GP_VS"),
  Ex3 = c("IS-VS", "IS-CE", "IS-Pareto", "WAMK-SIS", "Lasso", "SpAM", "RF-RFE", "OptiTree", "RF_VS", "GP_VS"), 
  Ex4 = c("IS-VS", "IS-CE", "IS-Pareto", "WAMK-SIS", "Lasso", "SpAM", "RF-RFE", "OptiTree", "RF_VS", "GP_VS")  
)

# Method file mapping
METHOD_FILES <- list(
  "IS-VS"     = "IS_VS.R",
  "IS-CE"     = "IS_CE.R",
  "IS-Pareto" = "IS_Pareto.R",
  "WAMK-SIS"  = "wamk.R",
  "Lasso"     = "Lasso.R",
  "SpAM"      = "SpAM.R",
  "RF-RFE"    = "randomForest.R",
  "OptiTree"  = "OptiTreeStrata.R",
  "RF_VS" = "randomForest_VS.R", 
  "GP_VS" = "GP_VS.R"
)

METHOD_FUNCTIONS <- list(
  "IS-VS"     = "run_is_vs",
  "IS-CE"     = "run_is_ce",
  "IS-Pareto" = "run_is_pareto",
  "WAMK-SIS"  = "run_wamk_sis",
  "Lasso"     = "run_lasso",
  "SpAM"      = "run_spam",
  "RF-RFE"    = "run_rf_rfe",
  "OptiTree"  = "run_optitree",
  "RF_VS" = "run_rf_vs",
  "GP_VS" = "run_gp_vs"
)

# =============================================================================
# Core Functions
# =============================================================================

#' Run a single method for a single example
#' 
#' @param example Integer 1-4 specifying which example
#' @param method Character string specifying which method
#' @param n_reps Number of replications
#' @param ncores Number of parallel cores
#' @param save_results Whether to save results to file
#' @return List with results
run_method <- function(example, 
                       method, 
                       n_reps = 25, 
                       ncores = 25, 
                       save_results = TRUE) {
  
  ex_name <- paste0("Ex", example)
  
  # Check method availability
if (!method %in% METHOD_AVAILABILITY[[ex_name]]) {
    warning(sprintf("Method '%s' not available for %s. Skipping.", method, ex_name))
    return(list(
      example = example,
      method = method,
      available = FALSE,
      error = "Method not available for this example"
    ))
  }
  
  message(sprintf("\n=== Running %s: %s ===", ex_name, method))
  message(sprintf("Replications: %d, Cores: %d", n_reps, ncores))
  
  # Construct method file path
  method_file <- file.path(ROOT, ex_name, "methods", 
                           paste0("ex", example, "_", METHOD_FILES[[method]]))
  
  if (!file.exists(method_file)) {
    warning(sprintf("Method file not found: %s", method_file))
    return(list(
      example = example,
      method = method,
      available = FALSE,
      error = paste("File not found:", method_file)
    ))
  }
  
  # Source and run
  source(method_file)
  run_fn <- get(METHOD_FUNCTIONS[[method]])
  
  start_time <- Sys.time()
  
  tryCatch({
    results <- run_fn(n_reps = n_reps, ncores = ncores, ROOT = ROOT)
    end_time <- Sys.time()
    elapsed <- end_time - start_time
    
    message(sprintf("Completed in %s", format(elapsed)))
    
    output <- list(
      example = example,
      method = method,
      available = TRUE,
      POE_est = results,
      n_reps = n_reps,
      elapsed_time = elapsed,
      timestamp = Sys.time()
    )
    
    if (save_results) {
      output_file <- file.path(RESULTS_DIR, 
                               sprintf("ex%d_%s.Rdata", example, gsub("-", "_", method)))
      save(output, file = output_file)
      message(sprintf("Saved to: %s", output_file))
    }
    
    return(output)
    
  }, error = function(e) {
    warning(sprintf("Method %s failed for %s: %s", method, ex_name, e$message))
    return(list(
      example = example,
      method = method,
      available = TRUE,
      error = e$message
    ))
  })
}


#' Run all methods for a single example
#' 
#' @param example Integer 1-4 specifying which example
#' @param methods Vector of method names (default: all core methods)
#' @param n_reps Number of replications
#' @param ncores Number of parallel cores
#' @param save_results Whether to save results to file
#' @return List with results for each method
run_example <- function(example, 
                        methods = CORE_METHODS, 
                        n_reps = 25, 
                        ncores = 25,
                        save_results = TRUE) {
  
  ex_name <- paste0("Ex", example)
  message(sprintf("\n########## Running Example %d (D=%d) ##########\n", 
                  example, EXAMPLE_PARAMS[[ex_name]]$D))
  
  results <- list()
  
  for (method in methods) {
    results[[method]] <- run_method(example, method, n_reps, ncores, save_results)
  }
  
  return(results)
}


#' Run all examples with all methods
#' 
#' @param examples Vector of example numbers (default: 1:4)
#' @param methods Vector of method names (default: all core methods)
#' @param n_reps Number of replications
#' @param ncores Number of parallel cores
#' @param save_results Whether to save results to file
#' @return Nested list with results for each example and method
run_all_examples <- function(examples = 1:4,
                             methods = CORE_METHODS,
                             n_reps = 25,
                             ncores = 25,
                             save_results = TRUE) {
  
  message("==========================================================")
  message("  IJDS Paper - Running All Examples")
  message(sprintf("  Examples: %s", paste(examples, collapse = ", ")))
  message(sprintf("  Methods: %s", paste(methods, collapse = ", ")))
  message(sprintf("  Replications: %d, Cores: %d", n_reps, ncores))
  message("==========================================================\n")
  
  all_results <- list()
  total_start <- Sys.time()
  
  for (ex in examples) {
    ex_name <- paste0("Ex", ex)
    all_results[[ex_name]] <- run_example(ex, methods, n_reps, ncores, save_results)
  }
  
  total_end <- Sys.time()
  total_elapsed <- total_end - total_start
  
  message("\n==========================================================")
  message(sprintf("  Total time: %s", format(total_elapsed)))
  message("==========================================================")
  
  # Save combined results
  if (save_results) {
    combined_file <- file.path(RESULTS_DIR, "all_results.Rdata")
    save(all_results, file = combined_file)
    message(sprintf("Combined results saved to: %s", combined_file))
  }
  
  return(all_results)
}


#' Extract probability estimates from results
#' 
#' @param result Single method result object
#' @return Vector of probability estimates
extract_estimates <- function(result) {
  if (!is.null(result$error) || !result$available) {
    return(NA)
  }
  
  poe <- result$POE_est
  
  # Handle different return formats
  if (is.list(poe)) {
    # Extract est_vec from each replication
    estimates <- sapply(poe, function(x) {
      if (is.list(x) && "est_vec" %in% names(x)) {
        mean(x$est_vec, na.rm = TRUE)
      } else if (is.list(x)) {
        # Try first element (often est_vec)
        if (length(x) > 0 && is.numeric(x[[1]])) {
          mean(x[[1]], na.rm = TRUE)
        } else {
          NA
        }
      } else if (is.numeric(x)) {
        mean(x, na.rm = TRUE)
      } else {
        NA
      }
    })
    return(estimates)
  } else if (is.numeric(poe)) {
    return(poe)
  } else {
    return(NA)
  }
}


#' Generate result table from all results
#' 
#' @param all_results Results from run_all_examples()
#' @param metric Which metric to compute: "mean", "sd", "rmse", or "all"
#' @return Data frame with results table
generate_result_table <- function(all_results, metric = "all") {
  
  examples <- names(all_results)
  methods <- CORE_METHODS
  
  # Initialize result matrices
  mean_mat <- matrix(NA, nrow = length(methods), ncol = length(examples),
                     dimnames = list(methods, examples))
  sd_mat <- matrix(NA, nrow = length(methods), ncol = length(examples),
                   dimnames = list(methods, examples))
  
  for (ex in examples) {
    for (method in methods) {
      if (!is.null(all_results[[ex]][[method]])) {
        estimates <- extract_estimates(all_results[[ex]][[method]])
        if (!all(is.na(estimates))) {
          mean_mat[method, ex] <- mean(estimates, na.rm = TRUE)
          sd_mat[method, ex] <- sd(estimates, na.rm = TRUE)
        }
      }
    }
  }
  
  if (metric == "mean") {
    return(as.data.frame(mean_mat))
  } else if (metric == "sd") {
    return(as.data.frame(sd_mat))
  } else {
    # Return combined table with mean (sd) format
    result_df <- data.frame(Method = methods)
    
    for (ex in examples) {
      result_df[[ex]] <- sprintf("%.4f (%.4f)", mean_mat[, ex], sd_mat[, ex])
    }
    
    return(result_df)
  }
}


#' Print formatted result table
#' 
#' @param result_table Data frame from generate_result_table()
print_result_table <- function(result_table) {
  message("\n==========================================================")
  message("  IJDS Paper - Main Result Table")
  message("  Format: Mean (SD)")
  message("==========================================================\n")
  
  print(result_table, row.names = FALSE)
  
  message("\n")
}


#' Load saved results
#' 
#' @param results_dir Directory containing saved results
#' @return Nested list with all results
load_results <- function(results_dir = RESULTS_DIR) {
  
  all_results <- list()
  
  for (ex in 1:4) {
    ex_name <- paste0("Ex", ex)
    all_results[[ex_name]] <- list()
    
    for (method in CORE_METHODS) {
      file_name <- sprintf("ex%d_%s.Rdata", ex, gsub("-", "_", method))
      file_path <- file.path(results_dir, file_name)
      
      if (file.exists(file_path)) {
        load(file_path)  # loads 'output'
        all_results[[ex_name]][[method]] <- output
      }
    }
  }
  
  return(all_results)
}


# =============================================================================
# Quick Start Examples
# =============================================================================

message("
==========================================================
  IJDS Paper - Main Runner Loaded
==========================================================

Quick Start:

  # Run a single method for one example
  result <- run_method(1, 'IS-VS', n_reps = 25, ncores = 25)

  # Run all methods for one example
  ex1_results <- run_example(1, n_reps = 25, ncores = 25)

  # Run everything (all 4 examples, all 8 methods)
  all_results <- run_all_examples(n_reps = 25, ncores = 25)

  # Generate and print result table
  table <- generate_result_table(all_results)
  print_result_table(table)

  # Load previously saved results
  all_results <- load_results()

==========================================================
")


# =============================================================================
# Command Line Interface
# =============================================================================

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  
  # Default parameters
  examples <- 1:4
  methods <- CORE_METHODS
  n_reps <- 25
  ncores <- 25
  
  # Parse arguments
  for (arg in args) {
    if (grepl("^--example=", arg)) {
      examples <- as.integer(strsplit(sub("^--example=", "", arg), ",")[[1]])
    } else if (grepl("^--method=", arg)) {
      methods <- strsplit(sub("^--method=", "", arg), ",")[[1]]
    } else if (grepl("^--reps=", arg)) {
      n_reps <- as.integer(sub("^--reps=", "", arg))
    } else if (grepl("^--cores=", arg)) {
      ncores <- as.integer(sub("^--cores=", "", arg))
    }
  }
  
  # Run
  all_results <- run_all_examples(examples, methods, n_reps, ncores)
  
  # Generate and print table
  table <- generate_result_table(all_results)
  print_result_table(table)
}
