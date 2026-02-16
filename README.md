
Dimension Reduction in Importance Sampling: Balancing Concentration and
Exploration through Variable Selection
==================================================================

This repository contains the code to reproduce the numerical experiments
reported in the paper. 


The code reproduces the results in Table 1.

----------------------------------------------------------------------
A. Reproducing Table 1 and 2 (All Methods Excluding ECL)
----------------------------------------------------------------------

All R-based methods are implemented in the `R_methods` folder.

1. Open RStudio and set working directory to R_methods folder
2. Run the following command to reproduce all R-based results:
    
   source("main.R")
   
   all_results <- run_all_examples(n_reps = 25, ncores = 25) 
   
   table <- generate_result_table(all_results)
   
   print_result_table(table)

    Here, `ncores` denotes the number of CPU cores used for parallel computation.
    Users can set this value according to the number of cores available on their machine.
    
4. The results will appear in the console.

Optional: running subsets of experiments

After sourcing `main.R`, you may run the following commands:

- Run a single method for one example:

  result <- run_method(1, "IS-VS", n_reps = 25, ncores = 25)

- Run all methods for one example:

  ex1_results <- run_example(1, n_reps = 25, ncores = 25)


Estimated runtime:
Approximately 14 hours using 25 CPU cores and 25 replications.

Available methods:
IS-VS, IS-CE, IS-Pareto, WAMK-SIS, Lasso, SpAM, RF-RFE, OptiTreeStrat, RF-VS, and GP-VS

----------------------------------------------------------------------
B. Reproducing the ECL Baseline 
----------------------------------------------------------------------

The ECL baseline is implemented in Python and provided separately in the
`Python_ECL` folder.

To reproduce the ECL results, open and run the following Jupyter notebooks:

- ecl_example1.ipynb   (Example 1, D = 4)
- ecl_example2.ipynb   (Example 2, D = 5)
- ecl_example3.ipynb   (Example 3, D = 10)

Each notebook is self-contained and can be run independently.

Estimated runtime:
Approximately 2 days per example.

----------------------------------------------------------------------
C. Reproducing Table 3 and 4
----------------------------------------------------------------------

The analysis of IS-VS is provided in the Ex1_analysis.rmd. 

----------------------------------------------------------------------
D. Reproducing the Case Study
----------------------------------------------------------------------

The case study is implemented in Matlab and provided separately in the
`Case_study` folder.

To reproduce the results, open and run the following Jupyter notebooks:

- Experiment_main.m   (WAMK-SIS)
  (Data generated: Ex2_n_exp_25_2023-8-3-23-54.mat)
  
- Experiment_main_CV.m   (IS-VS)
  (Data generated: Ex2_n_exp_25_2024-12-23-13-40.mat)

Each notebook is self-contained and can be run independently.

Estimated runtime:
Approximately 7 days.


----------------------------------------------------------------------
E. Software Requirements
----------------------------------------------------------------------

R version: 4.2.2

Required R packages include:
doParallel, truncnorm, splines, rootSolve, rmutil, pracma, mvtnorm,
cubature, caret, stats, lhs, LHD, R.utils, crch, loo, GauPro, DiceKriging, kernlab, 
SAM, glmnet, randomForest


Python version: 3.11

Required Python packages include:
numpy, pandas, pyDOE, scipy, scikit-learn, matplotlib, dill

----------------------------------------------------------------------
F. Notes
----------------------------------------------------------------------

The implementation of the ECL baseline method is adapted from publicly
available reference implementations associated with prior work, with
modifications to ensure consistency with the experimental setup used
in this paper.

----------------------------------------------------------------------
License
----------------------------------------------------------------------

MIT License
