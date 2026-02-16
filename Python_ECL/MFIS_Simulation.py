"""
Calculates MFIS estimates for failure probability using ECL
adaptive designs from the Hartmann6 experiment and space-filling designs.

@author: D. Austin Cole  austin.cole8@vt.edu
"""

from datetime import date
import numpy as np
import pandas as pd
from pyDOE import lhs
import scipy.stats as ss
from sklearn.gaussian_process import GaussianProcessRegressor as GPR
from sklearn.gaussian_process.kernels import RBF, WhiteKernel
import os
import sys
fileDir = os.getcwd()
sourcePath = os.path.join(fileDir, 'MFISPy')
sys.path.append(sourcePath)
from mfis import BiasingDistribution
from mfis import MultiFidelityIS
#from mfis.mv_independent_distribution import MVIndependentDistribution
from mfis.mv_dependent_distribution import MVDependentDistribution
import multiprocessing

sourcePath = fileDir
sys.path.append(sourcePath)
from eclGP import Numerical_Ex4
from eclGP import *
import matplotlib.pyplot as plt
import time
import dill

case = 'Ex4'
Model = Numerical_Ex4()
bounds = ((-5,5),(-5,5),(-5,5),(-5,5), (-5, 5))
threshold = 24.2
dim = len(bounds)
alpha = 0.01

N_HIGH_FIDELITY_INPUTS = 1200
n_init = 10*dim
n_select = 600 - n_init
MC_REPS = 25
batch_size = 10

mfis_probs = np.ones((MC_REPS, 3))
#gauss_kernel = 1.0 * RBF(length_scale= np.repeat(0.5, dim),length_scale_bounds=(1e-4, 1e4))
# Stochastic
gauss_kernel = 1.0 * RBF(length_scale= np.repeat(0.5, dim),length_scale_bounds=(1e-4, 1e4)) + WhiteKernel(noise_level=1e-2, noise_level_bounds=(1e-10, 1e1))

####################################
# 2. MFIS
####################################
unif_dist = ss.uniform(0,1)
m1 = 0
s1 = 1
norm_dist = ss.norm(loc=m1, scale=s1)

input_distribution = MVDependentDistribution(
    distributions=[sample_ex4])

def limit_state_function(y):
    return y - threshold

## Estimate failure probability (alpha)
'''
num_failures = np.zeros((100,))
for i in range(len(num_failures)):
    new_X = input_distribution.draw_samples(1000000)
    new_Y = hartmann_model.predict(new_X)
    num_failures[i] = np.sum(hartmann_limit_state_function(new_Y)>0)
    print(num_failures[i])
    print(i)
alpha = np.mean(num_failures)/len(new_X)
'''
mfis_probs = np.ones((25, 3))

ecl_designs_df = pd.read_csv(f"ecl_designs_Ex4_070724.csv",header=None)
ecl_designs = np.array(ecl_designs_df)

def ECL_experiment(j):
    start_time = time.time()
    ecl_gp = None
       
    ## ECL
    design_X = ecl_designs[:, ((dim+1)*j):((dim+1)*(j+1)-1)]
    design_Y = ecl_designs[:, (dim+1)*(j+1)-1]
    ecl_gp =  GPR(kernel=gauss_kernel, alpha=1e-6)
    ecl_gp.fit(design_X, design_Y)
    
    ## ECL-Batch
#     design_X_batch = ecl_designs_batch[:, ((dim+1)*j):((dim+1)*(j+1)-1)]
#     design_Y_batch = ecl_designs_batch[:, (dim+1)*(j+1)-1]
#     ecl_gp_batch =  GPR(kernel=gauss_kernel, alpha=1e-6)
#     ecl_gp_batch.fit(design_X_batch, design_Y_batch)
            
  
    # Initialize Biasing Distributions
    ecl_bd =  BiasingDistribution(trained_surrogate=ecl_gp,
                            limit_state=limit_state_function,
                            input_distribution=input_distribution)
    ecl_bd_ucb =  BiasingDistribution(trained_surrogate=ecl_gp,
                            limit_state=limit_state_function,
                            input_distribution=input_distribution)
    
    ## Fit Biasing Distributions
    ecl_failed_inputs = np.empty((0, dim))   
    ecl_failed_inputs_ucb = np.empty((0, dim))

    # Get sample outputs from GPs and classify failures  
    for k in range(100):
        sample_inputs = sample_n(100000)#input_distribution.draw_samples(100000)
        
        ecl_sample_outputs, ecl_sample_std = \
            ecl_gp.predict(sample_inputs, return_std=True)
        ecl_failed_inputs_new = sample_inputs[
            limit_state_function(ecl_sample_outputs.flatten()) > 0,:]
        ecl_failed_inputs = np.vstack((ecl_failed_inputs,
                                       ecl_failed_inputs_new))

        ecl_failed_inputs_ucb_new = sample_inputs[
            limit_state_function(
                ecl_sample_outputs.flatten() + 1.645*ecl_sample_std) > 0,:]
        ecl_failed_inputs_ucb = np.vstack((ecl_failed_inputs_ucb,
                                           ecl_failed_inputs_ucb_new))
        
        if (k % 100) == 0:
            print(k)

    if len(ecl_failed_inputs) < 1:
        mfis_probs[j, 0] = 0
    else:
        ecl_bd.fit_from_failed_inputs(ecl_failed_inputs,
                                   max_clusters=10, covariance_type='diag')
        ## Failure probability estimates
        XX_ecl = ecl_bd.draw_samples(N_HIGH_FIDELITY_INPUTS - n_init - n_select)
        hf_ecl_outputs = Model.predict(XX_ecl)
        multi_IS_ecl = \
            MultiFidelityIS(limit_state=limit_state_function,
                            biasing_distribution=ecl_bd,
                            input_distribution=input_distribution,
                            bounds=bounds)
        ecl_mfis_stats = multi_IS_ecl.get_failure_prob_estimate(XX_ecl,
                                                              hf_ecl_outputs)
        weights_is = multi_IS_ecl.calc_importance_weights(XX_ecl)
        mfis_probs[j, 0] = ecl_mfis_stats[0]
        
    if len(ecl_failed_inputs_ucb) < 1:
        mfis_probs[j, 1] = 0
    else:   
        ecl_bd_ucb.fit_from_failed_inputs(ecl_failed_inputs_ucb,
                           max_clusters=10, covariance_type='diag')
    
        XX_ecl_ucb = \
        ecl_bd_ucb.draw_samples(N_HIGH_FIDELITY_INPUTS - n_init - n_select)
 
        hf_ecl_ucb_outputs = Model.predict(XX_ecl_ucb)
        multi_IS_ecl_ucb = \
            MultiFidelityIS(limit_state=limit_state_function,
                            biasing_distribution=ecl_bd_ucb,
                            input_distribution=input_distribution,
                            bounds=bounds)
        ecl_ucb_mfis_stats = multi_IS_ecl_ucb.get_failure_prob_estimate(
            XX_ecl_ucb, hf_ecl_ucb_outputs)
        weights_ucb = multi_IS_ecl_ucb.calc_importance_weights(XX_ecl_ucb)
        mfis_probs[j, 1] = ecl_ucb_mfis_stats[0]
   


    ## MC

    XX_mc = sample_n(N_HIGH_FIDELITY_INPUTS)#input_distribution.draw_samples(N_HIGH_FIDELITY_INPUTS)

    mc_outputs = Model.predict(XX_mc)
    mc_failures = XX_mc[limit_state_function(mc_outputs) > 0,:]
    mfis_probs[j, 2] = len(mc_failures) / N_HIGH_FIDELITY_INPUTS

    print(j)
    # End timer
    end_time = time.time()
    # Calculate elapsed time
    elapsed_time = end_time - start_time
    print(j,'th Experiment',"Total Elapsed time: ", elapsed_time) 
    
    mfis_probs_df = pd.DataFrame(mfis_probs)
    mfis_probs_df = mfis_probs_df.rename(columns={0:'ECL', 1:'ECL(UCB)',
                                                  2:'MC'})
    return mfis_probs_df
    #mfis_probs_df.to_csv(f'data/{case}_mfis_estimates_ex4_ECL_experiment_try.csv',index=False)
#dill.dump_session(f'{case}_Stochastic_num_{N_HIGH_FIDELITY_INPUTS}_{MC_REPS}rep')

def MFIS(n_parallel):
    pool = multiprocessing.Pool(processes= n_parallel)
    results = pool.map(ECL_experiment, [j for j in range(5)])
    pool.close()
    pool.join()

if __name__ == "__main__":
    sample_estimation = MFIS(5)

    
