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
from mfis.mv_independent_distribution import MVIndependentDistribution
import multiprocessing
sourcePath = fileDir
sys.path.append(sourcePath)
from eclGP import Numerical_Ex1
from eclGP import *
import matplotlib.pyplot as plt
import time
import dill

def limit_state_function(y):
    return y - threshold

Model = Numerical_Ex2()
bounds = ((-5,5),(-5,5),(-5,5),(-5,5))
threshold = 18.99
dim = len(bounds)
alpha = 0.01

N_HIGH_FIDELITY_INPUTS = 1200
n_init = 10*dim
n_select = 600 - n_init
MC_REPS = 25
batch_size = 10

unif_dist = ss.uniform(0,1)
m1 = 0
s1 = 1
norm_dist = ss.norm(loc=m1, scale=s1)

input_distribution = MVIndependentDistribution(
    distributions=[norm_dist, norm_dist, norm_dist, norm_dist])

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
mfis_probs = np.ones((MC_REPS, 3))

ecl_designs_df = pd.read_csv(f"data/Ex2/ecl_designs_Ex2_122124.csv",header=None)
ecl_designs = np.array(ecl_designs_df)

# ecl_designs_batch_df = pd.read_csv(f"data/ecl_batch10_designs_{case}_{file_date}.csv",header=None)
# ecl_designs_batch = np.array(ecl_designs_batch_df)
# designs_list = [ecl_designs,ecl_designs_batch]

gauss_kernel = 1.0 * RBF(length_scale= np.repeat(0.5, dim),length_scale_bounds=(1e-4, 1e4)) + WhiteKernel(noise_level=1e-2, noise_level_bounds=(1e-10, 1e1))
  
def find_prob(j):
    # start_time = time.time()
    np.random.seed(j)
    ecl_gp = None
    mfis_probs = np.ones((1, 3))   
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
        sample_inputs = input_distribution.draw_samples(10000)
        
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
        mfis_probs[0, 0] = 0
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
        mfis_probs[0, 0] = ecl_mfis_stats[0]
       
    if len(ecl_failed_inputs_ucb) < 1:
        mfis_probs[0, 1] = 0
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
        mfis_probs[0, 1] = ecl_ucb_mfis_stats[0]
   


    ## MC

    XX_mc = input_distribution.draw_samples(N_HIGH_FIDELITY_INPUTS)

    mc_outputs = Model.predict(XX_mc)
    mc_failures = XX_mc[limit_state_function(mc_outputs) > 0,:]
    mfis_probs[0, 2] = len(mc_failures) / N_HIGH_FIDELITY_INPUTS

    # print(j)
    # End timer
    # end_time = time.time()
    # Calculate elapsed time
    # elapsed_time = end_time - start_time
    # print(j,'th Experiment',"Total Elapsed time: ", elapsed_time) 
    
    mfis_probs_df = pd.DataFrame(mfis_probs)
    mfis_probs_df = mfis_probs_df.rename(columns={0:'ECL', 1:'ECL(UCB)',
                                                  2:'MC'})
    # mfis_probs_df.to_csv(f'data/{case}_mfis_estimates_{file_date}.csv',index=False)
    # dill.dump_session(f'{case}_Stochastic_num_{N_HIGH_FIDELITY_INPUTS}_{MC_REPS}rep')

    return mfis_probs_df


def generate_initial(n_parallel):
    
    pool = multiprocessing.Pool(processes= n_parallel)
    prob = pool.map(find_prob, range(5))
    pool.close()
    pool.join()

    #sample_y = torch.tensor(sample_y, dtype=torch.double)
    #temp_loss = torch.tensor(temp_loss, dtype=torch.double)
    #train_con = torch.tensor(train_con, dtype=torch.double)

    #train_obj = torch.cat((sample_y.reshape(-1,1),temp_loss.reshape(-1,1)),1)
    return prob

if __name__ == "__main__":

    prob = generate_initial(25)

