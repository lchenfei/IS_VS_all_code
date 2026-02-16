function [] = SetDefaultParameters()
global p_var; %number of variablesclear
global target_quantile;
global N_exp; % N_exp = Total_number_of_experiments_for_each_estimator
global N_T; % N = Total computational resource (i.e., total number of sampling execution from the conditional distribution, Y | X = x_i.)
global LoadType; %k=1, Blade 1 root flapwise bending moment ; k=2, Blade 1 root edgewise bending moment


% Experiment setup
p_var = 5; %number of variablesclear
target_quantile = 2.06; %   2.06 for TipDeflection
LoadType = 3; % 1 for flapwise, 2 for edgewise, 3 for TipDeflection

N_exp = 25;%25; % N_exp = Total_number_of_experiments_for_each_estimator
N_T = 2000; % N = Total computational resource (i.e., total number of sampling execution from the conditional distribution, Y | X = x_i.)


end

