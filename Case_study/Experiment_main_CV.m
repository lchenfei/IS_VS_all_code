global savefilename;

global p_var; %number of variablesclear
global target_quantile;

global N_exp; % N_exp = Total_number_of_experiments_for_each_estimator
global N_T; % N = Total computational resource (i.e., total number of sampling execution from the conditional distribution, Y | X = x_i.)

global LoadType;

%% Setup
ifsave = 1; %save the worskpace if 1.
randseed = 1996;
rng(randseed);

n_parallel = 25; %25; 
poolobj = parpool('local',n_parallel);
spmd
s = RandStream.create('mrg32k3a','NumStreams',numlabs,'StreamIndices',labindex, 'Seed', 'shuffle');
RandStream.setGlobalStream(s);
end



%% Compute a quantile based on a large-scale CMC and a refinement by nonlinear solver
SetDefaultParameters();
%%

fprintf('\n Target_quantile is %g \n', target_quantile);

fprintf('WAMK-SIS Method\n');

%%-----------------------------------------
%% WAMK_SIS
%%-----------------------------------------
n_exp = N_exp;

%% CEM
final_est_list = nan(1,n_exp);
samples_CV_update = cell(1,n_exp);
responses_CV_update = cell(1,n_exp);
model_info  = cell(1,n_exp);
bandwidth_info = cell(1,n_exp);
weight_info = cell(1,n_exp);
weightv_info = cell(1, n_exp);
v_list = cell(1, n_exp);
est_vec_list = cell(1, n_exp);
vall_list = cell(1, n_exp);
evaluation_list = cell(1, n_exp); 

target_quantile_= target_quantile;N_T_=N_T;LoadType_=LoadType;
tstart = tic;

% Load samples generated from WAMK-SIS to make sure the same pilot samples.  
parfor i = 1:n_exp
    disp(i)
    [final_est_list(i), est_vec_list{i}, bandwidth_info{i}, weight_info{i}, weightv_info{i}, samples_CV_update{i}, responses_CV_update{i}, v_list{i}, vall_list{i}, evaluation_list{i}] = WAMK_SIS_Pareto_Measure_Evaluation_CV(samples, responses, target_quantile_, N_T_, LoadType_, i);
    fprintf('Finish inside parfor\n');
end
telapsed = toc(tstart);
fprintf('Finish nonparamSIS2 method\n');

fprintf('WAMK-SIS results: ')
fprintf('Elapsed Time = %g seconds = %g minutes = %g hours\n', telapsed, telapsed/60, telapsed/(60*60))
fprintf('Sample Mean  = %g\n', mean(final_est_list))
fprintf('Standard Err = %g\n', std(final_est_list))


%% Wrapup
delete(poolobj)
if (ifsave == 1)
    datetime = fix(clock);
    datetime_str = horzcat(num2str(datetime(1)),'-',num2str(datetime(2)),'-',num2str(datetime(3)),'-',num2str(datetime(4)),'-',num2str(datetime(5)));    
    saveFileName = horzcat('Ex2_n_exp_',num2str(n_exp),'_',datetime_str,'.mat');
    save(saveFileName);
end