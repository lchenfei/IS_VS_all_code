n_per_step = 8;
n_per_initstep = 8;
t = 1;
exp_ind = 1;
LoadType = 3;
N_load_type = 3;
n_exp = 25;
ifsave = 1;

n_parallel = 25; %25; 
poolobj = parpool('local',n_parallel);
spmd
s = RandStream.create('mrg32k3a','NumStreams',numlabs,'StreamIndices',labindex, 'Seed', 'shuffle');
RandStream.setGlobalStream(s);
end

X_all = cell(1, 25);
Y_all = cell(1, 25);
tstart = tic;
parfor i = 1:n_exp
    X = sample_from_original(8);
    X_all{i} = X;
    Y = zeros([n_per_step,1]);
    k_type_of_load_for_a_run_i_j = zeros([n_per_initstep,N_load_type]);
    exp_ind = i;
    for index_i = 1:8  %drange(1:N_long)   % index i for each x_i

                NewRandomSeed = round(-2147483648 + (2147483647-(-2147483648))*rand(1)); 
                run_ind = index_i*10^3+t*10^7+exp_ind;
                display(run_ind);
                k_type_of_load_for_a_run_i_j(index_i,:) = NREL_simulator_local_tmp_if_Unix(run_ind,X(index_i,1), X(index_i,2),X(index_i,3),X(index_i,4),X(index_i,5),NewRandomSeed, 2, 0);	%return value: [RootMyb1, RootMxb1]				
                Y(index_i)= k_type_of_load_for_a_run_i_j(index_i,LoadType);

    end
    Y_all{i} = Y;
end
telapsed = toc(tstart);
fprintf('Finish nonparamSIS2 method\n');

fprintf('WAMK-SIS results: ')
fprintf('Elapsed Time = %g seconds = %g minutes = %g hours\n', telapsed, telapsed/60, telapsed/(60*60))
%fprintf('Sample Mean  = %g\n', mean(final_est_list))
%fprintf('Standard Err = %g\n', std(final_est_list))


%% Wrapup
delete(poolobj)
if (ifsave == 1)
    datetime = fix(clock);
    datetime_str = horzcat(num2str(datetime(1)),'-',num2str(datetime(2)),'-',num2str(datetime(3)),'-',num2str(datetime(4)),'-',num2str(datetime(5)));    
    saveFileName = horzcat('Ex2_n_exp_',num2str(n_exp),'_',datetime_str,'.mat');
    save(saveFileName);
end