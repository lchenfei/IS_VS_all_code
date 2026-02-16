function [ final_est,bandwidth_cell, weight_cell, weight_cellv, X,Y, v_cell] = WAMK_SIS( target_quantile,N_T, LoadType,exp_ind)
%% Setup 
n_steps = 1;
n_per_initstep = 600; 
n_per_step = round(N_T./n_steps); 
bandwidth_cell = cell(1,n_steps);
weight_cell = cell(1,n_steps);
weight_cellv = cell(1, n_steps);
v_cell = cell(1, n_steps);
%% t = 0
X1 = 3.0 + (25-3).*rand(n_per_initstep,1);
X2 = 0.05+ 0.4.*rand(n_per_initstep,1);
X3 = -0.4 + 1.4.*rand(n_per_initstep,1);
X4 = -10 + 20 * rand(n_per_initstep,1);
X5 = 0.01 + 0.04 * rand(n_per_initstep,1);
X = [X1,X2,X3,X4,X5];
range = [3, 25,
        0.05, 0.45,
        -0.4, 1,
        -10, 10,
        0.01, 0.05];
hw = (1.0/(22.0*0.4*1.4*20*0.04)).*ones(n_per_initstep,1);
Wx = original_fx(X)./hw;

%% sample Y
N_load_type = 3;
k_type_of_load_for_a_run_i_j = zeros([n_per_initstep,N_load_type]);
Y = zeros([n_per_initstep,1]);
for index_i = 1:n_per_initstep  %drange(1:N_long)   % index i for each x_i
           
            NewRandomSeed = round(-2147483648 + (2147483647-(-2147483648))*rand(1)); 
            run_ind = index_i*10^3+exp_ind;
            display(run_ind);
            k_type_of_load_for_a_run_i_j(index_i,:) = NREL_simulator_local_tmp_if_Unix(run_ind,X(index_i,1),X(index_i,2),X(index_i,3),X(index_i,4),X(index_i,5), NewRandomSeed, 2, 0);    %return value: [RootMyb1, RootMxb1]             
            Y(index_i)= k_type_of_load_for_a_run_i_j(index_i,LoadType);
          
end


%% get initial bandwidth
Z = (Y > target_quantile)*1;


% lower_TI = @(x) mu_TI(x)-0.06;
% upper_TI = @(x) mu_TI(x)+0.06;
% get initial bandwidth
h1 = fminbnd(@(x) AMISE_uni(X(:,1),Z,x,3,25),0.1,2); 
h2 = fminbnd(@(x) AMISE_uni(X(:,2),Z,x,0.05,0.45),0.01,0.1); 
h3 = fminbnd(@(x) AMISE_uni(X(:,3),Z,x,-0.4,1),0.01,0.5); 
h4 = fminbnd(@(x) AMISE_uni(X(:,4),Z,x,-10,10),0.1,2); 
h5 = fminbnd(@(x) AMISE_uni(X(:,5),Z,x,0.01,0.05),0.01,0.02); 
h_init = [h1,h2,h3,h4,h5];

combo = nchoosek(1:5,2);
h_init_mat = h_init(combo);

%% t = 1 through n_steps
est_vec = [];
for t = 1:n_steps     
    %% get optimal bivariate bandwidth
    h_new_mat = [];
    for i = 1:10
        cb = combo(i,:);
        h_init = h_init_mat(i,:);
        h = get_h_V2_final(X(:,cb),Z,h_init(1),h_init(2),range(cb(1),1),range(cb(1),2),range(cb(2),1),range(cb(2),2)); 
        h_new = h(end,:);
        h_new_mat = [h_new_mat;h_new];
    end
    bandwidth_cell{t} = h_new_mat;
    %% Decide weights
    cost_vec = [];
    S_all = zeros(length(Z), 10);
    for i = 1:10
        cb = combo(i,:);
        S = S_bi(X(:,cb),X(:,cb),Z,h_new_mat(i,1),h_new_mat(i,2));
        cost = -nansum(Z.*log(S)+(1-Z).*log(1-S));
        cost_vec = [cost_vec;cost];
        S_all(:, i) = S;
    end

    

    weight = 1./cost_vec/sum(1./cost_vec);
    %disp(weight)
    if sum(isnan(weight))>0
        weight = repelem(1/10,10).';
    end
    weight_cell{t} = weight;
    %v_list = S_find_v_list_ce(weight, Z, S_all);
    v_list = 1:10;
    weight_v = weight/sum(weight(v_list));
    weight_cellv{t} = weight_v;
    v_cell{t} = v_list;
    disp(v_list)
    %% Sampling new X;
    X_new = sampling_X(n_per_step,X,Z,h_new_mat,combo,weight_v, v_list);
    S_new = S_add(X_new,X,Z,h_new_mat,combo,weight_v, v_list);
    H_new = sqrt(S_new);
    %% sample new Y
    k_type_of_load_for_a_run_i_j = zeros([n_per_step,N_load_type]);
    Y_new = zeros([n_per_step,1]);
    for index_i = 1:n_per_step  %drange(1:N_long)   % index i for each x_i

                NewRandomSeed = round(-2147483648 + (2147483647-(-2147483648))*rand(1)); 
                run_ind = index_i*10^3+t*10^7+exp_ind;
                display(run_ind);
                k_type_of_load_for_a_run_i_j(index_i,:) = NREL_simulator_local_tmp_if_Unix(run_ind,X_new(index_i,1), X_new(index_i,2),X_new(index_i,3),X_new(index_i,4),X_new(index_i,5),NewRandomSeed, 2, 0);	%return value: [RootMyb1, RootMxb1]				
                Y_new(index_i)= k_type_of_load_for_a_run_i_j(index_i,LoadType);

    end
    Y = [Y;Y_new];
    %%
    Z_new = (Y_new > target_quantile)*1;
    est_current = sum((1./H_new)/sum(1./H_new).*Z_new);
    est_vec = [est_vec;est_current];
    % aggregate data
    X = [X;X_new];
    Z = [Z;Z_new]; 
    h_init_mat = h_new_mat;
end
    final_est = mean(est_vec);
    fprintf('final est = %g \n', final_est); 
 
end

