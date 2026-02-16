exp_ind = 10;
X = samples{exp_ind}(1:600, :);
Y = responses{exp_ind}(1:600);

%% get initial bandwidth
Z = (Y > target_quantile)*1;

h_new_mat = bandwidth_info{10}{1};

combo = nchoosek(1:5,2);

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

[s_list, v_list] = S_find_v_list_ce_value_order_all(weight, Z, S_all);
    
[select_list, sample_screen] = Pareto_select(X, Z, h_new_mat, weight, v_list, combo, 200);

[v_count, train_list] = S_find_v_list_cv_index_measure_evaluation(X, Z, h_new_mat, select_list, 200);

v_position = find(select_list == v_count); 

sample_pareto = sample_screen((1 + (v_position - 1) * 200): (v_position * 200), :); 
