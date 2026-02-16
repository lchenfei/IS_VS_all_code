function[index, train_list] = S_find_v_list_cv_index_measure_evaluation(X, Z, h_new_mat, v_list_all, n_per_step)
m = length(v_list_all);
est_cell = zeros(10, m);
est_train = zeros(10, m);
index_all = [];
%n_per_step = length(X)/2;
for k = 1:10
   randomIndices = randperm(length(X));
   train_indices = randomIndices(1:length(X)/2);
   combo = nchoosek(1:5, 2);
   validation_indices = setdiff(randomIndices, train_indices);
   X_train = X(train_indices, :);
   X_validation = X(validation_indices, :);
   Z_train = Z(train_indices);
   Z_test = Z(validation_indices);
    %% Decide weights
    cost_vec = [];
    S_train = zeros(length(Z_train), 10);
    for i = 1:10
        cb = combo(i,:);
        S = S_bi(X_train(:,cb),X_train(:,cb),Z_train,h_new_mat(i,1),h_new_mat(i,2));
        cost = -nansum(Z_train.*log(S)+(1-Z_train).*log(1-S));
        cost_vec = [cost_vec;cost];
        S_train(:, i) = S;
    end
    weight = 1./cost_vec/sum(1./cost_vec);
    % [s_list_train, v_list_train] = S_find_v_list_ce_value_all(weight, Z_train, S_train);
    % index_all = [index_all, v_list_train];

    v_n = 0;

    for v_list_i = v_list_all
         v_n = v_n + 1;
         v_list = index_find(v_list_i);
         weight_v = weight/sum(weight(v_list));
         X_new = sampling_X(n_per_step,X_train,Z_train,h_new_mat,combo,weight_v, v_list);
         S_new = S_add(X_new,X_train,Z_train,h_new_mat,combo,weight_v, v_list);
         S_test = S_add(X_new,X_validation,Z_test,h_new_mat,combo,weight_v, v_list);
         H_new = sqrt(S_new);
         est_sample = sum((1./H_new)/sum(1./H_new).*S_new);
         est_train(k, v_n) = est_sample;
         est_current = sum((1./H_new)/sum(1./H_new).*S_test);
         est_cell(k, v_n) = est_current;
    end

end

std_train = [];

measure = [];
for i = v_list_all
    measure = [measure, length(index_find(i))];
end

train_list = find_train(est_train, measure);

std_est = [];
for i = train_list
    std_est = [std_est, std(est_cell(:, i))];
end
[s_min, index_min] = min(std_est);

index_position = train_list(index_min);
index = v_list_all(index_position);
% est_mat = cell2mat(est_cell);
% est_all_cell{num_iteration} = est_mat;
% std_all = [];
% for i = 1:1023
%     std_all = [std_all, std(est_mat(:, i))];
% end
% [value, index] = min(std_all);
% value_all = [value_all; value];
% index_all = [index_all; index];
end