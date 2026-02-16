function [select_list, sample_screen] = Pareto_select(X, Z, h_new_mat, weight, v_list, combo, n_per_step)
select_list = [];
sample_screen = [];
for w_count = v_list
    w_list = index_find(w_count);
    weight_w = weight/sum(weight(w_list));
    X_new = sampling_X(n_per_step,X,Z,h_new_mat,combo,weight_w, w_list);
    S_new = S_add(X_new,X,Z,h_new_mat,combo,weight, w_list);
    H_new = sqrt(S_new);
    w = log((1./H_new)/sum(1./H_new));
    [l, k_hat] = psislw(w, 1);
    
    if (k_hat < 0.7)
        select_list = [select_list, w_count];
        sample_screen = [sample_screen;X_new];
    end

    if(length(select_list) >= 8)
        break
    end
    
end

end