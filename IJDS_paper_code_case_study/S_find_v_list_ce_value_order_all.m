function [s_all, index_all] = S_find_v_list_ce_value_order_all(weight, Z, S_all)
    n = length(weight);
    i_all = [];
    ce_all = [];
    for i = 1:n
        index = nchoosek(1:n, i);
        ce_i = [];
        for j = 1:size(index, 1)
            ls = index(j, :);
            ce_i = [ce_i, S_find_ce(weight, Z, S_all, ls)];
            
        end
        ce_all = [ce_all, ce_i];
        [~, i_min] = min(ce_i);
        i_all = [i_all, i_min];
    end
    [s_values, s_index] = sort(ce_all);
    s_all = s_values;
    index_all = s_index;
%     [~, p_i] = min(ce_all);
%     i_index = i_all(p_i);
%     min_index = nchoosek(1:n, p_i);
%     min_index = min_index(i_index, :);
end