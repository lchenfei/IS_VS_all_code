function[train_list] = find_train(est_train, measure)

n = length(measure);
disp(measure);
SE = [];
for i = 1:n
    SE = [SE, std(est_train(:, i))];
end

train_list = [];
    
for i = 1:n
    upper_all = [];
    lower_all = [];
        
     for j = 1:n
            if measure(i) > measure(j) && i ~= j
                upper_all = [upper_all, (SE(i) - SE(j)) / (measure(j) - measure(i))];
            end
            if measure(i) < measure(j) && i ~= j
                lower_all = [lower_all, (SE(i) - SE(j)) / (measure(j) - measure(i))];
            end
            if measure(i) == measure(j) && i ~= j
                upper_all = [upper_all, (SE(j) - SE(i)) / (measure(i) - measure(j))];
            end
     end
        
        % Determine upper and lower bounds
     if isempty(upper_all)
            upper_bound = Inf;
     else
            upper_bound = min(upper_all);
     end
        
     if isempty(lower_all)
            lower_bound = 0;
     else
            lower_bound = max(max(lower_all), 0);
     end
        
        % Display bounds (equivalent to print in R)
        disp(['Lower bound for i = ', num2str(i), ': ', num2str(lower_bound)]);
        disp(['Upper bound for i = ', num2str(i), ': ', num2str(upper_bound)]);
        
        % Check if upper_bound > lower_bound
     if upper_bound > lower_bound
            train_list = [train_list, i];
     end

end
    
    % Display the final result
    disp('Indices satisfying the condition:');
    disp(train_list);

end