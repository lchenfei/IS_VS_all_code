numRange = 1:10;
frequencyCounts = zeros(length(numRange), 1);
for i = 1:length(v_list)
    % Iterate through each cell in the 10x1 cell array within the current cell
    for j = 1:length(v_list{i})
        % Count occurrences of each number and update frequency count
        for k = numRange
            frequencyCounts(k) = frequencyCounts(k) + sum(v_list{i}{j} == k);
        end
    end
end

est_vec = 0
for i = 1:length(v_list)
    % Iterate through each cell in the 10x1 cell array within the current cell
    for j = 1:10
        if isequal(v_list{i}{j}, [1, 2, 5]) == 1
            est_vec = [est_vec, est_vec_list{i}(j)];
        end
    end
    
end

disp(length(est_vec))
disp(mean(est_vec))
disp(std(est_vec))