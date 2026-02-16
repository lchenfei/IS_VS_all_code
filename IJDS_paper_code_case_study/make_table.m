v_list_all = [];
for i = 1:25
    v_list_all = [v_list_all, cellfun(@(x) mat2str(x), v_list{i}, 'UniformOutput', false)];
end

est_vec_list_all = [];
for i = 1:25
    for k = 1:10
        est_vec_list_all = [est_vec_list_all, est_vec_list{i}(k)];
    end
end

index_list_all = [];
for i = 1:25
    for k = 1:10
        index_list_all = [index_list_all, cellfun(@(x) mat2str(x), v_list{i}(k), 'UniformOutput', false)];
    end
end

%uniqueNumbers = unique(index_list_all);
%meanResponses = accumarray(index_list_all', est_vec_list_all', [], @mean);

% cellArray = index_list_all;
% responses = est_vec_list_all;
% uniqueArrays = cell(0);
% meanResponses = zeros(1, 0);
% responsesAll = zeros(0);
% stdErrors = zeros(1, 0);
% for i = 1:numel(cellArray)
%     currentArray = cellArray{i};
%     currentResponse = responses(i);
%     
%     % Check if the current array is already in the uniqueArrays list
%     idx = find(cellfun(@(x) isequal(x, currentArray), uniqueArrays));
%     if isempty(idx)
%         % If the array is not in the list, add it and initialize the mean response
%         uniqueArrays{end+1} = currentArray;
%         meanResponses(end+1) = currentResponse;
%         responsesAll(end+ 1) = currentResponse;
%         responsesAll
%     else
%         % If the array is already in the list, update the mean response
%         meanResponses(idx) = mean([meanResponses(idx), currentResponse]);
%         responsesAll(idx) = [responsesAll(idx), currentResponse];
%         responsesAll = [responses, currentResponse];
%     end
%     responsesAll
%     %stdErrors(idx) = std(responsesAll);
% end

% Use a containers.Map to group response values by array
responseGroups = containers.Map('KeyType', 'char', 'ValueType', 'any');

arrays = index_list_all;
responses = est_vec_list_all;

for i = 1:length(arrays)
    % Convert the current array to a string as a unique key
    key = mat2str(arrays{i});
    
    % If the key doesn't exist in the map, create a new empty array for it
    if ~isKey(responseGroups, key)
        responseGroups(key) = [];
    end
    
    % Append the corresponding response value to the array in the map
    responseGroups(key) = [responseGroups(key), responses(i)];
end

keys_all = keys(responseGroups);
for i = 1:length(keys_all)
    responseValues = responseGroups(keys_all{i});
    n = length(responseValues);
    stdError = std(responseValues);% / sqrt(n);
    
    % Output the standard error
    fprintf('Standard error for array %s: %f\n', keys_all{i}, stdError);

    mean_responses = mean(responseValues);
    fprintf('mean for array %s: %f\n', keys_all{i}, mean_responses);

    len_responses = length(responseValues);
    fprintf('length for array %s: %f\n', keys_all{i}, len_responses);
end
