kernel_length = [];
for i = 1:25
    for j = 1:10
        kernel_length = [kernel_length, length(v_list{i}{j})];
    end
end