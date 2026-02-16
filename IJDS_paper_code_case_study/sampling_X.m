function [X_new] = sampling_X(N,X,Z,h_mat,combo,weight, v_list)    
    N_acc = 0;
    X_new = [];
    while(N_acc < N)
        %Step 1: Generate X from q'(x)
        x = sample_from_original(1);
        
        %Step 2: Generate Y from U(0, q'(X)) 
        ori_pdf = original_fx(x);
        y = unifrnd(0,ori_pdf);
        kernel_IS_pdf = ori_pdf.* sqrt(S_add(x,X,Z,h_mat,combo,weight, v_list));
        
        if y<=kernel_IS_pdf
            X_new = [X_new;x];
            N_acc = N_acc+1;
        end
        
end

