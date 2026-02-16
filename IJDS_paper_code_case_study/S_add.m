% S_add = function(x_mat, X_mat, Z, h_mat, combo){
%   S_est = 0
%   for (c in 1:nrow(combo)){
%     cb = combo[c,]
%     h = h_mat[c,]
%     S_temp = S_bi(matrix(x_mat[,cb],ncol=2),X_mat[,cb],Z,h[1],h[2])
%     S_est = S_est+S_temp
%   }
%   S = S_est/nrow(combo)
%   S[is.na(S)] = 0
%   S[S==Inf] = 0
%   return(S)
% }

function [S] = S_add(x_mat,X_mat,Z,h_mat,combo,weight, v_list)
    S_est = 0;
    if isscalar(v_list)
       c = v_list;
       w = weight(c);
       cb = combo(c,:);
       h = h_mat(c,:);
       S_temp = S_bi(x_mat(:,cb),X_mat(:,cb),Z,h(1),h(2));
       S_est = S_est+w*S_temp;
    else
       for i = 1:length(v_list)
           c = v_list(i);
           w = weight(c);
           cb = combo(c,:);
           h = h_mat(c,:);
           S_temp = S_bi(x_mat(:,cb),X_mat(:,cb),Z,h(1),h(2));
           S_est = S_est+w*S_temp;
       end
    end
    S = S_est;
    S(isinf(S)|isnan(S)) = 0;
end