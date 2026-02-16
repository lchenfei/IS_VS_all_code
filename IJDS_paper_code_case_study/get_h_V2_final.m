%% Correponding R code
% get_h_V2_final = function(X_mat,Z,h1,h2,xmin=-5,xmax=5,ymin=-5,ymax=5){
%   h_mat = c(h1,h2)
%   for (i in 1:3){
%     h = get_h_V2(X_mat,Z,h1,h2,xmin=-5,xmax=5,ymin=-5,ymax=5)
%     h1 = h[1]
%     h2 = h[2]
%     h_mat = rbind(h_mat,h)
%   }
%   return(h_mat)
% }

function [h_mat] = get_h_V2_final(X_mat,Z,h1,h2,xmin,xmax,ymin,ymax)
    h_mat = [h1,h2];
    for i=1:3
        h = get_h_V2(X_mat,Z,h1,h2,xmin,xmax,ymin,ymax);
        h1 = h(1);
        h2 = h(2);
        h_mat = [h_mat;h];
    end
end

%% Corresponding R code
% get_h_V2 = function(X_mat,Z,h1,h2,xmin=-5,xmax=5,ymin=-5,ymax=5){
%   a = left_1_power4(X_mat,Z,h1,h2,xmin,xmax,ymin,ymax)
%   b = left_2_power4(X_mat,Z,h1,h2,xmin,xmax,ymin,ymax)
%   c = left_power2(X_mat,Z,h1,h2,xmin,xmax,ymin,ymax)
%   d = right_coef(X_mat,Z,h1,h2,xmin,xmax,ymin,ymax)
%   
%   h_out = optim(c(h1,h2),function(h) a*h[1]^4+b*h[2]^4+c*h[1]^2*h[2]^2+d/h[1]/h[2])
%   h_out = h_out$par
%   return(h_out)
% }

function [h_out]= get_h_V2(X_mat,Z,h1,h2,xmin,xmax,ymin,ymax)
  a = left_1_power4(X_mat,Z,h1,h2,xmin,xmax,ymin,ymax);
  b = left_2_power4(X_mat,Z,h1,h2,xmin,xmax,ymin,ymax);
  c = left_power2(X_mat,Z,h1,h2,xmin,xmax,ymin,ymax);
  d = right_coef(X_mat,Z,h1,h2,xmin,xmax,ymin,ymax);
  
  h_out = fminsearch(@(h) a*h(1)^4+b*h(2)^4+c*h(1)^2*h(2)^2+d/h(1)/h(2),[h1,h2]);
  interval = zeros(1,2);
  interval(1) = xmax - xmin;
  interval(2) = ymax - ymin;
  for i = 1:2
      if h_out(i) >= interval(i)
          h_out(i) = interval(i);
      end
  end    
end

%% Corresponding R code
% left_power2 = function(X_mat,Z,h1,h2,xmin=-5,xmax=5,ymin=-5,ymax=5){
%   l = integral2(function(x,y) left_power2_integrand(x,y,X_mat,Z,h1,h2),xmin,xmax,ymin,ymax)
%   return(l$Q)
% }

function [coef] = left_power2(X_mat,Z,h1,h2,xmin,xmax,ymin,ymax)
    coef = integral2(@(x,y) left_power2_integrand(x,y,X_mat,Z,h1,h2),xmin,xmax,ymin,ymax,'RelTol',1e-3);
end

%% Corresponding R code
% # left integral for h1^2*h2^2
% left_power2_integrand = function(x1,x2,X_mat,Z,h1,h2){
%   n = nrow(X_mat)
%   x_mat = matrix(c(x1,x2),ncol = 2)
%   coef = (coef_hsq(x_mat,X_mat,Z,h1,h2))
%   integrand = 2*coef[,1]*coef[,2]
%   return(integrand)
% }

function [integrand] = left_power2_integrand(x1,x2,X_mat,Z,h1,h2)
    shape = size(x1);
    N = size(X_mat,1);
    x1 = reshape(x1,[],1);
    x2 = reshape(x2,[],1);
    x_mat = [x1,x2];
    coef_all = coef_hsq(x_mat,X_mat,Z,h1,h2);
    integrand = 2*coef_all(:,1).*coef_all(:,2);
    integrand = reshape(integrand,shape);
end


%% Corresponding R code
% left_2_power4 = function(X_mat,Z,h1,h2,xmin=-5,xmax=5,ymin=-5,ymax=5){
%   l2 = integral2(function(x,y) left_power4_integrand(x,y,X_mat,Z,h1,h2,2),xmin,xmax,ymin,ymax)
%   return(l2$Q)
% }

function [coef] = left_2_power4(X_mat,Z,h1,h2,xmin,xmax,ymin,ymax)
    coef = integral2(@(x,y) left_power4_integrand(x,y,X_mat,Z,h1,h2,2),xmin,xmax,ymin,ymax,'RelTol',1e-3);
end

%% Corresponding R code
% left_1_power4 = function(X_mat,Z,h1,h2,xmin=-5,xmax=5,ymin=-5,ymax=5){
%   l1 = integral2(function(x,y) left_power4_integrand(x,y,X_mat,Z,h1,h2,1),xmin,xmax,ymin,ymax)
%   return(l1$Q)
% }

function [coef] = left_1_power4(X_mat,Z,h1,h2,xmin,xmax,ymin,ymax)
    coef = integral2(@(x,y) left_power4_integrand(x,y,X_mat,Z,h1,h2,1),xmin,xmax,ymin,ymax,'RelTol',1e-3);
end

%% Corresponding R code
% # left integral for h1^4 and h2^4
% left_power4_integrand = function(x1,x2,X_mat,Z,h1,h2,ind){
%   n = nrow(X_mat)
%   x_mat = matrix(c(x1,x2),ncol = 2)
%   coef = (coef_hsq(x_mat,X_mat,Z,h1,h2))[,ind]
%   integrand = coef^2
%   return(integrand)
% }

function [integrand] = left_power4_integrand(x1,x2,X_mat,Z,h1,h2,ind)
    shape = size(x1);
    N = size(X_mat,1);
    x1 = reshape(x1,[],1);
    x2 = reshape(x2,[],1);
    x_mat = [x1,x2];
    coef_all = coef_hsq(x_mat,X_mat,Z,h1,h2);
    coef = coef_all(:,ind);
    integrand = coef.^2;
    integrand = reshape(integrand,shape);
end









%% Corresponding R code
% right = function(X_mat,Z,h1,h2,xmin=-5,xmax=5,ymin=-5,ymax=5){
%   r = integral2(function(x,y) right_integrand(x,y,X_mat,Z,h1,h2),xmin,xmax,ymin,ymax)
%   return(r$Q)
% }
% # right coef
% right_coef = function(X_mat,Z,h1,h2,xmin=-5,xmax=5,ymin=-5,ymax=5){
%   n = nrow(X_mat)
%   r = right(X_mat,Z,h1,h2,xmin=-5,xmax=5,ymin=-5,ymax=5)
%   return(r/(4*pi)/n)
% }

function[coef] = right_coef(X_mat,Z,h1,h2,xmin,xmax,ymin,ymax)
    r = integral2(@(x,y) right_integrand(x,y,X_mat,Z,h1,h2),xmin,xmax,ymin,ymax,'RelTol',1e-3);
    N = size(X_mat,1);
    coef = r/(4*pi)/N;
end    


%% Corresponding R code
% # right integral
% right_integrand = function(x1,x2,X_mat,Z,h1,h2){
%   n = nrow(X_mat)
%   x_mat = matrix(c(x1,x2),ncol = 2)
%   S = S_bi(x_mat,X_mat,Z,h1,h2)
%   f = 1/(n*2*pi*h1*h2)*S2_bi(x_mat,X_mat,h1,h2)
%   int = S*(1-S)/f
%   int[is.na(int)]=0
%   int[int==Inf]=0
%   return(int)
% }

function [int] = right_integrand(x1,x2,X_mat,Z,h1,h2)
    shape = size(x1);
    N = size(X_mat,1);
    x1 = reshape(x1,[],1);
    x2 = reshape(x2,[],1);
    x_mat = [x1,x2];
    S = S_bi(x_mat,X_mat,Z,h1,h2);
    f = 1/(N*2*pi*h1*h2)*S2_bi(x_mat,X_mat,h1,h2);
    int = S.*(1-S)./f;
    int(isinf(int)|isnan(int)) = 0;
    int = reshape(int,shape);
end

%% Corresponding R code
% # calculate the coefficient of h1^2 and h2^2
% coef_hsq = function(x_mat,X_mat,Z,h1,h2){
%   n = nrow(X_mat)
%   coef = S_d1_bi(x_mat,X_mat,Z,h1,h2)*S2_d1_bi(x_mat,X_mat,h1,h2)/S2_bi(x_mat,X_mat,h1,h2) +
%     S_d2_bi(x_mat,X_mat,Z,h1,h2)/2
%   return(coef)
% }
function [coef] = coef_hsq(x_mat,X_mat,Z,h1,h2)
    S2_hat = S2_bi(x_mat,X_mat,h1,h2);
    S2_hat = repmat(S2_hat,[1,2]);
    coef = S_d1_bi(x_mat,X_mat,Z,h1,h2).*S2_d1_bi(x_mat,X_mat,h1,h2)./S2_hat+...
     S_d2_bi(x_mat,X_mat,Z,h1,h2)/2;
    coef(isinf(coef)|isnan(coef)) = 0;
end

%% Corresponding R code
% # calculate the second derivative of S
% S_d2_bi = function(x_mat,X_mat,Z,h1,h2){
%   S1_hat = S1_bi(x_mat,X_mat,Z,h1,h2)
%   S2_hat = S2_bi(x_mat,X_mat,h1,h2)
%   S1_hat_d1 = S1_d1_bi(x_mat,X_mat,Z,h1,h2)
%   S2_hat_d1 = S2_d1_bi(x_mat,X_mat,h1,h2)
%   S1_hat_d2 = S1_d2_bi(x_mat,X_mat,Z,h1,h2)
%   S2_hat_d2 = S2_d2_bi(x_mat,X_mat,h1,h2)
%   
%   
%   n = (S1_hat_d2*S2_hat-S2_hat_d2*S1_hat)*S2_hat^2-
%     2*S2_hat*S2_hat_d1*(S1_hat_d1*S2_hat-S2_hat_d1*S1_hat)
%   
%   S_hat_d2 = n/(S2_hat^4)
%   S_hat_d2[is.na(S_hat_d2)] = 0
%   S_hat_d2[S_hat_d2==Inf] = 0
%   
%   return(S_hat_d2)
% }

function [S_d2_vec] = S_d2_bi(x_mat,X_mat,Z,h1,h2)
  S1_hat = S1_bi(x_mat,X_mat,Z,h1,h2);
  S2_hat = S2_bi(x_mat,X_mat,h1,h2);
  S1_hat_d1 = S1_d1_bi(x_mat,X_mat,Z,h1,h2);
  S2_hat_d1 = S2_d1_bi(x_mat,X_mat,h1,h2);
  S1_hat_d2 = S1_d2_bi(x_mat,X_mat,Z,h1,h2);
  S2_hat_d2 = S2_d2_bi(x_mat,X_mat,h1,h2);
  
  S1_hat = repmat(S1_hat,[1,2]);
  S2_hat = repmat(S2_hat,[1,2]);
  
  n = (S1_hat_d2.*S2_hat-S2_hat_d2.*S1_hat).*S2_hat.^2-...
    2*S2_hat.*S2_hat_d1.*(S1_hat_d1.*S2_hat-S2_hat_d1.*S1_hat);
  
  S_d2_vec = n./(S2_hat.^4);
  S_d2_vec(isinf(S_d2_vec)|isnan(S_d2_vec)) = 0;
  
end














%% Correponding R code
% # calculate the second derivative of S2
% S2_d2_bi = function(x_mat,X_mat,h1,h2){
%   grad_mat = matrix(NA, ncol = 2, nrow = nrow(x_mat))
%   for (i in 1:nrow(x_mat)){
%     x = x_mat[i,]
%     X1 = X_mat[,1]
%     X2 = X_mat[,2]
%     grad1 = sum(exp(-(X1-x[1])^2/(2*h1^2)-
%                       (X2-x[2])^2/(2*h2^2))*((X1-x[1])^2/h1^4-1/h1^2))
%     grad2 = sum(exp(-(X1-x[1])^2/(2*h1^2)-
%                       (X2-x[2])^2/(2*h2^2))*((X2-x[2])^2/h2^4-1/h2^2))
%     grad_mat[i,] = c(grad1,grad2)
%   }
%   return(grad_mat)
% }

function [grad_mat] = S2_d2_bi(x_mat,X_mat,h1,h2)
    N = size(x_mat,1);
    grad_mat = zeros(N,2);
    for i = 1:N
        x1 = x_mat(i,1);
        x2 = x_mat(i,2);
        X1 = X_mat(:,1);
        X2 = X_mat(:,2);
        grad1 = sum(exp(-(X1-x1).^2/(2*h1^2)-...
                       (X2-x2).^2/(2*h2^2)).*((X1-x1).^2/h1^4-1/h1^2));
        
        grad2 = sum(exp(-(X1-x1).^2/(2*h1^2)-...
                       (X2-x2).^2/(2*h2^2)).*((X2-x2).^2/h2^4-1/h2^2));
        grad_mat(i,:) = [grad1,grad2];
    end
end


%% Corresponding R code
% # calculate the second derivative of S1
% S1_d2_bi = function(x_mat,X_mat,Z,h1,h2){
%   grad_mat = matrix(NA, ncol = 2, nrow = nrow(x_mat))
%   for (i in 1:nrow(x_mat)){
%     x = x_mat[i,]
%     X1 = X_mat[,1]
%     X2 = X_mat[,2]
%     grad1 = sum(exp(-(X1-x[1])^2/(2*h1^2)-
%                       (X2-x[2])^2/(2*h2^2))*((X1-x[1])^2/h1^4-1/h1^2)*Z)
%     grad2 = sum(exp(-(X1-x[1])^2/(2*h1^2)-
%                       (X2-x[2])^2/(2*h2^2))*((X2-x[2])^2/h2^4-1/h2^2)*Z)
%     grad_mat[i,] = c(grad1,grad2)
%   }
%   return(grad_mat)
% }

function [grad_mat] = S1_d2_bi(x_mat,X_mat,Z,h1,h2)
    N = size(x_mat,1);
    grad_mat = zeros(N,2);
    for i = 1:N
        x1 = x_mat(i,1);
        x2 = x_mat(i,2);
        X1 = X_mat(:,1);
        X2 = X_mat(:,2);
        grad1 = sum(exp(-(X1-x1).^2/(2*h1^2)-...
                       (X2-x2).^2/(2*h2^2)).*((X1-x1).^2/h1^4-1/h1^2).*Z);
        
        grad2 = sum(exp(-(X1-x1).^2/(2*h1^2)-...
                       (X2-x2).^2/(2*h2^2)).*((X2-x2).^2/h2^4-1/h2^2).*Z);
        grad_mat(i,:) = [grad1,grad2];
    end
end


%% Corresponding R code
% # calculate the first derivative of S
% S_d1_bi = function(x_mat,X_mat,Z,h1,h2){
%   S_d1 = (S1_d1_bi(x_mat,X_mat,Z,h1,h2)*S2_bi(x_mat,X_mat,h1,h2)-
%           S2_d1_bi(x_mat,X_mat,h1,h2)*S1_bi(x_mat,X_mat,Z,h1,h2))/
%         (S2_bi(x_mat,X_mat,h1,h2))^2
%   S_d1[is.na(S_d1)] = 0
%   S_d1[S_d1==Inf] = 0
%   return(S_d1)
% }
function [S_d1_vec] = S_d1_bi(x_mat,X_mat,Z,h1,h2)
    S1_hat = S1_bi(x_mat,X_mat,Z,h1,h2);
    S2_hat = S2_bi(x_mat,X_mat,h1,h2);
    S1_hat_d1 = S1_d1_bi(x_mat,X_mat,Z,h1,h2);
    S2_hat_d1 = S2_d1_bi(x_mat,X_mat,h1,h2);
    
    S1_hat = repmat(S1_hat,[1,2]);
    S2_hat = repmat(S2_hat,[1,2]);
    
    S_d1_vec = (S1_hat_d1.*S2_hat-S2_hat_d1.*S1_hat)./S2_hat.^2;
    S_d1_vec(isinf(S_d1_vec)|isnan(S_d1_vec)) = 0;
end


%% Corresponding R code
% # Calculate the the first derivative of S2
% S2_d1_bi = function(x_mat,X_mat,h1,h2){
%   grad_mat = matrix(NA, ncol = 2, nrow = nrow(x_mat))
%   for (i in 1:nrow(x_mat)){
%     x = x_mat[i,]
%     X1 = X_mat[,1]
%     X2 = X_mat[,2]
%     grad1 = sum(exp(-(X1-x[1])^2/(2*h1^2)-
%                       (X2-x[2])^2/(2*h2^2))*(X1-x[1])/h1^2)
%     grad2 = sum(exp(-(X1-x[1])^2/(2*h1^2)-
%                       (X2-x[2])^2/(2*h2^2))*(X2-x[2])/h2^2)
%     grad_mat[i,] = c(grad1,grad2)
%   }
%   return(grad_mat)
% }
function [grad_mat] = S2_d1_bi(x_mat,X_mat,h1,h2)
    N = size(x_mat,1);
    grad_mat = zeros(N,2);
    for i=1:N
        x1 = x_mat(i,1);
        x2 = x_mat(i,2);
        X1 = X_mat(:,1);
        X2 = X_mat(:,2);
        grad1 = sum(exp(-(X1-x1).^2/(2*h1^2)- ...
                    (X2-x2).^2/(2*h2^2)).*(X1-x1)/h1^2);
        grad2 = sum(exp(-(X1-x1).^2/(2*h1^2)- ...
                    (X2-x2).^2/(2*h2^2)).*(X2-x2)/h2^2);
        grad_mat(i,:) = [grad1,grad2];
    end
end

%% Corresponding R code
% # calculate the first derivative of S1(x)
% S1_d1_bi = function(x_mat,X_mat,Z,h1,h2){
%   grad_mat = matrix(NA, ncol = 2, nrow = nrow(x_mat))
%   for (i in 1:nrow(x_mat)){
%     x = x_mat[i,]
%     X1 = X_mat[,1]
%     X2 = X_mat[,2]
%     grad1 = sum(exp(-(X1-x[1])^2/(2*h1^2)-
%                      (X2-x[2])^2/(2*h2^2))*Z*(X1-x[1])/h1^2)
%     grad2 = sum(exp(-(X1-x[1])^2/(2*h1^2)-
%                       (X2-x[2])^2/(2*h2^2))*Z*(X2-x[2])/h2^2)
%     grad_mat[i,] = c(grad1,grad2)
%   }
%   return(grad_mat)
% }

function [grad_mat] = S1_d1_bi(x_mat,X_mat,Z,h1,h2)
    N = size(x_mat,1);
    grad_mat = zeros(N,2);
    for i=1:N
        x1 = x_mat(i,1);
        x2 = x_mat(i,2);
        X1 = X_mat(:,1);
        X2 = X_mat(:,2);
        grad1 = sum(exp(-(X1-x1).^2/(2*h1^2)- ...
                    (X2-x2).^2/(2*h2^2)).*Z.*(X1-x1)/h1^2);
        grad2 = sum(exp(-(X1-x1).^2/(2*h1^2)- ...
                    (X2-x2).^2/(2*h2^2)).*Z.*(X2-x2)/h2^2);
        grad_mat(i,:) = [grad1,grad2];
    end
end

%% Corresponding R code
% # The function of S
% S_bi = function(x_mat,X_mat,Z,h1,h2){
%   S1_hat = S1_bi(x_mat,X_mat,Z,h1,h2)
%   S2_hat = S2_bi(x_mat,X_mat,h1,h2)
%   S_hat = S1_hat/S2_hat
%   S_hat[is.na(S_hat)] = 0
%   S_hat[S_hat==Inf] = 0
%   return(S_hat)
% }
function [S_vec] = S_bi(x_mat,X_mat,Z,h1,h2)
    % X is X, Z is (Y>l), h is bandwidth
    S1_vec = S1_bi(x_mat,X_mat,Z,h1,h2);
    S2_vec = S2_bi(x_mat,X_mat,h1,h2);
    S_vec = S1_vec./S2_vec;
    S_vec(isinf(S_vec)|isnan(S_vec)) = 0;
end

%% Correponding R code
% # The numerator of S(x)
% S1_bi = function(x_mat,X_mat,Z,h1,h2){
%   S1_b_vec = c()
%   for (i in 1:nrow(x_mat)){
%     x = x_mat[i,]
%     X1 = X_mat[,1]
%     X2 = X_mat[,2]
%     S1_b = sum(exp(-(X1-x[1])^2/(2*h1^2)-
%                  (X2-x[2])^2/(2*h2^2))*Z)
%     S1_b_vec = c(S1_b_vec,S1_b)
%   }
%   return(S1_b_vec)
% }

function [S1_vec] = S1_bi(x_mat,X_mat,Z,h1,h2)
    N = size(x_mat,1); % number of rows in matrix
    S1_vec = [];
    for i=1:N
        x1 = x_mat(i,1);
        x2 = x_mat(i,2);
        X1 = X_mat(:,1);
        X2 = X_mat(:,2);
        S1 = sum(exp(-(X1-x1).^2/(2*h1^2)-(X2-x2).^2/(2*h2^2)).*Z);
        S1_vec = [S1_vec;S1];
    end
end
%% Correponding R code
% # The denominator of S(x)
% S2_bi = function(x_mat,X_mat,h1,h2){
%   S2_b_vec = c()
%   for (i in 1:nrow(x_mat)){
%     x = x_mat[i,]
%     X1 = X_mat[,1]
%     X2 = X_mat[,2]
%     S2_b = sum(exp(-(X1-x[1])^2/(2*h1^2)-
%                      (X2-x[2])^2/(2*h2^2)))
%     S2_b_vec = c(S2_b_vec,S2_b)
%   }
%   return(S2_b_vec)
% }

function [S2_vec] = S2_bi(x_mat,X_mat,h1,h2)
    N = size(x_mat,1); % number of rows in matrix
    S2_vec = [];
    for i=1:N
        x1 = x_mat(i,1);
        x2 = x_mat(i,2);
        X1 = X_mat(:,1);
        X2 = X_mat(:,2);
        S2 = sum(exp(-(X1-x1).^2/(2*h1^2)-(X2-x2).^2/(2*h2^2)));
        S2_vec = [S2_vec;S2];
    end
end








