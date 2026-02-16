% 
% 
% function [h] = get_best_bandwidth(X,Z)
%     fun = @(x)AMISE(X,Z,x);
%     h = fminbnd(fun,0.05,0.3);
% end


%% Corresponding R code
% AMWISE = function(X,Z,h,left,right){
%   R_K = 1/(2*sqrt(pi))
%   n = length(X)
%   b = get_b_w(X,Z,h,left,right)
%   a = get_a_w(X,Z,h,left,right)
%   A = h^4*a+R_K*b/(n*h)
%   return(A)
% }

function [A] = AMISE_uni(X,Z,h,left,right)
    R_K = 1/(2*sqrt(pi));
    n = length(X);
    b = get_b(X,Z,h,left,right);
    a = get_a(X,Z,h,left,right);  
    A = h^4*a+R_K*b/(n*h);
end




%% Corresponding R code
% get_b_w = function(X,Z,h,left,right){
%   b = integral(function(x) integrand_b_w(x,X,Z,h), left,right)
%   return(b)
% }

function [b] = get_b(X,Z,h,left,right)
    integrand = @(x)integrand_b(x,X,Z,h);
    b = quadgk(integrand, left,right);
end




%% Corresponding R code
% get_a_w = function(X,Z,h,left,right){
%   a = integral(function(x) integrand_a_w(x,X,Z,h), left,right)
%   return(a)
% }
function [a] = get_a(X,Z,h,left,right)
    integrand = @(x)integrand_a(x,X,Z,h);
    a = quadgk(integrand, left,right);
end



%% Corresponding R code 
% ######## get b
% integrand_b_w = function(x_vec,X,Z,h){
%   n = length(X)
%   int_b = (S(x_vec,X,Z,h))^1.5*(1-S(x_vec,X,Z,h))*dnorm(x_vec)/S2(x_vec,X,Z,h)*sqrt(2*pi)*n*h
%   int_b[is.na(int_b)] = 0
%   int_b[int_b==Inf] = 0
%   return(int_b)
% }
function [int_b] = integrand_b(x_vec,X,Z,h)
    % f0 means the original density function of X
    n = length(X);
%     f = reshape(f0(x_vec),length(x_vec),1);
    int_b = (S(x_vec,X,Z,h)).*(1-S(x_vec,X,Z,h))./S2(x_vec,X,h)*sqrt(2*pi)*n*h;
    int_b(isinf(int_b)|isnan(int_b)) = 0;
    int_b = int_b';
end


%% Corresponding R code
% integrand_a_w = function(x_vec,X,Z,h){
%   int_a = (S_d2(x_vec,X,Z,h)/2+S_d1(x_vec,X,Z,h)*S2_d1(x_vec,X,Z,h)/S2(x_vec,X,Z,h))^2*dnorm(x_vec)*sqrt(S(x_vec,X,Z,h))
%   int_a[is.na(int_a)] = 0
%   int_a[int_a==Inf] = 0
%   return(int_a)
% }

function [int_a] = integrand_a(x_vec,X,Z,h)
    % f0 means the original density function of X
    %f = reshape(f0(x_vec),length(x_vec),1);
    int_a = (S_d2(x_vec,X,Z,h)/2+S_d1(x_vec,X,Z,h).*S2_d1(x_vec,X,h)./S2(x_vec,X,h)).^2;
    int_a(isinf(int_a)|isnan(int_a)) = 0;
    int_a = int_a';
end


%% corresponding R code
% S = function(x_vec,X,Z,h){
%   S1_hat = S1(x_vec,X,Z,h)
%   S2_hat = S2(x_vec,X,Z,h)
%   S_hat = S1_hat/S2_hat
%   S_hat[is.na(S_hat)] = 0
%   S_hat[S_hat==Inf] = 0
%   return(S_hat)
% }
function [S_vec] = S(x_vec,X,Z,h)
    % X is X, Z is (Y>l), h is bandwidth
    S1_vec = S1(x_vec,X,Z,h);
    S2_vec = S2(x_vec,X,h);
    S_vec = S1_vec./S2_vec;
    S_vec(isinf(S_vec)|isnan(S_vec)) = 0;
end
    
%% corresponding R code
% # numerator of estimate of S(x)
% S1 = function(x_vec,X,Z,h){
%   S1_vec = c()
%   # the length of X and Z should be same
%   for (i in 1:length(x_vec)){
%     x = x_vec[i]
%     S1 = sum(exp(-(X-x)^2/(2*h^2))*Z)
%     S1_vec = c(S1_vec,S1)
%   }
%   return(S1_vec)
% }
function [S1_vec] = S1(x_vec, X, Z, h)
    % X is X, Z is (Y>l), h is bandwidth
    N = length(x_vec);
    S1_vec=[];

    for i=1:N
        diff = x_vec(i)-X;
        S1_temp = sum(exp(-diff.^2./(2*h^2)).*Z);
        S1_vec = [S1_vec;S1_temp];
    end
end
%% Corresponding R code
% # The denominator of S(x)
% S2 = function(x_vec,X,Z,h){
%   S2_vec = c()
%   # the length of X and Z should be same
%   for (i in 1:length(x_vec)){
%     x = x_vec[i]
%     S2 = sum(exp(-(X-x)^2/(2*h^2)))
%     S2_vec = c(S2_vec,S2)
%   }
%   return(S2_vec)
% }
function [S2_vec] = S2(x_vec, X, h)
    % X is X, Z is (Y>l), h is bandwidth
    N = length(x_vec);
    S2_vec=[];

    for i=1:N
        diff = x_vec(i)-X;
        S2_temp = sum(exp(-diff.^2./(2*h^2)));
        S2_vec = [S2_vec;S2_temp];
    end
end


%% Corresponding R code
% # calculate the first derivative of S1(x)
% S1_d1 = function(x_vec,X,Z,h){
%   S1_d1_vec = c()
%   # the length of X and Z should be same
%   for (i in 1:length(x_vec)){
%     x = x_vec[i]
%     S1_d1 = sum(exp(-(X-x)^2/(2*h^2))*Z*(X-x)/h^2)
%     S1_d1_vec = c(S1_d1_vec,S1_d1)
%   }
%   return(S1_d1_vec)
% }

function [S1_d1_vec] = S1_d1(x_vec, X, Z, h)
    N = length(x_vec);
    S1_d1_vec=[];

    for i=1:N
        diff = X-x_vec(i);
        S1_d1_temp = sum(exp(-diff.^2./(2*h^2)).*Z.*diff./h^2);
        S1_d1_vec = [S1_d1_vec;S1_d1_temp];
    end
end

%% Corresponding R code
% # calculate the second derivative of S1(x)
% S1_d2 = function(x_vec,X,Z,h){
%   S1_d2_vec = c()
%   # the length of X and Z should be same
%   for (i in 1:length(x_vec)){
%     x = x_vec[i]
%     S1_d2 = sum( (exp(-(X-x)^2/(2*h^2))*(X-x)^2/h^4-exp(-(X-x)^2/(2*h^2))/h^2)*Z )
%     S1_d2_vec = c(S1_d2_vec,S1_d2)
%   }
%   return(S1_d2_vec)
% }

function [S1_d2_vec] = S1_d2(x_vec, X, Z, h)
    N = length(x_vec);
    S1_d2_vec=[];

    for i=1:N
        diff = X-x_vec(i);
        S1_d2_temp = sum(exp(-diff.^2./(2*h^2)).*(diff.^2./h^4-1/h^2).*Z);
        S1_d2_vec = [S1_d2_vec;S1_d2_temp];
    end
end

%% Corresponding R code
% # calculate the first derivative of S2(x)
% S2_d1 = function(x_vec,X,Z,h){
%   S2_d1_vec = c()
%   # the length of X and Z should be same
%   for (i in 1:length(x_vec)){
%     x = x_vec[i]
%     S2_d1 = sum(exp(-(X-x)^2/(2*h^2))*(X-x)/h^2)
%     S2_d1_vec = c(S2_d1_vec,S2_d1)
%   }
%   return(S2_d1_vec)
% }
function [S2_d1_vec] = S2_d1(x_vec, X, h)
    N = length(x_vec);
    S2_d1_vec=[];

    for i=1:N
        diff = X-x_vec(i);
        S2_d1_temp = sum(exp(-diff.^2./(2*h^2)).*diff./h^2);
        S2_d1_vec = [S2_d1_vec;S2_d1_temp];
    end
end

% %% Corresponding R code
% # calculate the second derivative of S2(x)
% S2_d2 = function(x_vec,X,Z,h){
%   S2_d2_vec = c()
%   # the length of X and Z should be same
%   for (i in 1:length(x_vec)){
%     x = x_vec[i]
%     S2_d2 = sum( (exp(-(X-x)^2/(2*h^2))*(X-x)^2/h^4-exp(-(X-x)^2/(2*h^2))/h^2))
%     S2_d2_vec = c(S2_d2_vec,S2_d2)
%   }
%   return(S2_d2_vec)
% }
function [S2_d2_vec] = S2_d2(x_vec, X, h)
    N = length(x_vec);
    S2_d2_vec=[];

    for i=1:N
        diff = X-x_vec(i);
        S2_d2_temp = sum(exp(-diff.^2./(2*h^2)).*(diff.^2./h^4-1/h^2));
        S2_d2_vec = [S2_d2_vec;S2_d2_temp];
    end
end

%% Corresponding R code
% # calculate the first derivative of S(x)
% S_d1 = function(x_vec,X,Z,h){
%   S1_hat = S1(x_vec,X,Z,h)
%   S2_hat = S2(x_vec,X,Z,h)
%   S1_hat_d1 = S1_d1(x_vec,X,Z,h)
%   S2_hat_d1 = S2_d1(x_vec,X,Z,h)
%   n = (S1_hat_d1*S2_hat-S2_hat_d1*S1_hat)
%   S_hat_d1 = n/(S2_hat^2)
%   S_hat_d1[is.na(S_hat_d1)] = 0
%   S_hat_d1[S_hat_d1==Inf] = 0
%   return(S_hat_d1)
% }
function [S_d1_vec] = S_d1(x_vec,X,Z,h)
    S1_hat = S1(x_vec,X,Z,h);
    S2_hat = S2(x_vec,X,h);
    S1_hat_d1 = S1_d1(x_vec,X,Z,h);
    S2_hat_d1 = S2_d1(x_vec,X,h);
    
    S_d1_vec = (S1_hat_d1.*S2_hat-S2_hat_d1.*S1_hat)./S2_hat.^2;
    S_d1_vec(isinf(S_d1_vec)|isnan(S_d1_vec)) = 0;
end

%% Corresponding R code
% # calculate the first derivative of S(x)
% S_d2 = function(x_vec,X,Z,h){
%   S1_hat = S1(x_vec,X,Z,h)
%   S2_hat = S2(x_vec,X,Z,h)
%   S1_hat_d1 = S1_d1(x_vec,X,Z,h)
%   S2_hat_d1 = S2_d1(x_vec,X,Z,h)
%   S1_hat_d2 = S1_d2(x_vec,X,Z,h)
%   S2_hat_d2 = S2_d2(x_vec,X,Z,h)
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
function [S_d2_vec] = S_d2(x_vec,X,Z,h)
    S1_hat = S1(x_vec,X,Z,h);
    S2_hat = S2(x_vec,X,h);
    S1_hat_d1 = S1_d1(x_vec,X,Z,h);
    S2_hat_d1 = S2_d1(x_vec,X,h);
    S1_hat_d2 = S1_d2(x_vec,X,Z,h);
    S2_hat_d2 = S2_d2(x_vec,X,h);
    
    n = (S1_hat_d2.*S2_hat-S2_hat_d2.*S1_hat).*S2_hat.^2-2*S2_hat.*S2_hat_d1.*(S1_hat_d1.*S2_hat-S2_hat_d1.*S1_hat);
    S_d2_vec = n./S2_hat.^4;
    S_d2_vec(isinf(S_d2_vec)|isnan(S_d2_vec)) = 0;

end








