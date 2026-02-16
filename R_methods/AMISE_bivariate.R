
library(pracma)
library(mvtnorm)
#-------------------------------- Data Generation function
# mu = function(X){
#   m = 20*(1-exp(-0.2*sqrt(1/2*rowSums(X^2))))+
#         (exp(1)-exp(1/2*rowSums(cos(2*pi*X))))
#   return(m)
# }
# 
# get_Y = function(X){
#   n = nrow(X)
#   Y = rnorm(n,mu(X),1)
#   return(Y)
# }
# 
# S_true = function(X,l){
#   S = 1-pnorm(l,mu(X),1)
#   return(S)
# }

# l = 9.37
# X = rmvnorm(5000,mean = c(0,0), diag(4,2))
# Y = get_Y(X)
# Z = (Y>l)*1
# 
# mean(Y>l)
# 
# 
# x = rmvnorm(10, mean = c(0,0), diag(2,2))


# The numerator of S(x)
S1_bi = function(x_mat,X_mat,Z,h1,h2){
  S1_b_vec = c()
  for (i in 1:nrow(x_mat)){
    x = x_mat[i,]
    X1 = X_mat[,1]
    X2 = X_mat[,2]
    S1_b = sum(exp(-(X1-x[1])^2/(2*h1^2)-
                 (X2-x[2])^2/(2*h2^2))*Z)
    S1_b_vec = c(S1_b_vec,S1_b)
  }
  return(S1_b_vec)
}

# The denominator of S(x)
S2_bi = function(x_mat,X_mat,h1,h2){
  S2_b_vec = c()
  for (i in 1:nrow(x_mat)){
    x = x_mat[i,]
    X1 = X_mat[,1]
    X2 = X_mat[,2]
    S2_b = sum(exp(-(X1-x[1])^2/(2*h1^2)-
                     (X2-x[2])^2/(2*h2^2)))
    S2_b_vec = c(S2_b_vec,S2_b)
  }
  return(S2_b_vec)
}

# The function of S
S_bi = function(x_mat,X_mat,Z,h1,h2){
  S1_hat = S1_bi(x_mat,X_mat,Z,h1,h2)
  S2_hat = S2_bi(x_mat,X_mat,h1,h2)
  S_hat = S1_hat/S2_hat
  S_hat[is.na(S_hat)] = 0
  S_hat[S_hat==Inf] = 0
  return(S_hat)
}



# calculate the first derivative of S1(x)
S1_d1_bi = function(x_mat,X_mat,Z,h1,h2){
  grad_mat = matrix(NA, ncol = 2, nrow = nrow(x_mat))
  for (i in 1:nrow(x_mat)){
    x = x_mat[i,]
    X1 = X_mat[,1]
    X2 = X_mat[,2]
    grad1 = sum(exp(-(X1-x[1])^2/(2*h1^2)-
                     (X2-x[2])^2/(2*h2^2))*Z*(X1-x[1])/h1^2)
    grad2 = sum(exp(-(X1-x[1])^2/(2*h1^2)-
                      (X2-x[2])^2/(2*h2^2))*Z*(X2-x[2])/h2^2)
    grad_mat[i,] = c(grad1,grad2)
  }
  return(grad_mat)
}

# Calculate the the first derivative of S2
S2_d1_bi = function(x_mat,X_mat,h1,h2){
  grad_mat = matrix(NA, ncol = 2, nrow = nrow(x_mat))
  for (i in 1:nrow(x_mat)){
    x = x_mat[i,]
    X1 = X_mat[,1]
    X2 = X_mat[,2]
    grad1 = sum(exp(-(X1-x[1])^2/(2*h1^2)-
                      (X2-x[2])^2/(2*h2^2))*(X1-x[1])/h1^2)
    grad2 = sum(exp(-(X1-x[1])^2/(2*h1^2)-
                      (X2-x[2])^2/(2*h2^2))*(X2-x[2])/h2^2)
    grad_mat[i,] = c(grad1,grad2)
  }
  return(grad_mat)
}

# calculate the first derivative of S
S_d1_bi = function(x_mat,X_mat,Z,h1,h2){
  S_d1 = (S1_d1_bi(x_mat,X_mat,Z,h1,h2)*S2_bi(x_mat,X_mat,h1,h2)-
          S2_d1_bi(x_mat,X_mat,h1,h2)*S1_bi(x_mat,X_mat,Z,h1,h2))/
        (S2_bi(x_mat,X_mat,h1,h2))^2
  S_d1[is.na(S_d1)] = 0
  S_d1[S_d1==Inf] = 0
  return(S_d1)
}

# calculate the second derivative of S1
S1_d2_bi = function(x_mat,X_mat,Z,h1,h2){
  grad_mat = matrix(NA, ncol = 2, nrow = nrow(x_mat))
  for (i in 1:nrow(x_mat)){
    x = x_mat[i,]
    X1 = X_mat[,1]
    X2 = X_mat[,2]
    grad1 = sum(exp(-(X1-x[1])^2/(2*h1^2)-
                      (X2-x[2])^2/(2*h2^2))*((X1-x[1])^2/h1^4-1/h1^2)*Z)
    grad2 = sum(exp(-(X1-x[1])^2/(2*h1^2)-
                      (X2-x[2])^2/(2*h2^2))*((X2-x[2])^2/h2^4-1/h2^2)*Z)
    grad_mat[i,] = c(grad1,grad2)
  }
  return(grad_mat)
}

# calculate the second derivative of S2
S2_d2_bi = function(x_mat,X_mat,h1,h2){
  grad_mat = matrix(NA, ncol = 2, nrow = nrow(x_mat))
  for (i in 1:nrow(x_mat)){
    x = x_mat[i,]
    X1 = X_mat[,1]
    X2 = X_mat[,2]
    grad1 = sum(exp(-(X1-x[1])^2/(2*h1^2)-
                      (X2-x[2])^2/(2*h2^2))*((X1-x[1])^2/h1^4-1/h1^2))
    grad2 = sum(exp(-(X1-x[1])^2/(2*h1^2)-
                      (X2-x[2])^2/(2*h2^2))*((X2-x[2])^2/h2^4-1/h2^2))
    grad_mat[i,] = c(grad1,grad2)
  }
  return(grad_mat)
}

# calculate the second derivative of S
S_d2_bi = function(x_mat,X_mat,Z,h1,h2){
  S1_hat = S1_bi(x_mat,X_mat,Z,h1,h2)
  S2_hat = S2_bi(x_mat,X_mat,h1,h2)
  S1_hat_d1 = S1_d1_bi(x_mat,X_mat,Z,h1,h2)
  S2_hat_d1 = S2_d1_bi(x_mat,X_mat,h1,h2)
  S1_hat_d2 = S1_d2_bi(x_mat,X_mat,Z,h1,h2)
  S2_hat_d2 = S2_d2_bi(x_mat,X_mat,h1,h2)
  
  
  n = (S1_hat_d2*S2_hat-S2_hat_d2*S1_hat)*S2_hat^2-
    2*S2_hat*S2_hat_d1*(S1_hat_d1*S2_hat-S2_hat_d1*S1_hat)
  
  S_hat_d2 = n/(S2_hat^4)
  S_hat_d2[is.na(S_hat_d2)] = 0
  S_hat_d2[S_hat_d2==Inf] = 0
  
  return(S_hat_d2)
}

# calculate the coefficient of h1^2 and h2^2
coef_hsq = function(x_mat,X_mat,Z,h1,h2){
  n = nrow(X_mat)
  coef = S_d1_bi(x_mat,X_mat,Z,h1,h2)*S2_d1_bi(x_mat,X_mat,h1,h2)/S2_bi(x_mat,X_mat,h1,h2) +
    S_d2_bi(x_mat,X_mat,Z,h1,h2)/2
  coef[is.na(coef)]=0
  coef[coef==Inf]=0
  return(coef)
}

#----------- calculate the right coef
# right integral
right_integrand = function(x1,x2,X_mat,Z,h1,h2){
  n = nrow(X_mat)
  x_mat = matrix(c(x1,x2),ncol = 2)
  S = S_bi(x_mat,X_mat,Z,h1,h2)
  f = 1/(n*2*pi*h1*h2)*S2_bi(x_mat,X_mat,h1,h2)
  int = S*(1-S)/f
  int[is.na(int)]=0
  int[int==Inf]=0
  return(int)
}
right = function(X_mat,Z,h1,h2,xmin=-5,xmax=5,ymin=-5,ymax=5){
  r = integral2(function(x,y) right_integrand(x,y,X_mat,Z,h1,h2),xmin,xmax,ymin,ymax)
  return(r$Q)
}
# right coef
right_coef = function(X_mat,Z,h1,h2,xmin=-5,xmax=5,ymin=-5,ymax=5){
  n = nrow(X_mat)
  r = right(X_mat,Z,h1,h2,xmin,xmax,ymin,ymax)
  return(r/(4*pi)/n)
}


#----------- calculate the left coef
# left integral for h1^4 and h2^4
left_power4_integrand = function(x1,x2,X_mat,Z,h1,h2,ind){
  n = nrow(X_mat)
  x_mat = matrix(c(x1,x2),ncol = 2)
  coef = (coef_hsq(x_mat,X_mat,Z,h1,h2))[,ind]
  integrand = coef^2
  return(integrand)
}

left_1_power4 = function(X_mat,Z,h1,h2,xmin=-5,xmax=5,ymin=-5,ymax=5){
  l1 = integral2(function(x,y) left_power4_integrand(x,y,X_mat,Z,h1,h2,1),xmin,xmax,ymin,ymax)
  return(l1$Q)
}

left_2_power4 = function(X_mat,Z,h1,h2,xmin=-5,xmax=5,ymin=-5,ymax=5){
  l2 = integral2(function(x,y) left_power4_integrand(x,y,X_mat,Z,h1,h2,2),xmin,xmax,ymin,ymax)
  return(l2$Q)
}

# left integral for h1^2*h2^2
left_power2_integrand = function(x1,x2,X_mat,Z,h1,h2){
  n = nrow(X_mat)
  x_mat = matrix(c(x1,x2),ncol = 2)
  coef = (coef_hsq(x_mat,X_mat,Z,h1,h2))
  integrand = 2*coef[,1]*coef[,2]
  return(integrand)
}

left_power2 = function(X_mat,Z,h1,h2,xmin=-5,xmax=5,ymin=-5,ymax=5){
  l = integral2(function(x,y) left_power2_integrand(x,y,X_mat,Z,h1,h2),xmin,xmax,ymin,ymax)
  return(l$Q)
}

# Version 2
get_h_V2 = function(X_mat,Z,h1,h2,xmin=-5,xmax=5,ymin=-5,ymax=5){
  a = left_1_power4(X_mat,Z,h1,h2,xmin,xmax,ymin,ymax)
  b = left_2_power4(X_mat,Z,h1,h2,xmin,xmax,ymin,ymax)
  c = left_power2(X_mat,Z,h1,h2,xmin,xmax,ymin,ymax)
  d = right_coef(X_mat,Z,h1,h2,xmin,xmax,ymin,ymax)
  
  h_out = optim(c(h1,h2),function(h) a*h[1]^4+b*h[2]^4+c*h[1]^2*h[2]^2+d/h[1]/h[2], lower = c(0.01,0.01))
  h_out = h_out$par
  return(h_out)
}

get_h_V2_final = function(X_mat,Z,h1,h2,xmin=-4,xmax=4,ymin=-4,ymax=4){
  h_mat = c(h1,h2)
  for (i in 1:3){
    #print(paste('Bandwidth selection, Iteration',i))
    h = get_h_V2(X_mat,Z,h1,h2,xmin,xmax,ymin,ymax)
    h1 = h[1]
    h2 = h[2]
    h_mat = rbind(h_mat,h)
  }
  return(h_mat)
}

# # Version 1
# get_h1 = function(X_mat,Z,h1,h2,xmin=-5,xmax=5,ymin=-5,ymax=5){
#   a = left_1_power4(X_mat,Z,h1,h2,xmin,xmax,ymin,ymax)
#   b = left_2_power4(X_mat,Z,h1,h2,xmin,xmax,ymin,ymax)
#   c = left_power2(X_mat,Z,h1,h2,xmin,xmax,ymin,ymax)
#   d = right_coef(X_mat,Z,h1,h2,xmin,xmax,ymin,ymax)
#   h1_out = optimise(function(x) a*x^4+b*h2^4+c*h2^2*x^2+d/h2/x, interval = c(0,1))
#   return(h1_out$minimum)
# }
# 
# 
# get_h2 = function(X_mat,Z,h1,h2,xmin=-5,xmax=5,ymin=-5,ymax=5){
#   a = left_1_power4(X_mat,Z,h1,h2,xmin,xmax,ymin,ymax)
#   b = left_2_power4(X_mat,Z,h1,h2,xmin,xmax,ymin,ymax)
#   c = left_power2(X_mat,Z,h1,h2,xmin,xmax,ymin,ymax)
#   d = right_coef(X_mat,Z,h1,h2,xmin,xmax,ymin,ymax)
#   h2_out = optimise(function(x) a*h1^4+b*x^4+c*h1^2*x^2+d/h1/x, interval = c(0,1))
#   return(h2_out$minimum)
# }
# 
# 
# get_h = function(X_mat,Z,h1,h2,xmin=-5,xmax=5,ymin=-5,ymax=5){
#   h1_vec = c(h1)
#   h2_vec = c(h2)
#   for (i in 1:5){
#     print(paste('Bandwidth selection, Iteration',i))
#     h1 = get_h1(X_mat,Z,h1,h2,xmin=-5,xmax=5,ymin=-5,ymax=5)
#     h2 = get_h2(X_mat,Z,h1,h2,xmin=-5,xmax=5,ymin=-5,ymax=5)
#     h1_vec = c(h1_vec,h1)
#     h2_vec = c(h2_vec,h2)
#   }
#   return(list(h1=h1_vec,h2=h2_vec))
# }


source('AMISE_univariate.R')


# optimize(function(x) AMWISE(X[,2],Z,x,-5,5),interval = c(0.1,1))
# 
# h_vec = seq(0.1,1,0.1)
# am_vec = c()
# for ( h in h_vec){
#   print(h)
#   am = AMWISE(X[,1],Z,0.03,-5,5)
#   am_vec = c(am_vec,am)
# }
# 
# 
# 
# 
# 
# start_time <- Sys.time()
# H = get_h_V2_final(X,Z,0.48,0.475)
# end_time <- Sys.time()
# 
# plot(H[,2], main = 'h2', xlab = 'iteration', ylab = 'h2', ylim = c(0,0.5))
# plot(H[,1], main = 'h1', xlab = 'iteration', ylab = 'h1', ylim = c(0,0.5))
# 
# H = get_h_V2(X,Z,0.1,0.1)
# 
# start_time <- Sys.time()
# H = get_h(X,Z,0.1644670,0.1714809)
# end_time <- Sys.time()
# 
# end_time - start_time
# 
# h1_track = c(0.4770000,0.3604373,0.3110729,0.2761167,0.2482267,0.2260600,
#              0.2089364,0.1947224,0.1830624,0.1725432,0.1644670,0.1562161,
#              0.1485114,0.1416300,0.1347885,0.1292984)
# 
# h2_track = c(0.4680000,0.3633209,0.3127600,0.2783315,0.2507016,0.2289158,
#              0.2117851,0.1988452,0.1874178,0.1778671,0.1714809,0.1663842,
#              0.1612576,0.1570763,0.1519913,0.1484187)
# 
# plot(1:16,h1_track,ylim = c(0,0.5))
# plot(1:16,h2_track, ylim = c(0,0.5))
# 
# 
# 
# 
# right_coef(X,Z,0.477,0.468)
# 
# coef_hsq(x,X,Z,0.477,0.468)
# 
# S_est = S_bi(x,X,Z,0.477,0.468)



