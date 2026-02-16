library(cubature)
# numerator of estimate of S(x)
S1 = function(x_vec,X,Z,h){
  S1_vec = c()
  # the length of X and Z should be same
  for (i in 1:length(x_vec)){
    x = x_vec[i]
    S1 = sum(exp(-(X-x)^2/(2*h^2))*Z)
    S1_vec = c(S1_vec,S1)
  }
  return(S1_vec)
}

# calculate the first derivative of S1(x)
S1_d1 = function(x_vec,X,Z,h){
  S1_d1_vec = c()
  # the length of X and Z should be same
  for (i in 1:length(x_vec)){
    x = x_vec[i]
    S1_d1 = sum(exp(-(X-x)^2/(2*h^2))*Z*(X-x)/h^2)
    S1_d1_vec = c(S1_d1_vec,S1_d1)
  }
  return(S1_d1_vec)
}
# calculate the second derivative of S1(x)
S1_d2 = function(x_vec,X,Z,h){
  S1_d2_vec = c()
  # the length of X and Z should be same
  for (i in 1:length(x_vec)){
    x = x_vec[i]
    S1_d2 = sum( (exp(-(X-x)^2/(2*h^2))*(X-x)^2/h^4-exp(-(X-x)^2/(2*h^2))/h^2)*Z )
    S1_d2_vec = c(S1_d2_vec,S1_d2)
  }
  return(S1_d2_vec)
}

# The denominator of S(x)
S2 = function(x_vec,X,Z,h){
  S2_vec = c()
  # the length of X and Z should be same
  for (i in 1:length(x_vec)){
    x = x_vec[i]
    S2 = sum(exp(-(X-x)^2/(2*h^2)))
    S2_vec = c(S2_vec,S2)
  }
  return(S2_vec)
}

# calculate the first derivative of S2(x)
S2_d1 = function(x_vec,X,Z,h){
  S2_d1_vec = c()
  # the length of X and Z should be same
  for (i in 1:length(x_vec)){
    x = x_vec[i]
    S2_d1 = sum(exp(-(X-x)^2/(2*h^2))*(X-x)/h^2)
    S2_d1_vec = c(S2_d1_vec,S2_d1)
  }
  return(S2_d1_vec)
}

# calculate the second derivative of S2(x)
S2_d2 = function(x_vec,X,Z,h){
  S2_d2_vec = c()
  # the length of X and Z should be same
  for (i in 1:length(x_vec)){
    x = x_vec[i]
    S2_d2 = sum( (exp(-(X-x)^2/(2*h^2))*(X-x)^2/h^4-exp(-(X-x)^2/(2*h^2))/h^2))
    S2_d2_vec = c(S2_d2_vec,S2_d2)
  }
  return(S2_d2_vec)
}


S = function(x_vec,X,Z,h){
  S1_hat = S1(x_vec,X,Z,h)
  S2_hat = S2(x_vec,X,Z,h)
  S_hat = S1_hat/S2_hat
  S_hat[is.na(S_hat)] = 0
  S_hat[S_hat==Inf] = 0
  return(S_hat)
}

Practical_S = function(x_vec,X,Z,h){
  S = c()
  for (x in x_vec){
    if (x<(-5)|x>5){
      S_new = 0
    }else{
      diff = (X-x)
      S_new = sum(Z*dnorm(diff,sd=h))/sum(dnorm(diff,sd=h))
    }
    S = c(S, S_new)
  }
  return(S)
}
# calculate the first derivative of S(x)
S_d1 = function(x_vec,X,Z,h){
  S1_hat = S1(x_vec,X,Z,h)
  S2_hat = S2(x_vec,X,Z,h)
  S1_hat_d1 = S1_d1(x_vec,X,Z,h)
  S2_hat_d1 = S2_d1(x_vec,X,Z,h)
  n = (S1_hat_d1*S2_hat-S2_hat_d1*S1_hat)
  S_hat_d1 = n/(S2_hat^2)
  S_hat_d1[is.na(S_hat_d1)] = 0
  S_hat_d1[S_hat_d1==Inf] = 0
  return(S_hat_d1)
}

# calculate the first derivative of S(x)
S_d2 = function(x_vec,X,Z,h){
  S1_hat = S1(x_vec,X,Z,h)
  S2_hat = S2(x_vec,X,Z,h)
  S1_hat_d1 = S1_d1(x_vec,X,Z,h)
  S2_hat_d1 = S2_d1(x_vec,X,Z,h)
  S1_hat_d2 = S1_d2(x_vec,X,Z,h)
  S2_hat_d2 = S2_d2(x_vec,X,Z,h)
  
  
  n = (S1_hat_d2*S2_hat-S2_hat_d2*S1_hat)*S2_hat^2-
    2*S2_hat*S2_hat_d1*(S1_hat_d1*S2_hat-S2_hat_d1*S1_hat)
  
  S_hat_d2 = n/(S2_hat^4)
  S_hat_d2[is.na(S_hat_d2)] = 0
  S_hat_d2[S_hat_d2==Inf] = 0
  
  return(S_hat_d2)
}



integrand_a = function(x_vec,X,Z,h){
  int_a = (S_d2(x_vec,X,Z,h)/2+S_d1(x_vec,X,Z,h)*S2_d1(x_vec,X,Z,h)/S2(x_vec,X,Z,h))^2
  int_a[is.na(int_a)] = 0
  int_a[int_a==Inf] = 0
  return(int_a)
}

get_a = function(X,Z,h,left=-Inf,right=Inf){
  #a = integral(function(x) integrand_a_w(x,X,Z,h), left,right)
  a = adaptIntegrate(function(x) integrand_a(x,X,Z,h), lowerLimit = left, upperLimit = right)
  a = a$integral
  return(a)
}

######## get b
integrand_b = function(x_vec,X,Z,h){
  n = length(X)
  int_b = S(x_vec,X,Z,h)*(1-S(x_vec,X,Z,h))/S2(x_vec,X,Z,h)*sqrt(2*pi)*n*h
  int_b[is.na(int_b)] = 0
  int_b[int_b==Inf] = 0
  return(int_b)
}

get_b = function(X,Z,h,left=-Inf,right=Inf){
  #b = integral(function(x) integrand_b_w(x,X,Z,h), left,right)
  b = adaptIntegrate(function(x) integrand_b(x,X,Z,h), lowerLimit = left, upperLimit = right)
  b = b$integral
  return(b)
}

get_bandwidth = function(h_init,X,Z,left=-Inf, right=Inf, iter=5){
  R_K = 1/(2*sqrt(pi))
  n = length(X)
  h = h_init
  for ( i in 1:iter){
    b = get_b(X,Z,h,left,right)
    a = get_a(X,Z,h,left,right)
    h = (R_K*b/(4*n*a))^{1/5}
    print(h)
  }
  return(h)
}




AMISE = function(X,Z,h,left=-Inf,right=Inf){
  R_K = 1/(2*sqrt(pi))
  n = length(X)
  b = get_b(X,Z,h,left,right)
  a = get_a(X,Z,h,left,right)
  A = h^4*a+R_K*b/(n*h)
  return(A)
}





