library(mvtnorm)
library(cubature)
# function to get X
get_X = function(num,mu_lst,sigma_lst,pi){
  K = length(pi)
  components = sample(1:K, prob = pi, size = num, replace = T)
  sample_X = c()
  for (n in 1:K){
    num_K = sum(components==n)
    if (num_K==0){
      next
    }
    sample_X = rbind(sample_X,rmvnorm(num_K, mean = mu_lst[[n]], sigma = sigma_lst[[n]]))
  }
  return(sample_X)
}
# fucntion to get Y


#example 10 -extension of example 3 and 7 and 9
mu = function(X){
  X1 = X[,1];  X2 = X[,2];  X3 = X[,3];  X4 = X[,4];  X5 = X[,5] 
  m = 65 - 40 * exp(-0.2 * sqrt((X1^2 + X2^2)/2)) - 
    20 * exp(-0.2 * abs(X1)) - 
    5 * exp(-0.2 * sqrt((X2^2 + X3^2)/2))  - 
    0.1 * exp(-0.2 * sqrt((X4^2 + X5^2)/2)) - 
    exp(cos(2 * pi * X1 * X2)) - exp(cos(2 * pi * X1 * X3)) - exp(cos(2 * pi * X2 * X3))
  return(m)
}




get_Y = function(X){
  n = nrow(X)
  Y = rnorm(n,mu(X),1)
  return(Y)
}



# sampling X using acceptance-rejection algorithm
sampling_X = function(num,X,Z,h_mat,combo,weight){
  D = ncol(X)
  x = c()
  nrow(x)
  i=0
  while ((nrow(x)<num)||is.null(x)){
    #print(i)
    i = i+1
    x_new = rmvnorm(1,rep(0,D))
    u = runif(1,min=0,max=dmvnorm(x_new))
    if (u<=sqrt(S_add(x_new,X,Z,h_mat,combo,weight))*dmvnorm(x_new)){
      x = rbind(x,x_new)
    }
  }
  return(x)
}

sampling_X_new = function(num,X,Z,h_mat,combo,weight, i_list){
  D = ncol(X)
  x = c()
  nrow(x)
  i=0
  #print(combo)
  #set.seed(1)
  while ((nrow(x)<num)||is.null(x)){
    #print(i)
    i = i+1
    #x_new = rmvnorm(1,rep(0,D))
    mu = c(0, 0, 0)
    Sigma = diag(3)
    x_i = rmvnorm(1, mu, Sigma)
    X1 = x_i[, 1]
    X2 = rnorm(1, X1, 1)
    X3 = rnorm(1, X1, 1)
    X4 = x_i[, 2]
    X5 = x_i[, 3]
    x_new = cbind(X1, X2, X3, X4, X5)
    u = runif(1,min=0,max=dmvnorm(x_i) * dmvnorm(X2, X1, diag(1)) * dmvnorm(X3, X1, diag(1)))
    #print(x_new)
    #print(x_new[combo[i_list, ]])
    #print(dmvnorm(x_new))
    #print(dmvnorm(x_new[combo[i_list, ]]) * dmvnorm(x_new[-combo[i_list, ]]))
    #print(sqrt(S_add2(x_new,X,Z,h_mat,combo,weight, i_list)))
    #print(x_new)
    #print(combo)
    if (u<=sqrt(S_add2(x_new,X,Z,h_mat,combo,weight, i_list))*dmvnorm(x_i) * dmvnorm(X2, X1, diag(1)) *  dmvnorm(X3, X1, diag(1))){
      x = rbind(x,x_new)
    }
  }
  #print(x)
  
  # x_new = rmvnorm(num, rep(0, D))
  # for(i in 1:D){
  #   if(i %in% combo[i_list, ]){
  #     x_new[, i] = x[, i] 
  #   }
  # }
  # print(x)
  return(x)
}

sampling_X_new_GP = function(num, X, Y, l, cb){
  D = ncol(X)
  x = c()
  s = c()
  i=0
  X_train = as.matrix(X[, cb])
  # gp_model <- GauPro(X_train, Y)
  invisible(capture.output({train_model <- km(formula = ~1,     # Constant mean
                                              design = as.matrix(X_train), # Input variables
                                              response = Y, # Output variable
                                              covtype = "gauss", nugget.estim = TRUE) })) # Kernel type
  
  while ((nrow(x)<num)||is.null(x)){
    #print(i)
    i = i+1
    mu = c(0, 0, 0)
    Sigma = diag(3)
    x_i = rmvnorm(1, mu, Sigma)
    X1 = x_i[, 1]
    X2 = rnorm(1, X1, 1)
    X3 = rnorm(1, X1, 1)
    X4 = x_i[, 2]
    X5 = x_i[, 3]
    x_new = cbind(X1, X2, X3, X4, X5)
    u = runif(1,min=0,max=dmvnorm(x_i) * dmvnorm(X2, X1, diag(1)) * dmvnorm(X3, X1, diag(1)))
    
    
    
    Y_predict = predict(train_model, matrix(x_new[, cb], nrow = 1), type = "UK")
    
    # print(Y_predict)
    
    mu = Y_predict$mean
    sigma = Y_predict$sd
    
    # print(mu)
    # print(sigma)
    
    
    
    S_predict = 1 - pnorm((l - mu)/ sigma)
    
    # print(S_predict)
    
    if (u<=sqrt(S_predict)*dmvnorm(x_i) * dmvnorm(X2, X1, diag(1)) *  dmvnorm(X3, X1, diag(1))){
      x = rbind(x,x_new)
      s = rbind(s, S_predict)
    }
  }
  # x_new = rmvnorm(num, rep(0, D))
  # for(i in 1:D){
  #   if(i %in% combo[i_list, ]){
  #     x_new[, i] = x[, i]
  #   }
  # }
  return(list(x, s))
}

sampling_X_new_GP_Classification = function(num, X, Z, l, cb){
  D = ncol(X)
  x = c()
  s = c()
  i=0
  X_train = as.matrix(X[, cb])
  
  # print(dim(X_train))
  
  gp_class_model <- gausspr(x = X_train, y = Z, kernel = "rbfdot")
  # S_temp <- predict(gp_class_model, X_train, type = "probabilities")
  
  while ((nrow(x)<num)||is.null(x)){
    #print(i)
    i = i+1
    # x_new = rmvnorm(1,rep(0,D))
    
    mu = c(0, 0, 0)
    Sigma = diag(3)
    x_i = rmvnorm(1, mu, Sigma)
    X1 = x_i[, 1]
    X2 = rnorm(1, X1, 1)
    X3 = rnorm(1, X1, 1)
    X4 = x_i[, 2]
    X5 = x_i[, 3]
    x_new = cbind(X1, X2, X3, X4, X5)
    u = runif(1,min=0,max=dmvnorm(x_i) * dmvnorm(X2, X1, diag(1)) * dmvnorm(X3, X1, diag(1)))
    
    # print(x_new)
    # u = runif(1,min=0,max=dmvnorm(x_new))
    
    # print(matrix(x_new[, cb], nrow = 1))
    
    S_predict = predict(gp_class_model, matrix(x_new[, cb], nrow = 1), type = "probabilities")
    
    # print(Y_predict)
    
    
    # print(S_predict)
    
    if(S_predict < 0){
      next
    }
    
    if (u<=sqrt(S_predict)*dmvnorm(x_new)){
      x = rbind(x,x_new)
      s = rbind(s, S_predict)
    }
  }
  # x_new = rmvnorm(num, rep(0, D))
  # for(i in 1:D){
  #   if(i %in% combo[i_list, ]){
  #     x_new[, i] = x[, i]
  #   }
  # }
  return(list(x, s))
}

sampling_X_new_randomForest = function(num, X, Z, l, cb){
  D = ncol(X)
  x = c()
  s = c()
  i=0
  
  # print(cb)
  
  X_train = as.data.frame(X)
  
  # print(X_train)
  
  z_factor = factor(Z, levels = c(FALSE, TRUE))
  rf_model <- randomForest(z_factor ~ ., data = X_train[, cb, drop = FALSE], 
                           importance = TRUE, ntree = 500)
  
  while ((nrow(x)<num)||is.null(x)){
    #print(i)
    i = i+1
    mu = c(0, 0, 0)
    Sigma = diag(3)
    x_i = rmvnorm(1, mu, Sigma)
    X1 = x_i[, 1]
    X2 = rnorm(1, X1, 1)
    X3 = rnorm(1, X1, 1)
    X4 = x_i[, 2]
    X5 = x_i[, 3]
    x_new = cbind(X1, X2, X3, X4, X5)
    u = runif(1,min=0,max=dmvnorm(x_i) * dmvnorm(X2, X1, diag(1)) * dmvnorm(X3, X1, diag(1)))
    colnames(x_new) = c("V1", "V2", "V3", "V4", "V5")
    
    # u = runif(1,min=0,max=dmvnorm(x_new))
    
    # print(matrix(x_new[, cb], nrow = 1))
    # print(x_new[, cb])
    S_predict = predict(rf_model, x_new[, cb, drop = FALSE], type = "prob")[, "TRUE"]
    
    # print(Y_predict)
    
    
    # print(S_predict)
    
    if(S_predict < 0){
      next
    }
  
    
    if (u<=sqrt(S_predict)*dmvnorm(x_new) * dmvnorm(X2, X1, diag(1)) * dmvnorm(X3, X1, diag(1))){
      x = rbind(x,x_new)
      s = rbind(s, S_predict)
    }
  }
  # x_new = rmvnorm(num, rep(0, D))
  # for(i in 1:D){
  #   if(i %in% combo[i_list, ]){
  #     x_new[, i] = x[, i]
  #   }
  # }
  return(list(x, s))
}

get_Est = function(Y,H,l){
  prob = sum((Y>l)*1/H/sum(1/H))
  return(prob)
}
# # calculate C
# integrand_C = function(x,X_mat,Z,h1,h2){
#   x_mat = matrix(x,ncol = 2)
#   int = sqrt(S_bi(x_mat,X_mat,Z,h1,h2))*dmvnorm(x_mat)
#   return(int)
# }
# 
# get_C = function(X_mat,Z,h1,h2){
#   C = (adaptIntegrate(function(x) integrand_C(x,X_mat,Z,h1,h2),
#                       lowerLimit = c(-Inf,-Inf),
#                       upperLimit = c(Inf,Inf)))$integral
#   return(C)
# }





