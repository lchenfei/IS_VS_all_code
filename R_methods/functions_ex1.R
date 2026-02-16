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
mu = function(X){
  X1 = X[,1]
  X2 = X[,2]
  X3 = X[,3]
  X4 = X[,4]
  m = 40*(1-exp(-0.2*sqrt((X1^2+X2^2)/2)))+
    20*(1-exp(-0.2*abs(X1)))+
    5*(1-exp(-0.2*sqrt((X2^2+X3^2+X4^2)/3)))-
    exp(cos(2*pi*X1*X3))-
    exp(cos(2*pi*X2*X3))-
    exp(cos(2*pi*X1*X2))
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
  while ((nrow(x)<num)||is.null(x)){
    #print(i)
    i = i+1
    x_new = rmvnorm(1,rep(0,D))
    u = runif(1,min=0,max=dmvnorm(x_new))
    if (u<=sqrt(S_add2(x_new,X,Z,h_mat,combo,weight, i_list))*dmvnorm(x_new)){
      x = rbind(x,x_new)
    }
  }
  # x_new = rmvnorm(num, rep(0, D))
  # for(i in 1:D){
  #   if(i %in% combo[i_list, ]){
  #     x_new[, i] = x[, i]
  #   }
  # }
  return(x)
}


# sampling_X_new_GP = function(num, X, Y, l, cb){
#   D = ncol(X)
#   x = c()
#   s = c()
#   i=0
#   X_train = as.matrix(X[, cb])
#   # gp_model <- GauPro(X_train, Y)
#   invisible(capture.output({train_model <- km(formula = ~1,     # Constant mean
#                                               design = as.matrix(X_train), # Input variables
#                                               response = Y, # Output variable
#                                               covtype = "gauss", nugget.estim = TRUE) })) # Kernel type
#   
#   while ((nrow(x)<num)||is.null(x)){
#     #print(i)
#     i = i+1
#     x_new = rmvnorm(1,rep(0,D))
#     # print(x_new)
#     u = runif(1,min=0,max=dmvnorm(x_new))
#     
#     
#     
#     Y_predict = predict(train_model, matrix(x_new[, cb], nrow = 1), type = "UK")
#     
#     # print(Y_predict)
#     
#     mu = Y_predict$mean
#     sigma = Y_predict$sd
#     
#     # print(mu)
#     # print(sigma)
#     
#     
#     
#     S_predict = 1 - pnorm((l - mu)/ sigma)
#     
#     # print(S_predict)
#     
#     if (u<=sqrt(S_predict)*dmvnorm(x_new)){
#       x = rbind(x,x_new)
#       s = rbind(s, S_predict)
#     }
#   }
#   # x_new = rmvnorm(num, rep(0, D))
#   # for(i in 1:D){
#   #   if(i %in% combo[i_list, ]){
#   #     x_new[, i] = x[, i]
#   #   }
#   # }
#   return(list(x, s))
# }

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
    x_new = rmvnorm(1,rep(0,D))
    # print(x_new)
    u = runif(1,min=0,max=dmvnorm(x_new))
    
    
    
    Y_predict = predict(train_model, matrix(x_new[, cb], nrow = 1), type = "UK")
    
    # print(Y_predict)
    
    mu = Y_predict$mean
    sigma = Y_predict$sd
    
    # print(mu)
    # print(sigma)
    
    
    
    S_predict = 1 - pnorm((l - mu)/ sigma)
    
    # print(S_predict)
    
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
    x_new = rmvnorm(1,rep(0,D))
    # print(x_new)
    u = runif(1,min=0,max=dmvnorm(x_new))
    
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
    x_new = rmvnorm(1,rep(0,D))
    x_new = as.data.frame(x_new)
    colnames(x_new) = c("V1", "V2", "V3", "V4")

    u = runif(1,min=0,max=dmvnorm(x_new))
    
    # print(matrix(x_new[, cb], nrow = 1))
    # print(x_new[, cb])
    S_predict = predict(rf_model, x_new[, cb, drop = FALSE], type = "prob")[, "TRUE"]
    
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


get_Est = function(Y,H,l){
  prob = sum((Y>l)*1/H/sum(1/H))
  return(prob)
}
