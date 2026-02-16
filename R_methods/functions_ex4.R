library(mvtnorm)
library(cubature)
source('Additive_S.R')
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
  X6 = X[,6];  X7 = X[,7];  X8 = X[,8];  X9 = X[,9];  X10 = X[,10]
  m = 40*(1-exp(-0.2*sqrt((X1^2+X2^2)/2)))+
    20*(1-exp(-0.2*abs(X1)))+
    5*(1-exp(-0.2*sqrt((X2^2+X3^2)/2)))- 0.1*exp(-0.2 * sqrt((X4^2+X5^2)/2)) -
    0.1*exp(-0.2 * sqrt((X6^2+X7^2)/2)) - 0.1*exp(-0.2 * sqrt((X8^2+X9^2+X10^2)/3)) - 
    exp(cos(2*pi*X1*X3))-
    exp(cos(2*pi*X2*X3))-
    exp(cos(2*pi*X1*X2))
  for (i  in 11:50){
    m = m - 0.01 * exp(-0.2 * abs(X[, i]))
  }
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
  #set.seed(1)
  while ((nrow(x)<num)||is.null(x)){
    #print(i)
    i = i+1
    x_new = rmvnorm(1,rep(0,D))
    u = runif(1,min=0,max=dmvnorm(x_new))
    #print(x_new)
    #print(x_new[combo[i_list, ]])
    #print(dmvnorm(x_new))
    #print(dmvnorm(x_new[combo[i_list, ]]) * dmvnorm(x_new[-combo[i_list, ]]))
    #print(sqrt(S_add2(x_new,X,Z,h_mat,combo,weight, i_list)))
    # if (length(i_list) >1){
    #  print(S_add2(x_new,X,Z,h_mat,combo,weight, i_list))
    # }
    if (u<=sqrt(S_add2(x_new,X,Z,h_mat,combo,weight, i_list))*dmvnorm(x_new)){
      x = rbind(x,x_new)
      # print(x)
    }
  }
  #print(x)
  
  # x_new = rmvnorm(num, rep(0, D))
  # for(i in 1:D){
  #   if(i %in% combo[i_list, ]){
  #     x_new[, i] = x[, i] 
  #   }
  # }
  #print(x_new)
  return(x)
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

S_find_ce = function(weight, z, S_all, ls){
  weight_ls = weight/sum(weight[ls])
  S_predicted = rep(0, length(z))
  for(i in ls){
    S_predicted = S_predicted + S_all[, i] * weight_ls[i]
  }
  # sum_ce = 0
  # for(i in 1:length(S_predicted)){
  #   if(z[i] == 1){
  #     sum_ce = sum_ce - log2(S_predicted[i])
  #   }
  #   else{
  #     sum_ce = sum_ce - log2(1 - S_predicted[i])
  #   }
  # }
  sum_ce = -sum(z * log(S_predicted) + (1 - z) * log(1 - S_predicted), na.rm = T)
  return(sum_ce)
}
# num, X, Y, l, cb
sampling_ce_order_greedy_GP = function(num, X, z, l, k){
  n = dim(X)[1]
  i_all = c()
  ce_all = c()
  h_s_pareto = c()
  sample_x_pareto_all = c()
  sample_s_pareto_all = c()
  v_list = 1:50
  outer_break = FALSE
  for(j in 1:n){
    combination = combn(1:50, j , simplify = FALSE)
    for(i_row in 1:length(combination)){
        i = combination[[i_row]]
        sample_new = sampling_X_new_GP(num,X,z,l, i)
        # print(i)
        sample_X_new = sample_new[[1]]
        sample_S_new = sample_new[[2]]
        sample_w = (1/sqrt(sample_S_new))/(sum(1/sqrt(sample_S_new)))
        diagnostic_pareto = psis(sample_w, r_eff = 1)
        k_score = diagnostic_pareto$diagnostics$pareto_k
        # print(k_score)
        if(k_score < 0.7){
          h_s_pareto = c(h_s_pareto, list(i))
          sample_x_pareto_all = c(sample_x_pareto_all, list(sample_X_new))
          sample_s_pareto_all = c(sample_s_pareto_all, list(sample_S_new))
        }
        if(length(h_s_pareto) >= k){
          outer_break = TRUE
          break
        }
    }
    if(outer_break){
      break
    }
  }
    
  return(c(list(h_s_pareto), list(sample_x_pareto_all), list(sample_s_pareto_all)))
}


sampling_ce_order_greedy = function(num, X, z, h_mat, combo, weight, S_all, k){
  n = length(weight)
  i_all = c()
  ce_all = c()
  h_s_pareto = c()
  sample_x_pareto_all = c()
  for(i in 1: n){
    ce_i = c()
    if(i == 1){
      index = t(combn(seq(1:n), i))
      ce_i = c()
      for(j in 1:nrow(index)){
        ls = index[j, ]
        ce_i = c(ce_i, S_find_ce(weight, z, S_all, ls))
      }
      ce =  min(ce_i)
      i_a = which.min(ce_i)
      
      in_index = i_a
      ce_all = c(ce_all, ce)
      i_all = c(i_all, i_a)
      weight_i = weight/sum(weight[i_a])
      sample_X_new = sampling_X_new(num,X,z,h_mat,combo,weight_i, i_a)
      
      sample_S_new = S_add2(sample_X_new,X, z, h_mat,combo, weight_i, i_a)
      sample_w = (1/sqrt(sample_S_new))/(sum(1/sqrt(sample_S_new)))
      diagnostic_pareto = psis(sample_w, r_eff = 1)
      k_score = diagnostic_pareto$diagnostics$pareto_k
      print(k_score)
      if(k_score < 0.7){
        h_s_pareto = c(h_s_pareto, list(i_a))
        sample_x_pareto_all = c(sample_x_pareto_all, list(sample_X_new))
      }
      if(length(h_s_pareto) >= k){
        break
      }
    }
    else{
      #index = t(combn(seq(1:n), i))
      in_or_not_index = c()
      in_or_not_j = c()
      # for(j in 1:nrow(index)){
      #   if(in_or_not(i_a, index[j, ])){
      #     in_or_not_index = c(in_or_not_index, i_a)
      #     in_or_not_j = c(in_or_not_j, j)
      #     ls = index[j, ]
      #     ce_i = c(ce_i, S_find_ce(weight, z, S_all, ls))
      #   }
      # }
      for(j in 1:n){
        # print(j)
        if(in_or_not(j, i_a) == F){
          in_or_not_j = c(in_or_not_j, j)
          ls = c(i_a, j)
          #print(ls)
          ce_i = c(ce_i, S_find_ce(weight, z, S_all, ls))
        }
      }

      ce_update =  min(ce_i)
      ce_all = c(ce_all, ce_update)
      in_index_update = in_or_not_j[which.min(ce_i)]
      i_a_update = c(i_a, in_index_update) #index[in_index_update, ]
      i_all = c(i_all, list(i_a_update))
      weight_i = weight/sum(weight[i_a_update])
      sample_X_new = sampling_X_new(num,X,z,h_mat,combo,weight_i, i_a_update)
      # print(sample_X_new)
      
      sample_S_new = S_add2(sample_X_new,X, z, h_mat,combo, weight_i, i_a_update)
      sample_w = (1/sqrt(sample_S_new))/(sum(1/sqrt(sample_S_new)))
      diagnostic_pareto = psis(sample_w, r_eff = 1)
      k_score = diagnostic_pareto$diagnostics$pareto_k
      print(k_score)
      if(k_score < 0.7){
        h_s_pareto = c(h_s_pareto, list(i_a_update))
        sample_x_pareto_all = c(sample_x_pareto_all, list(sample_X_new))
      }
      if(length(h_s_pareto) >= k){
        break
      }
      #print(i_all)
      #if(ce_update >= ce){
      #  min_index = t(combn(seq(1:n), i - 1))[in_index, ]
      #  return(min_index)
      #}
      in_index = in_index_update
      ce = ce_update
      i_a = i_a_update
    }
    
  }
  # print(length(ce_all))
  #min_index = t(combn(seq(1:n), i))[in_index, ]
  # i_list = find_estimator(45, ce_all, i_all)
  return(c(list(h_s_pareto), list(sample_x_pareto_all)))
}

sampling_ce_order_greedy_2 = function(num, X, z, h_mat, combo, weight, S_all, k){
  n = length(weight)
  i_all = c()
  ce_all = c()
  h_s_pareto = c()
  sample_x_pareto_all = c()
  outer_break = FALSE
  for(i in 1: n){
    ce_i = c()
    if(i == 1){
      index = t(combn(seq(1:n), i))
      ce_i = c()
      i_count = c()
      for(j in 1:nrow(index)){
        ls = index[j, ]
        ce_i = c(ce_i, S_find_ce(weight, z, S_all, ls))
        i_count = c(i_count, ls)
      }
      ce =  min(ce_i)
      i_order = order(ce_i)
      k_score = 1
      m = 1
      while(k_score > 0.7){
      in_a = i_order[m] # which.min(ce_i)
      i_a = i_count[in_a]
      ce_all = c(ce_all, ce)
      i_all = c(i_all, i_a)
      weight_i = weight/sum(weight[i_a])
      sample_X_new = sampling_X_new(num,X,z,h_mat,combo,weight_i, i_a)
      
      sample_S_new = S_add2(sample_X_new,X, z, h_mat,combo, weight_i, i_a)
      sample_w = (1/sqrt(sample_S_new))/(sum(1/sqrt(sample_S_new)))
      diagnostic_pareto = psis(sample_w, r_eff = 1)
      k_score = diagnostic_pareto$diagnostics$pareto_k
      print(k_score)
      if(k_score < 0.7){
        h_s_pareto = c(h_s_pareto, list(i_a))
        sample_x_pareto_all = c(sample_x_pareto_all, list(sample_X_new))
      }
      if(length(h_s_pareto) >= k){
        break
      }
      m = m + 1
      }
    }
    else{
      #index = t(combn(seq(1:n), i))
      in_or_not_index = c()
      in_or_not_j = c()
      # for(j in 1:nrow(index)){
      #   if(in_or_not(i_a, index[j, ])){
      #     in_or_not_index = c(in_or_not_index, i_a)
      #     in_or_not_j = c(in_or_not_j, j)
      #     ls = index[j, ]
      #     ce_i = c(ce_i, S_find_ce(weight, z, S_all, ls))
      #   }
      # }
      for(j in 1:n){
        # print(j)
        if(in_or_not(j, i_a) == F){
          in_or_not_j = c(in_or_not_j, j)
          ls = c(i_a, j)
          #print(ls)
          ce_i = c(ce_i, S_find_ce(weight, z, S_all, ls))
        }
      }
      
      ce_update =  min(ce_i)
      ce_all = c(ce_all, ce_update)
      i_order = order(ce_i)
      k_score = 1
      m = 1
      while(k_score > 0.7){
      i_k = i_order[m]
      in_index_update = in_or_not_j[i_k]
      i_a_update = c(i_a, in_index_update) #index[in_index_update, ]
      i_all = c(i_all, list(i_a_update))
      weight_i = weight/sum(weight[i_a_update])
      sample_X_new = sampling_X_new(num,X,z,h_mat,combo,weight_i, i_a_update)
      # print(sample_X_new)
      
      sample_S_new = S_add2(sample_X_new,X, z, h_mat,combo, weight_i, i_a_update)
      sample_w = (1/sqrt(sample_S_new))/(sum(1/sqrt(sample_S_new)))
      diagnostic_pareto = psis(sample_w, r_eff = 1)
      k_score = diagnostic_pareto$diagnostics$pareto_k
      print(k_score)
      if(k_score < 0.7){
        h_s_pareto = c(h_s_pareto, list(i_a_update))
        sample_x_pareto_all = c(sample_x_pareto_all, list(sample_X_new))
      }
      if(length(h_s_pareto) >= k){
        outer_break = TRUE
        break
      }
      m = m + 1
      print(m)
      }
      #print(i_all)
      #if(ce_update >= ce){
      #  min_index = t(combn(seq(1:n), i - 1))[in_index, ]
      #  return(min_index)
      #}
      in_index = in_index_update
      ce = ce_update
      i_a = i_a_update
      print(h_s_pareto)
      if(outer_break){
        break
      }
    }
    
  }
  # print(length(ce_all))
  #min_index = t(combn(seq(1:n), i))[in_index, ]
  # i_list = find_estimator(45, ce_all, i_all)
  return(c(list(h_s_pareto), list(sample_x_pareto_all)))
}


sampling_X_new_randomForest = function(num, X, Z, l, cb){
  D = ncol(X)
  x = c()
  s = c()
  i=0
  
  X_train = as.data.frame(X)
  
  z_factor = factor(Z, levels = c(FALSE, TRUE))
  rf_model <- randomForest(z_factor ~ ., data = X_train[, cb, drop = FALSE], 
                           importance = TRUE, ntree = 500)
  
  while ((nrow(x)<num)||is.null(x)){
    #print(i)
    i = i+1
    x_new = rmvnorm(1,rep(0,D))
    x_new = as.data.frame(x_new)
    # colnames(x_new) = c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10")
    
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
