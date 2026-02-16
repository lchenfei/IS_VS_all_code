S_find_index_validation <- function(X, Z, h_new_mat, v_list_all, weight, Nm) {
  m <- length(v_list_all)
  est_cell <- matrix(0, 10, m)
  est_validation <- matrix(0, 10, m)
  index_all <- numeric()
  n_per_step <- Nm
  for (k in 1:10) {
    randomIndices <- sample(dim(X)[1])
    train_indices <- randomIndices[1:((dim(X)[1]) / 2)]
    combo <- t(combn(seq(1:dim(X)[2]), 2))
    validation_indices <- setdiff(randomIndices, train_indices)
    X_train <- X[train_indices, ]
    X_validation <- X[validation_indices, ]
    Z_train <- Z[train_indices]
    Z_test <- Z[validation_indices]
    # Decide weights
    # cost_vec <- numeric()
    # S_train <- matrix(0, length(Z_train), 10)
    # for (i in 1:10) {
    #   cb <- combo[,i]
    #   S <- S_bi(X_train[,cb], X_train[,cb], Z_train, h_new_mat[i,1], h_new_mat[i,2])
    #   cost <- -sum(Z_train * log(S) + (1 - Z_train) * log(1 - S), na.rm = TRUE)
    #   cost_vec <- c(cost_vec, cost)
    #   S_train[, i] <- S
    # }
    # weight <- 1 / cost_vec / sum(1 / cost_vec)
    # list_values <- S_find_v_list_ce_value_all(weight, Z_train, S_train)
    # s_list_train <- list_values$s_list
    # v_list_train <- list_values$v_list
    # index_all <- c(index_all, v_list_train)
    #print(X_train)
    #print(Z_train)
    v_n <- 0
    
    for (v_list in v_list_all) {
      v_n <- v_n + 1
      weight_v <- weight / sum(weight[v_list])
      #print(weight_v)
      X_new <- sampling_X_new(n_per_step, X_train, Z_train, h_new_mat, combo, weight_v, v_list)
      #print(X_new)
      # index = runif(n_per_step/2, 1, n_per_step)
      # X_sample = X_new[index, ]
      
      # S_sample <- S_add2(X_sample, X_train, Z_train, h_new_mat, combo, weight_v, v_list)
      #print(S_new)
      # S_other <- S_add2(X_sample, X_train, Z_train, h_new_mat, combo, weight_v, v_list)
      # H_sample <- sqrt(S_sample)
      
      
      S_new <- S_add2(X_new, X_train, Z_train, h_new_mat, combo, weight_v, v_list)
      #print(S_new)
      S_test <- S_add2(X_new, X_validation, Z_test, h_new_mat, combo, weight_v, v_list)
      H_new <- sqrt(S_new)
      
      est_current <- sum((1 / H_new) / sum(1 / H_new) * S_new)
      est_cell[k, v_n] <- est_current
      
      est_all <- sum((1 / H_new) / sum(1 / H_new) * S_test)
      est_validation[k, v_n] <- est_all
      
      print(est_cell)
    }
  }
  std_train <- apply(est_cell, 2, sd) 
  
  print(std_train)
  
  # i = 1
  # index_list = c()
  # while(i <= m){
  #   k = 0
  #   for(j in 1:i){
  #     if(std_train[j] > std_train[i]){
  #       k = k + 1
  #     }
  #   }
  #   if(k == i - 1){
  #     index_list = c(index_list, i)
  #   }
  #   i = i + 1
  # }
  
  k_length = c()
  
  for(i in 1:length(v_list_all)){
    k_length = c(k_length, length(v_list_all[[i]]))
  }
  
  print(k_length)
  
  index_list = find_alpha(std_train, k_length)
  
  print("index list is that")
  print(index_list)
  
  std_est <- apply(est_validation, 2, sd) 
  
  # print(index_list)
  
  # print(std_est)
  
  index_min = which.min(std_est[index_list])
  
  index_test <- index_list[index_min]
  
  index <- v_list_all[index_test]

  
  return(c(index, list(index_list)))
}

find_alpha <- function(SE, B){
  l = c()
  for(i in 1:length(SE)){
    
      upper_all = c()
      lower_all = c()
      for(j in 1:length(SE)){
        if(B[i] > B[j] & i != j){
          upper_all = c(upper_all, (SE[i] - SE[j])/(B[j] - B[i]))
        }
        if(B[i] < B[j] & i != j){
          lower_all = c(lower_all, (SE[i] - SE[j])/(B[j] - B[i]))
        }
        if(B[i] == B[j] & i != j){
          upper_all = c(upper_all, (SE[j] - SE[i])/(B[i] - B[j]))
        }
        
        # print(lower_all)
        # print(upper_all)
      }
      if(length(upper_all) == 0){
        upper_bound = Inf
      }else{
        upper_bound = min(upper_all)
      }
      
      if(length(lower_all) == 0){
        lower_bound = 0
      }else{
        lower_bound = max(max(lower_all), 0)
      }
      
      
      
    
    
    print(lower_bound)
    print(upper_bound)
    
    if(upper_bound > lower_bound){
      l  = c(l, i)
    }
    
  }
  print(l)
  return (l)
}

S_find_index_validation_GP <- function(X, Z, v_list_all, l, Nm) {
  m <- length(v_list_all)
  est_cell <- matrix(0, 10, m)
  est_validation <- matrix(0, 10, m)
  index_all <- numeric()
  n_per_step <- Nm
  for (k in 1:10) {
    randomIndices <- sample(dim(X)[1])
    train_indices <- randomIndices[1:((dim(X)[1]) / 2)]
    combo <- t(combn(seq(1:dim(X)[2]), 2))
    validation_indices <- setdiff(randomIndices, train_indices)
    X_train <- X[train_indices, ]
    X_validation <- X[validation_indices, ]
    Z_train <- Z[train_indices]
    Z_test <- Z[validation_indices]
    # Decide weights
    # cost_vec <- numeric()
    # S_train <- matrix(0, length(Z_train), 10)
    # for (i in 1:10) {
    #   cb <- combo[,i]
    #   S <- S_bi(X_train[,cb], X_train[,cb], Z_train, h_new_mat[i,1], h_new_mat[i,2])
    #   cost <- -sum(Z_train * log(S) + (1 - Z_train) * log(1 - S), na.rm = TRUE)
    #   cost_vec <- c(cost_vec, cost)
    #   S_train[, i] <- S
    # }
    # weight <- 1 / cost_vec / sum(1 / cost_vec)
    # list_values <- S_find_v_list_ce_value_all(weight, Z_train, S_train)
    # s_list_train <- list_values$s_list
    # v_list_train <- list_values$v_list
    # index_all <- c(index_all, v_list_train)
    print(length(X_train))
    print(Z_train)
    v_n <- 0
    
    for (v_list in v_list_all) {
      v_n <- v_n + 1
      
      # (num, X, Y, l, cb)
      
      print(v_list)
      
      print(X_train[1, ])

      new <- sampling_X_new_GP(n_per_step, X_train, Z_train, l, v_list)

      X_new = new[[1]]
      
      S_new <- new[[2]] # S_add2(X_new, X_train, Z_train, h_new_mat, combo, weight_v, v_list)
      
      print("x_new")
      
      X_train_1 = as.matrix(X_train[, v_list])
      print(X_train)
      gp_model <- GauPro(X_train_1, Z_train)
      train_model <- km(formula = ~1,     # Constant mean
                        design = X_train_1, # Input variables
                        response = Z_train, # Output variable
                        covtype = "gauss", nugget.estim = TRUE)  # Kernel type
      
      print(X_new[, v_list])
      
      Y_predict = predict(train_model, matrix(X_new[, v_list], nrow = dim(X_new)[[1]]), type = "UK")
      
      mu = Y_predict$mean
      sigma = Y_predict$sd
      
      # print(mu)
      # print(sigma)

      S_test <- pnorm((mu - l) / sigma)  # S_add2(X_new, X_validation, Z_test, h_new_mat, combo, weight_v, v_list)
      H_new <- sqrt(S_new)
      
      est_current <- sum((1 / H_new) / sum(1 / H_new) * S_new)
      est_cell[k, v_n] <- est_current
      
      est_all <- sum((1 / H_new) / sum(1 / H_new) * S_test)
      est_validation[k, v_n] <- est_all
      
      print(est_cell)
    }
  }
  std_train <- apply(est_cell, 2, sd) 
  
  print(std_train)
  

  
  k_length = c()
  
  for(i in 1:length(v_list_all)){
    k_length = c(k_length, length(v_list_all[[i]]))
  }
  
  print(k_length)
  
  index_list = find_alpha(std_train, k_length)
  
  print("index list is that")
  print(index_list)
  
  std_est <- apply(est_validation, 2, sd) 
  
  # print(index_list)
  
  # print(std_est)
  
  index_min = which.min(std_est[index_list])
  
  index_test <- index_list[index_min]
  
  index <- v_list_all[index_test]
  
  
  return(c(index, list(index_list)))
}

S_find_index_validation_GP_Classification <- function(X, Z, v_list_all, l, Nm) {
  m <- length(v_list_all)
  est_cell <- matrix(0, 10, m)
  est_validation <- matrix(0, 10, m)
  index_all <- numeric()
  n_per_step <- Nm
  for (k in 1:10) {
    randomIndices <- sample(dim(X)[1])
    train_indices <- randomIndices[1:((dim(X)[1]) / 2)]
    combo <- t(combn(seq(1:dim(X)[2]), 2))
    validation_indices <- setdiff(randomIndices, train_indices)
    X_train <- X[train_indices, ]
    X_validation <- X[validation_indices, ]
    Z_train <- Z[train_indices]
    Z_test <- Z[validation_indices]
    # Decide weights
    # cost_vec <- numeric()
    # S_train <- matrix(0, length(Z_train), 10)
    # for (i in 1:10) {
    #   cb <- combo[,i]
    #   S <- S_bi(X_train[,cb], X_train[,cb], Z_train, h_new_mat[i,1], h_new_mat[i,2])
    #   cost <- -sum(Z_train * log(S) + (1 - Z_train) * log(1 - S), na.rm = TRUE)
    #   cost_vec <- c(cost_vec, cost)
    #   S_train[, i] <- S
    # }
    # weight <- 1 / cost_vec / sum(1 / cost_vec)
    # list_values <- S_find_v_list_ce_value_all(weight, Z_train, S_train)
    # s_list_train <- list_values$s_list
    # v_list_train <- list_values$v_list
    # index_all <- c(index_all, v_list_train)
    # print(length(X_train))
    # print(Z_train)
    v_n <- 0
    
    for (v_list in v_list_all) {
      v_n <- v_n + 1
      
      # (num, X, Y, l, cb)
      
      # print(v_list)
      
      # print(X_train[1, ])
      
      new <- sampling_X_new_GP_Classification(n_per_step, X_train, Z_train, l, v_list)
      
      X_new = new[[1]]
      
      S_new <- new[[2]] # S_add2(X_new, X_train, Z_train, h_new_mat, combo, weight_v, v_list)
      
      # print("x_new")
      
      gp_class_model <- gausspr(x = X_validation[, v_list], y = Z_test, kernel = "rbfdot")
      S_test <- predict(gp_class_model, matrix(X_new[, v_list], nrow = dim(X_new)[[1]]), type = "probabilities")
      
      # print("S new is")
      # print(S_new)
      # print(S_test)
     
      H_new <- sqrt(S_new)
      
      est_current <- sum((1 / H_new) / sum(1 / H_new) * S_new)
      est_cell[k, v_n] <- est_current
      
      est_all <- sum((1 / H_new) / sum(1 / H_new) * S_test)
      est_validation[k, v_n] <- est_all
      
      print(est_cell)
    }
  }
  std_train <- apply(est_cell, 2, sd) 
  
  print(std_train)
  
  
  
  k_length = c()
  
  for(i in 1:length(v_list_all)){
    k_length = c(k_length, length(v_list_all[[i]]))
  }
  
  # print(k_length)
  
  index_list = find_alpha(std_train, k_length)
  
  print("index list is that")
  print(index_list)
  
  std_est <- apply(est_validation, 2, sd) 
  
  # print(index_list)
  
  # print(std_est)
  
  index_min = which.min(std_est[index_list])
  
  index_test <- index_list[index_min]
  
  index <- v_list_all[index_test]
  
  
  return(c(index, list(index_list)))
}

S_find_index_validation_randomForest <- function(X, Z, v_list_all, l, Nm) {
  m <- length(v_list_all)
  est_cell <- matrix(0, 10, m)
  est_validation <- matrix(0, 10, m)
  index_all <- numeric()
  n_per_step <- Nm
  for (k in 1:10) {
    randomIndices <- sample(dim(X)[1])
    train_indices <- randomIndices[1:((dim(X)[1]) / 2)]
    combo <- t(combn(seq(1:dim(X)[2]), 2))
    validation_indices <- setdiff(randomIndices, train_indices)
    X_train <- X[train_indices, ]
    X_validation <- X[validation_indices, ]
    Z_train <- Z[train_indices]
    Z_test <- Z[validation_indices]

    v_n <- 0
    
    for (v_list in v_list_all) {
      v_n <- v_n + 1
      
      new <- sampling_X_new_randomForest(n_per_step, X_train, Z_train, l, v_list)
      
      X_new = new[[1]]
      
      S_new <- new[[2]] # S_add2(X_new, X_train, Z_train, h_new_mat, combo, weight_v, v_list)
      
      X_new = as.data.frame(X_new)
      
      X_validation = as.data.frame(X_validation)
      
      z_test_factor = factor(Z_test, levels = c(FALSE, TRUE))
      
      rf_model <- randomForest(z_test_factor ~ ., data = X_validation[, v_list, drop = FALSE], 
                               importance = TRUE, ntree = 500)
      
      S_test <- predict(rf_model, X_new[, v_list, drop = FALSE], type = "prob")[, "TRUE"]
      
      H_new <- sqrt(S_new)
      
      est_current <- sum((1 / H_new) / sum(1 / H_new) * S_new)
      est_cell[k, v_n] <- est_current
      
      est_all <- sum((1 / H_new) / sum(1 / H_new) * S_test)
      est_validation[k, v_n] <- est_all
      
      print(est_cell)
    }
  }
  std_train <- apply(est_cell, 2, sd) 
  
  print(std_train)
  
  
  
  k_length = c()
  
  for(i in 1:length(v_list_all)){
    k_length = c(k_length, length(v_list_all[[i]]))
  }
  
  # print(k_length)
  
  index_list = find_alpha(std_train, k_length)
  
  print("index list is that")
  print(index_list)
  
  std_est <- apply(est_validation, 2, sd) 
  
  # print(index_list)
  
  # print(std_est)
  
  index_min = which.min(std_est[index_list])
  
  index_test <- index_list[index_min]
  
  index <- v_list_all[index_test]
  
  
  return(c(index, list(index_list)))
}


find_alpha <- function(SE, B){
  l = c()
  for(i in 1:length(SE)){
    
    upper_all = c()
    lower_all = c()
    for(j in 1:length(SE)){
      if(B[i] > B[j] & i != j){
        upper_all = c(upper_all, (SE[i] - SE[j])/(B[j] - B[i]))
      }
      if(B[i] < B[j] & i != j){
        lower_all = c(lower_all, (SE[i] - SE[j])/(B[j] - B[i]))
      }
      if(B[i] == B[j] & i != j){
        upper_all = c(upper_all, (SE[j] - SE[i])/(B[i] - B[j]))
      }
      
      # print(lower_all)
      # print(upper_all)
    }
    if(length(upper_all) == 0){
      upper_bound = Inf
    }else{
      upper_bound = min(upper_all)
    }
    
    if(length(lower_all) == 0){
      lower_bound = 0
    }else{
      lower_bound = max(max(lower_all), 0)
    }
    
    
    
    
    
    print(lower_bound)
    print(upper_bound)
    
    if(upper_bound > lower_bound){
      l  = c(l, i)
    }
    
  }
  print(l)
  return (l)
}