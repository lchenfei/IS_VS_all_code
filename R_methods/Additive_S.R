
source('AMISE_bivariate.R')

S_add = function(x_mat, X_mat, Z, h_mat, combo, weight){
  S_est = 0
  for (c in 1:nrow(combo)){
    cb = combo[c,]
    h = h_mat[c,]
    w = weight[c]
    S_temp = S_bi(matrix(x_mat[,cb],ncol=2),X_mat[,cb],Z,h[1],h[2])
    S_est = S_est+w*S_temp
  }
  S = S_est
  S[is.na(S)] = 0
  S[S==Inf] = 0
  return(S)
}

S_add2 = function(x_mat, X_mat, Z, h_mat, combo, weight, i_list){
  S_est = 0
  for (i in 1:length(i_list)){
    c = i_list[i]
    #print(c)
    cb = combo[c,]
    h = h_mat[c,]
    w = weight[c]
    #print(w)
    S_temp = S_bi(matrix(x_mat[,cb],ncol=2),X_mat[,cb],Z,h[1],h[2])
    S_est = S_est+w*S_temp
  }
  S = S_est
  S[is.na(S)] = 0
  S[S==Inf] = 0
  return(S)
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

S_find_var = function(weight, z, S_all, ls){
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
  # sum_ce = -sum(z * log(S_predicted) + (1 - z) * log(1 - S_predicted), na.rm = T)
  return(var(S_predicted))
}

S_find_var_difference = function(weight, z, S_all, ls_1, ls_2){
  if (length(ls_2) == 0){
    return(S_find_var(weight, z, S_all, ls_1))
  }
  S_difference = S_find_var(weight, z, S_all, ls_1) - S_find_var(weight, z, S_all, ls_2)
  return(S_difference)
}

index_which = function(weight, index){
  n = length(weight)
  m = 0
  for(i in 1:n){
    index_i = t(combn(seq(1:n), i))
    m = m + nrow(index_i)
    if(index <= m){
      ncount = m 
      ncount = ncount - nrow(index_i)
      return(index_i[index - ncount, ])
    }
  }
}

S_find_h_s_ce = function(weight, z, S_all){
  n = length(weight)
  i_all = c()
  ce_all = c()
  count = 0
  for(i in 1:n){
    index = t(combn(seq(1:n), i))
    ce_i = c()
    ls_i = c()
    for(j in 1:nrow(index)){
      ls = index[j, ]
      count = count + 1
      ce_i = c(ce_i, S_find_ce(weight, z, S_all, ls))
      ls_i = c(ls_i, count)
    }
    ce_all = c(ce_all, ce_i)
    i_all = c(i_all, ls_i)
  }
  sorted_ce_10 = sort(ce_all)[1:10]
  order_ce = order(ce_all)
  count_index = i_all[order_ce]
  min_index = c()
  for(i in 1:10){
    min_index = c(min_index, list(index_which(weight, count_index[i])))
  }
  
  #smallest CE
  #p_i = which.min(ce_all)
  #i_index = i_all[p_i]
  #min_index = t(combn(seq(1:n), p_i))[i_index, ]
  return(min_index)
}

S_find_h_s_small_ce = function(weight, z, S_all){
  n = length(weight)
  i_all = c()
  ce_all = c()
  count = 0
  for(i in 1:n){
    index = t(combn(seq(1:n), i))
    ce_i = c()
    ls_i = c()
    for(j in 1:nrow(index)){
      ls = index[j, ]
      count = count + 1
      ce_i = c(ce_i, S_find_ce(weight, z, S_all, ls))
      ls_i = c(ls_i, count)
    }
    ce_all = c(ce_all, ce_i)
    i_all = c(i_all, ls_i)
  }
  sorted_ce_10 = sort(ce_all)[1:10]
  order_ce = order(ce_all)
  count_index = i_all[order_ce]
  min_index = c()
  for(i in 1:10){
    min_index = c(min_index, list(index_which(weight, count_index[i])))
  }
  
  # smallest CE
  p_i = which.min(ce_all)
  i_index = i_all[p_i]
  min_index = t(combn(seq(1:n), p_i))[i_index, ]
  return(min_index)
}

S_find_h_s_ce_order = function(weight, z, S_all){
  n = length(weight)
  i_all = c()
  ce_all = c()
  count = 0
  for(i in 1:n){
    index = t(combn(seq(1:n), i))
    ce_i = c()
    ls_i = c()
    for(j in 1:nrow(index)){
      ls = index[j, ]
      count = count + 1
      ce_i = c(ce_i, S_find_ce(weight, z, S_all, ls))
      ls_i = c(ls_i, count)
    }
    ce_all = c(ce_all, ce_i)
    i_all = c(i_all, ls_i)
  }
  sorted_ce = sort(ce_all)#[1:10]
  order_ce = order(ce_all)
  count_index = i_all[order_ce]
  min_index = c()
  for(i in 1:length(sorted_ce)){
    min_index = c(min_index, list(index_which(weight, count_index[i])))
  }
  
  #smallest CE
  #p_i = which.min(ce_all)
  #i_index = i_all[p_i]
  #min_index = t(combn(seq(1:n), p_i))[i_index, ]
  return(min_index)
}

in_or_not = function(l1, l2){
  for(i in 1:length(l1)){
    j = l1[i]
    if( j %in% l2 == FALSE){
      return(FALSE)
    } 
  }
  return(TRUE)
}

find_estimator = function(n, ce_all, i_all){
  order_all = order(ce_all)
  i_list = c()
  for(number in 1:n){
    i = order_all[number]
    i_list = c(i_list, list(i_all[[i]]))
  }
  return(i_list)
}

S_find_h_s_ce_order_greedy = function(weight, z, S_all){
  n = length(weight)
  i_all = c()
  ce_all = c()
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
        if(in_or_not(j, i_a) == F){
          in_or_not_j = c(in_or_not_j, j)
          ls = c(i_a, j)
          #print(ls)
          ce_i = c(ce_i, S_find_ce(weight, z, S_all, ls))
        }
      }
      #print(ce_i)
      ce_update =  min(ce_i)
      ce_all = c(ce_all, ce_update)
      in_index_update = in_or_not_j[which.min(ce_i)]
      i_a_update = c(i_a, in_index_update) #index[in_index_update, ]
      i_all = c(i_all, list(i_a_update))
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
  i_list = find_estimator(45, ce_all, i_all)
  return(i_list)
}

S_find_h_s_var_order = function(weight, z, S_all){
  n = length(weight)
  i_all = c()
  ce_all = c()
  count = 0
  for(i in 1:n){
    index = t(combn(seq(1:n), i))
    ce_i = c()
    ls_i = c()
    for(j in 1:nrow(index)){
      ls = index[j, ]
      count = count + 1
      ce_i = c(ce_i, S_find_var(weight, z, S_all, ls))
      ls_i = c(ls_i, count)
    }
    ce_all = c(ce_all, ce_i)
    i_all = c(i_all, ls_i)
  }
  sorted_ce = sort(ce_all, decreasing = TRUE)#[1:10]
  order_ce = order(ce_all, decreasing = TRUE)
  count_index = i_all[order_ce]
  min_index = c()
  for(i in 1:length(sorted_ce)){
    min_index = c(min_index, list(index_which(weight, count_index[i])))
  }
  
  #smallest CE
  #p_i = which.min(ce_all)
  #i_index = i_all[p_i]
  #min_index = t(combn(seq(1:n), p_i))[i_index, ]
  return(min_index)
}

S_find_h_s_shapley_order = function(weight, z, S_all){
  n = length(weight)
  i_all = c()
  ce_all = c()
  count = 0
  ls_previous = c()
  
  ls = c()
  
  left_element = 1:n
  
  ls_list = c()

  while(length(left_element) > 0){
  for(i in 1:n){
    S_difference_list = c()
    if (!(i %in%  ls)){
      ls_add = c(ls, i)
      S_difference = S_find_var_difference(weight, z, S_all, ls_add, ls_previous)
      S_difference_list = c(S_difference_list, S_difference)
    }
  }
  index_min = which.min(S_difference_list)
  index_min = left_element[index_min]
  ls = c(ls, index_min)
  ls_previous = ls
  left_element = left_element[!left_element == index_min]
  ls_list = c(ls_list, list(ls))
  }
  #smallest CE
  #p_i = which.min(ce_all)
  #i_index = i_all[p_i]
  #min_index = t(combn(seq(1:n), p_i))[i_index, ]
  return(ls_list)
}

find_ce_order_greedy = function(num, X, z, h_mat, combo, weight, S_all){
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
      if(ce_update > min(ce_all)){
        break
      }
      else{
      ce_all = c(ce_all, ce_update)
      in_index_update = in_or_not_j[which.min(ce_i)]
      i_a_update = c(i_a, in_index_update) #index[in_index_update, ]
      i_all = c(i_all, list(i_a_update))
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
  
  min_index = which.min(ce_all)
  h_s_pareto = i_all[[min_index]]
  # i_list = find_estimator(45, ce_all, i_all)
  return(h_s_pareto)
}
