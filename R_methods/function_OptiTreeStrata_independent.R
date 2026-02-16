## This code is to generate Figure 4 in the manuscript
# Function to check and install missing packages
install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

required_packages <- c("truncnorm", "mvtnorm", "caret", 'rootSolve','stats','lhs','LHD','R.utils','crch')
install_if_missing(required_packages)
library(truncnorm)
library(mvtnorm)
library(caret)
library(rootSolve)
library(stats)
library(lhs)
library(LHD)
library(R.utils)
library(crch)  # For truncated T
source('AMISE_univariate.R')


#########################################################################
#1. Example functions
#########################################################################
## Marginal X distribution


x_dist <- function(x){
  dist = switch(case,
                '1' = dnorm(x,0, 1),
                '2' = dnorm(x,0, 1),
                '3' = dnorm(x,0, 1),
                'Heterogeneous' = dnorm(x,0, 1),
                'high_dimension' = dnorm(x,0, 1),
                'Gamma_1' = dgamma(x,5,5),
                'Gamma_2' = dgamma(x,5,5),
                'Gamma_3' = dgamma(x,5,5),
                't_trun_1' = dtt(x,df=1,left=min_x,right=max_x),
                't_trun_2' = dtt(x,df=1,left=min_x,right=max_x),
                't_trun_3' = dtt(x,df=1,left=min_x,right=max_x))
  return(dist)
}

## Make Y depending on examples
Make_Y <- function(X){
  Y = switch(case,
             "1" = rnorm(nrow(X), mean = mu_condi1(X),sd = 1),
             "2" = rnorm(nrow(X), mean = mu_condi2(X),sd = 1),
             "3" = rnorm(nrow(X), mean = mu_condi3(X),sd = 1),
             "4" = rnorm(nrow(X), mean = mu_condi4(X),sd = 1),
             'Heterogeneous' = rnorm(nrow(X), sd = abs(X[,1])),
             'high_dimension' = rnorm(nrow(X), mean = mu_condi_high(X),sd = 1),
             "Gamma_1" = rnorm(nrow(X), mean = mu_condi1(X),sd = 1),
             "Gamma_2" = rnorm(nrow(X), mean = mu_condi2(X),sd = 1),
             "Gamma_3" = rnorm(nrow(X), mean = mu_condi3(X),sd = 1),
             "t_trun_1" = rnorm(nrow(X), mean = mu_condi1(X),sd = 1),
             "t_trun_2" = rnorm(nrow(X), mean = mu_condi2(X),sd = 1),
             "t_trun_3" = rnorm(nrow(X), mean = mu_condi3(X),sd = 1)) 
  return(Y)
}

## Make Z depending on the interests
Make_Z <- function(Y,l){
  Z = switch(Interest,
             "Failure_prob" = ifelse(Y > l,1,0),
             "Shortfall" = ifelse(Y > l, Y,0),
             'Expectation' = Y)
  
  return(Z)
}

## Y|X ~ normal(mu_condi, sigma_condi)
#example 1
mu_condi1 = function(X){
  65 - 40 * exp(-0.2 * sqrt((X[,1]^2 + X[,2]^2)/2)) - 20 * exp(-0.2 * abs(X[,1])) -
    5 * exp(-0.2 * sqrt((X[,2]^2 + X[,3]^2)/2)) -
    exp(cos(2*pi*X[,1]*X[,2])) - exp(cos(2*pi*X[,2]*X[,3])) - exp(cos(2*pi*X[,3]*X[,1])) -
    exp(cos(2*pi*X[,1]*X[,2]*X[,3]))
}
sigma_condi = 1

#example 2
mu_condi2 = function(X){
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


#example 3 
mu_condi3 = function(X){
  X1 = X[,1];  X2 = X[,2];  X3 = X[,3];  X4 = X[,4];  X5 = X[,5] 
  X6 = X[,6];  X7 = X[,7];  X8 = X[,8];  X9 = X[,9];  X10 = X[,10]
  m = 40*(1-exp(-0.2*sqrt((X1^2+X2^2)/2)))+
    20*(1-exp(-0.2*abs(X1)))+
    5*(1-exp(-0.2*sqrt((X2^2+X3^2)/2)))- 0.1*exp(-0.2 * sqrt((X4^2+X5^2)/2)) -
    0.1*exp(-0.2 * sqrt((X6^2+X7^2)/2)) - 0.1*exp(-0.2 * sqrt((X8^2+X9^2+X10^2)/3)) - 
    exp(cos(2*pi*X1*X3))-
    exp(cos(2*pi*X2*X3))-
    exp(cos(2*pi*X1*X2))
  return(m)
}

#example 4 
mu_condi4 = function(X){
  X1 = X[,1];  X2 = X[,2];  X3 = X[,3];  X4 = X[,4];  X5 = X[,5] 
  m = 65 - 40 * exp(-0.2 * sqrt((X1^2 + X2^2)/2)) - 
    20 * exp(-0.2 * abs(X1)) - 
    5 * exp(-0.2 * sqrt((X2^2 + X3^2)/2))  - 
    0.1 * exp(-0.2 * sqrt((X4^2 + X5^2)/2)) - 
    exp(cos(2 * pi * X1 * X2)) - exp(cos(2 * pi * X1 * X3)) - exp(cos(2 * pi * X2 * X3))
  return(m)
}

#high-dimension
mu_condi_high = function(X){
  X1 = X[,1];  X2 = X[,2];  X3 = X[,3];  X4 = X[,4];  X5 = X[,5] 
  X6 = X[,6];  X7 = X[,7];  X8 = X[,8];  X9 = X[,9];  X10 = X[,10]
  m = 40*(1-exp(-0.2*sqrt((X1^2+X2^2)/2)))+
    20*(1-exp(-0.2*abs(X1)))+
    5*(1-exp(-0.2*sqrt((X2^2+X3^2)/2)))- 0.1*exp(-0.2 * sqrt((X4^2+X5^2)/2)) -
    0.1*exp(-0.2 * sqrt((X6^2+X7^2)/2)) - 0.1*exp(-0.2 * sqrt((X8^2+X9^2+X10^2)/3)) - 
    rowSums(0.01*exp(-0.2*abs(X[,11:50]))) -
    exp(cos(2*pi*X1*X3))-
    exp(cos(2*pi*X2*X3))-
    exp(cos(2*pi*X1*X2))
  return(m)
}

# Heterogeneous
Heterogeneous  = function(x){
  X1 = X[,1];  X2 = X[,2];  X3 = X[,3];
  return(X1^2 + 0.1*X2 + 0.1*X3)
}


#########################################################################
#2. Tree-based Stratified Sampling
#########################################################################

## Draw Pilot_sample sample and make data.frame##
Draw_Pilot_sample = function(n_allocation, p, min_x, max_x){
  X <- matrix(runif(n_allocation * p,min_x,max_x),ncol=p)
  Y = Make_Y(X)
  sample = data.frame(X=X, Y=Y)
  if(p == 1){colnames(sample) = c('X.1','Y')}
  return(sample)
}  

## Draw Pilot_sample using LHS
Draw_Pilot_sample_LHS = function(n_allocation, p, min_x, max_x){
  X = matrix(0, ncol = p, nrow = n_allocation)
  for (j in 1:p){
    X[,j] = qunif(randomLHS(n_allocation, 1),min_x,max_x)
  }
  Y = Make_Y(X)
  sample = data.frame(X=X, Y=Y)
  if(p == 1){colnames(sample) = c('X.1','Y')}
  return(sample)
}

## Draw Pilot_sample using OA-LHS
Draw_Pilot_sample_OA_LHS = function(n_allocation, p, min_x, max_x){
  prime_number = switch(as.character(n_allocation),'100' = 5, '170' = 13,'200' = 7, '340' = 13,'500' = 11, '1000' = 11)
  lambda = switch(as.character(n_allocation),'100' = 4,'170' = 1 ,'200' = 4,'340' = 2, '500' = 4, '1000' = 8)
  OA = create_OA(prime_number,lambda,p)
  D <- OA2LHD(OA)
  col_shuff <- sample(1:p)
  X = matrix(0, ncol = p, nrow = nrow(D))
  for (j in col_shuff){
    X[,j] = qunif(D[,j]/nrow(D) - runif(nrow(D))/nrow(D),min_x,max_x)
  }
  Y = Make_Y(X)
  sample = data.frame(X=X, Y=Y)
  if(p == 1){colnames(sample) = c('X.1','Y')}
  return(sample)
}  

create_OA = function(q, lambda,columns){
  Temp = vector(mode = "list", length = lambda)
  OA = matrix(ncol=columns)
  for(i in 1:lambda){
    Temp[[i]] = createBose(q, columns, bRandom = TRUE) +1
    OA = rbind(OA,Temp[[i]]) 
  }
  OA = OA[-1,]
  return(OA)
}



## Decide the bandwidth by minimizing AMISE
find_bandwidth = function(Data){
  P = dim(Data)[2] - 2
  h_hat = rep(0,P)
  if(Interest == 'Expectation'){
    for(p in 1:P){
      h_hat[p] = bw.nrd(Data[,p])
    }
  }else{
    for(p in 1:P){
      upper = range(Data[,p])[2] - range(Data[,p])[1]
      interval = cut(c(0.01, upper/10),breaks = 2)
      interval = levels(interval)
      result_AMISE = sapply(interval,function(i){
        lower = as.numeric(substring(i,2,unlist(gregexpr(',',i))-1))
        upper = as.numeric(substring(i,unlist(gregexpr(',',i))+1,unlist(gregexpr(']',i))-1))
        return((optimize(function(x) AMISE(Data[,p],Data$Z,x,min(Data[,p]),max(Data[,p])),interval = c(lower,upper))))
      })
      h_hat[p] = as.numeric(result_AMISE['minimum',][which.min(unlist(result_AMISE['objective',]))])
    }
  }
  return(h_hat)
}

find_bandwidth_iteration = function(Data,h_hat){
  P = dim(Data)[2] - 2
  limit_time = ifelse(P >=10, 600, 120)
  h_hat_iter = rep(0,P)
  if(Interest == 'Expectation'){
    for(p in 1:P){
      h_hat_iter[p] = bw.nrd(Data[,p])
    }
  }else{
    for(p in 1:P){
      Upper = range(Data[,p])[2] - range(Data[,p])[1]
      lower = h_hat[p] * 0.75
      upper = ifelse(h_hat[p] * 1.25 < Upper, h_hat[p] * 1.25, Upper)
      #result_AMISE = optimize(function(x) AMISE(Data[,p],Data$Z,x,min(Data[,p]),max(Data[,p])),interval = c(lower,upper))
      Optimize_AMISE_limit_time = tryCatch({
        withTimeout(optimize(function(x) AMISE(Data[,p],Data$Z,x,min(Data[,p]),max(Data[,p])),interval = c(lower,upper)), timeout = limit_time)
      }, error = function(e){
        if(grepl("reached elapsed time limit",e$message)){
          h_hat[p]
        }else{h_hat[p]}
      })
      if(is.list(Optimize_AMISE_limit_time) == TRUE){h_hat_iter[p] = Optimize_AMISE_limit_time$minimum
      }else{h_hat_iter[p] = Optimize_AMISE_limit_time}
    }
  }
  return(h_hat_iter)
}
# If optimization takes too long, stop and return the previous value
interruptor <- function(FUN,args, time.limit, alternative){
  
  results <- 
    tryCatch({
      withTimeout({FUN(args)}, timeout=time.limit)
    }, error = function(e){
      if(grepl("reached elapsed time limit",e$message)){
        alternative
      }else{paste(e$message,"EXTRACTERROR")}
    })
  
  if(grepl("EXTRACTERROR",results)){
    print(gsub("EXTRACTERROR","",results))
    results <- NULL
  } 
  
  return(results)
}   

## Calculate each stratum's values for univariate input variable
Calculate_strata_value = function(X,Z,h_hat,a){
  # X should be a vector
  
  # fit the Gaussian kernel regression
  hx = function(x){S(x,X,Z,h_hat)} #Gaussian kernel from AMISE.r / h_hat means the optimal width of gaussian kernel
  Expectation_stratum = function(x){hx(x) * x_dist(x)}
  tx = function(x){S(x,X,Z^2,h_hat)}
  Expectation2_stratum = function(x){tx(x) * x_dist(x)}
  
  p_strata = c(); mu_strata = c(); sd_strata = c();s_strata = c();
  strata_break = cut(X,a)
  strata = split(X,strata_break)
  I = length(a)-1
  
  for(i in 1:I){
    p_strata[i]  = as.numeric(integrate(x_dist,a[i],a[i+1])[1]) #find Pi from x_dist
    mu_strata[i] = as.numeric(integrate(Expectation_stratum,a[i],a[i+1])[1]) / p_strata[i] #Using kernel density
    if(Interest == 'Failure_prob'){
      s_strata[i] = mu_strata[i]
    }else{
      s_strata[i] = as.numeric(integrate(Expectation2_stratum,a[i],a[i+1])[1]) / p_strata[i] #Using kernel density
    }
    sd_strata[i] = sqrt(s_strata[i] - mu_strata[i]^2)
    p_strata[i]  = p_strata[i] / as.numeric(integrate(x_dist,a[1],a[length(a)])[1]) #find conditional Pi 
  }
  k = 1 / sd_strata
  k[k == Inf | k == - Inf] = 100000 #for preventing errors
  return(data.frame(p_strata,mu_strata, s_strata, sd_strata,k))
}

## Function to find a stratum boundary
find_range = function(temp,split_num,v){
  interval = as.character(temp[, dim(temp)[2]-split_num+v][1])
  lower = as.numeric(substring(interval,2,unlist(gregexpr(',',interval))-1))
  upper = as.numeric(substring(interval,unlist(gregexpr(',',interval))+1,unlist(gregexpr(']',interval))-1))
  return(c(lower,upper))
}

## Find breakpoints 
Find_a = function(X,Z,h_hat,a){
  a_list = list(); a_list[[1]] = a; obj_list = list()
  # Kernel regression surrogates for h(x) and t(x)
  # h_hat means the optimal width of gaussian kernel
  hx = function(x){S(x,X,Z,h_hat)} # surrogate for h(x) in Theorem 1"
  tx = function(x){S(x,X,Z^2,h_hat)} # surrogate for t(x) in Theorem 1" 
  
  for(i in 1:50){
    strata_value = Calculate_strata_value(X,Z,h_hat,a)
    p_strata = strata_value[,1];  mu_strata = strata_value[,2]
    s_strata = strata_value[,3];  sd_strata = strata_value[,4]
    k = strata_value[,5];    I = length(a)-1
    obj_list[[i]] = sum(p_strata * sd_strata) # checking if the algorithm is working correctly
    for(i in 1:(I-1)){
      equation <- function(x) {tx(x) * (k[i]-k[i+1]) - 2 * hx(x) * (k[i]*mu_strata[i]-k[i+1]*mu_strata[i+1]) + (k[i]*s_strata[i] - k[i+1]*s_strata[i+1])}
      All <- uniroot.all(equation, c(min(X), max(X)))
      if(sum(All) == 0){
        a[i+1] = a[i+1]
      }else{
        a[i+1] = All[which.min(abs(All - a[i+1]))]
      }
    }  
    a_list[[i+1]] = a
    if(sum(abs(a_list[[i+1]] - a_list[[i]])) < 0.001){
      return(list('p_strata'=p_strata,'mu_strata'=mu_strata,'s_strata'=s_strata,'sd_strata'=s_strata,'k'=k,'breaking_points'=a))
      break
    }
  }
  return(list('p_strata'=p_strata,'mu_strata'=mu_strata,'s_strata'=s_strata,'sd_strata'=s_strata,'k'=k,'breaking_points'=a))
}


## Calculate all input variables' stratum values
Search_variable = function(temp,Z,I,p){
  obj_value = c(); a_variable = list(); p_strata = list(); mu_strata = list(); sd_strata = list()
  for(v in 1:p){
    a_list = list(); obj_list = list()
    #Set initial breaking points
    #checking the variable already split
    if(sum(paste('split_',colnames(temp)[v],sep='') == colnames(temp))>0){
      variable = max(which((paste('split_',colnames(temp)[v],sep='') == colnames(temp))==TRUE))
      interval = as.character(temp[,variable][1])
      lower = as.numeric(substring(interval,2,unlist(gregexpr(',',interval))-1))
      upper = as.numeric(substring(interval,unlist(gregexpr(',',interval))+1,unlist(gregexpr(']',interval))-1))
      a = seq(lower,upper,length.out=I+1); a_list[[1]] = a # equi-distant
      
    }else{
      a = seq(min_x,max_x,length.out = I + 1); a_list[[1]] = a # min_x, max_x are provided from main file
    }
    
    result = Find_a(temp[,v],temp$Z,h_hat[v],a) # h_hat from entire input space
    obj_value[v] = sum(result[[1]] * result[[4]])
    p_strata[[v]] = result[[1]]; mu_strata[[v]] = result[[2]]; sd_strata[[v]] = result[[4]];
    a_variable[[v]] = sort(result[[6]])
  }
  summary_result = list('obj_value' = obj_value, 'breaking_points' = a_variable, 
                        'p_strata' = p_strata,'sd_strata' = sd_strata, 'mu_strata' = mu_strata )
  return(summary_result)
}



## Decide the splitting variable and split the stratum
Split_data = function(stratum_data,Z,I,p){
  if(nrow(stratum_data)<= 1){
    #print('The number of sample within the node is fewer than 1')
    return(list(stratum_data))
  }else if(mean(Z) == 0 | mean(Z) == 1 | sd(Z) == 0){
    #print('This node is already homogeneous')
    return(list(stratum_data))
  }else{
    result = Search_variable(stratum_data,Z,I,p)
    v = which.min(result[['obj_value']])
    a_v = result[['breaking_points']][which.min(result[['obj_value']])][[1]]
    if(sum(table(cut(stratum_data[,v],a_v,include.lowest = T)) == 0) >= 1){
      return(list(stratum_data))
    }
    stratum_data = cbind(stratum_data,cut(stratum_data[,v],a_v,include.lowest = T))
    colnames(stratum_data)[dim(stratum_data)[2]] = paste('split_X.',v,sep='')
    return(split(stratum_data,stratum_data[,dim(stratum_data)[2]]))
  }
}


## Calculate within multivariate Tree structure by using sample mean
Calculate_Mstrata_value = function(Tree,depth,l){
  p_strata = c(); mu_strata = c(); sd_strata = c();s_strata = c();
  for(i in 1:length(Tree[[depth]])){
    temp = Tree[[depth]][[i]]
    if(is.null(temp)){
      p_strata[i] = 0 ; mu_strata[i] = 0; 
      s_strata[i] = 0  ; sd_strata[i]= 0
    }else{
      p_temp = rep(1,dim(temp)[2])
      temp_Pilot = Tree[[1]][paste('X.',substring(colnames(temp)[dim(temp)[2]],9,9),sep='')][,1] #select variable
      n_num <- nrow(temp)
      
      ## finding the joint probability
      split_num <- dim(temp)[2]-dim(Tree[[1]])[2]
      for(ii in 1:split_num){  #number of split
        interval = as.character(temp[, dim(temp)[2]-split_num+ii][1])
        lower = as.numeric(substring(interval,2,unlist(gregexpr(',',interval))-1))
        upper = as.numeric(substring(interval,unlist(gregexpr(',',interval))+1,unlist(gregexpr(']',interval))-1))
        p_temp[dim(temp)[2]-split_num+ii] = integrate(x_dist,lower,upper)$value
      }
      ## Checking the same splitting variable
      if(sum(table(colnames(temp)) > 1) >= 1){
        repeated_variable = table(colnames(temp)) > 1
        repeated_variable = names(table(colnames(temp))[repeated_variable])
        for(ii in 1:length(repeated_variable)){
          repnum_variable = which(colnames(temp) == repeated_variable[ii])
          last_variable = max(repnum_variable)
          interval = as.character(temp[,last_variable][1])
          lower = as.numeric(substring(interval,2,unlist(gregexpr(',',interval))-1))
          upper = as.numeric(substring(interval,unlist(gregexpr(',',interval))+1,unlist(gregexpr(']',interval))-1))
          p_temp[last_variable] = integrate(x_dist,lower,upper)$value
          p_temp[repnum_variable[-length(repnum_variable)]] = 1
        }
      }
      p_strata[i] = prod(p_temp)  
      mu_strata[i] = mean(temp$Z) 
      s_strata[i] = mean((temp$Z)^2)
      sd_strata[i] = sqrt(n_num / (n_num-1) * (s_strata[i] - mu_strata[i]^2))
      if(mu_strata[i] == 0 | mu_strata[i] == 1){sd_strata[i] = 0}
      #if(n_num == 1){sd_strata[i] = 0}
    }
    
  }
  return(data.frame(p_strata,mu_strata, s_strata, sd_strata))
}



## Find Stratification using Tree structure
Tree_stratification = function(Data,n_simulation,nj_min,l,I,m,l_RR){
  
  ### Input
  alpha_candi = seq(0.000,0.01,length.out = 101) # set of alpha
  ## Divide Training and Validation sets. 
  N = nrow(Data); homogeneous = 1; 
  while(homogeneous == 1){
    Vali_num = sample(1:N, 0.3*N,replace = FALSE)
    homogeneous = ifelse(mean(Data[-Vali_num,]$Z) == 0,1,0)
  }
  Validation = Data[Vali_num,]; Data = Data[-Vali_num,]
  
  ### Initialization
  k = 0 # index
  RR = 1 # Variance Reduction rate
  Strata =  vector(mode = "list", length = 100) # To save tree strcture
  Strata[[k + 1]] = Data;
  
  
  split_num = 1; 
  obj = c(); # Save objective function value to calculate RR
  node_save = list(); # To name strata structure (node)
  Tree_p_save = list(); Tree_sd_save = list(); n_optimal_save = list()
  
  obj_vali_save = c() # To save obj value as splits go
  
  ###################################
  # Splitting strata
  ###################################
  ## Initialization for a splitting
  Strata[[split_num+1]] = append(Strata[[split_num+1]],Split_data(Strata[[split_num]],Strata[[split_num]]$Z,I,p))
  node = cbind(rep(split_num+1,length(Strata[[split_num+1]])), 1:length(Strata[[split_num+1]]))
  node_save[[split_num]] = node
  
  
  Tree_p = c(); Tree_sd = c(); 
  for(i in 1:nrow(node)){
    Strata_result = Calculate_Mstrata_value(Strata,node[i,1],l)
    Tree_p[i] = Strata_result[node[i,2],1]; Tree_sd[i] = Strata_result[node[i,2],4]
  }
  Tree_p_save[[split_num]] = Tree_p; Tree_sd_save[[split_num]] = Tree_sd
  obj[split_num] = sum(Tree_p * Tree_sd)
  #print(paste0('Split ',split_num,' times ','obj:',round(obj[split_num],5)))
  
  while((RR >= l_RR) & (split_num < floor(1 / l_RR - 1))){
    k = k + 1
    
    ## Identify splitting node and split
    temp_obj = c(); temp_candi = list()
    # conduct the univariate finding breaking points algorithm and calculate obj value
    for(ii in 1:nrow(node)){
      i_range = 1:nrow(node); i_range = i_range[-ii]
      temp_p_others = c(); temp_sd_others = c();
      for(i in i_range){
        temp_result = Calculate_Mstrata_value(Strata,node[i,1],l)
        temp_p_others[i] = temp_result[node[i,2],1] 
        temp_sd_others[i] = temp_result[node[i,2],4] 
      }
      temp_candi[[ii]] = Split_data(Strata[[node[ii,1]]][[node[ii,2]]],Strata[[node[ii,1]]][[node[ii,2]]]$Z,I,p)
      
      # for calculating strata values
      temp_split = vector(mode = "list", length = 100)
      temp_split[1:(node[ii,1]-1)] = Strata[1:(node[ii,1]-1)]
      temp_split[[node[ii,1]]] = append(temp_split[[node[ii,1]]],temp_candi[[ii]]) 
      temp_strata_value = Calculate_Mstrata_value(temp_split,node[ii,1],l)
      temp_p = temp_strata_value[,1]; temp_sd = temp_strata_value[,4]
      temp_obj[ii] = sum(temp_p * temp_sd) + sum(temp_p_others * temp_sd_others,na.rm = TRUE)
    }
    
    # choose the splitting node 
    split_nodenum = which.min(temp_obj)
    
    # update node information and Split
    if(length(temp_candi[[split_nodenum]]) == 1){ 
      #if splitting does not occur because the stratum is already homogeneous, stop growing a tree
      break
    }else{
      if(split_nodenum == nrow(node)){
        node = rbind(node[1:split_nodenum,], node[split_nodenum,])
      }else{
        node = rbind(node[1:split_nodenum,], node[split_nodenum,], node[(split_nodenum+1):nrow(node),])
      }
      node[split_nodenum:(split_nodenum+1),] = matrix(c(rep(node[split_nodenum,1] + 1,2), 2*node[split_nodenum,2]-1, 2*node[split_nodenum,2]),nrow=2)
      
      # Split the splitting node
      Strata[[node[split_nodenum,1]]][[node[split_nodenum,2]]] = temp_candi[[split_nodenum]][[1]]
      Strata[[node[(split_nodenum+1),1]]][[node[(split_nodenum+1),2]]] = temp_candi[[split_nodenum]][[2]]
      
      split_num = split_num + 1
      node_save[[split_num]] = node
    }
    
    # Calculate RR
    Tree_p = c(); Tree_sd = c(); 
    for(i in 1:nrow(node)){
      Strata_result = Calculate_Mstrata_value(Strata,node[i,1],l)
      Tree_p[i] = Strata_result[node[i,2],1]; Tree_sd[i] = Strata_result[node[i,2],4]
    }
    Tree_p_save[[split_num]] = Tree_p; Tree_sd_save[[split_num]] = Tree_sd
    obj[split_num] = sum(Tree_p * Tree_sd)
    RR = (obj[split_num-1] - obj[split_num])/obj[split_num-1]
    #print(paste0('Split ',split_num,' times ','obj:',round(obj[split_num],5),' RR:',round(RR,2)))
  }
  
  
  ###################################
  # Finding the best stratification
  ###################################
  
  ## Find candidate best sub-strata that minimize cost-complexity criterion
  Strata_alpha = c(); test = c()
  for(a in 1:100){
    Strata_alpha[a] = which.min(obj + alpha_candi[a]  * (2:(length(obj)+1)))
  }
  
  ## Assess the best sub-strata using customized Cross-validation
  for(S_alpha in unique(Strata_alpha)){
    ## Load node information
    node = node_save[[S_alpha]]
    
    ## Compute optimal allocation
    Tree_p = Tree_p_save[[S_alpha]]; Tree_sd = Tree_sd_save[[S_alpha]]; 
    if(sum(Tree_sd) == 0){ #When all nodes are homogeneous
      Best_strata_num = S_alpha
      alpha = min(alpha_candi[Strata_alpha == Best_strata_num])
      
      Best_strata = list()
      node = node_save[[Best_strata_num]]
      for(i in 1:nrow(node)){Best_strata[[i]] = Strata[[node[i,1]]][[node[i,2]]]}
      
      n_optimal = round(n_simulation / m * Tree_p)
      # Set minimum allocation to prevent sparse stratum
      lack_ind = n_optimal < nj_min  
      if( abs(sum(n_optimal[lack_ind] - nj_min)) >= 1){
        n_optimal = rep(nj_min, length(n_optimal))
        if(length(n_optimal) * nj_min < (n_simulation / m)){
          n_optimal[!lack_ind] = n_optimal[!lack_ind] + round(((n_simulation / m) - length(n_optimal) * nj_min) * (Tree_p[!lack_ind]) / sum(Tree_p[!lack_ind]))
        }
        
        # adjust rounding error by changing the maximum allocation
        if(sum(n_optimal) !=  round(n_simulation / m) ){
          n_optimal[which.max(n_optimal)] = n_optimal[which.max(n_optimal)]-(sum(n_optimal) -  n_simulation / m )
        }
      }
      
    }else{
      n_optimal = round(n_simulation / m * (Tree_p*Tree_sd) / sum(Tree_p*Tree_sd))
      # Set minimum allocation to prevent sparse stratum
      lack_ind = n_optimal < nj_min  
      if( abs(sum(n_optimal[lack_ind] - nj_min)) >= 1){
        n_optimal = rep(nj_min, length(n_optimal))
        if(length(n_optimal) * nj_min < (n_simulation / m)){
          n_optimal[!lack_ind] = n_optimal[!lack_ind] + round(((n_simulation / m) - length(n_optimal) * nj_min) * (Tree_p[!lack_ind]*Tree_sd[!lack_ind]) / sum(Tree_p[!lack_ind]*Tree_sd[!lack_ind]))
        }
        
        # adjust rounding error by changing the maximum allocation
        if(sum(n_optimal) !=  round(n_simulation / m) ){
          n_optimal[which.max(n_optimal)] = n_optimal[which.max(n_optimal)]-(sum(n_optimal) -  n_simulation / m )
        }
      }
    }
    
    n_optimal_save[[S_alpha]] = n_optimal
    
    ## Assess the best sub-strata using customized Cross-validation
    obj_vali_resampling = c();obj_vali_temp = matrix(rep(0, 10* dim(node)[1]),nrow = 10)
    for(i in 1:nrow(node)){
      # Select stratum
      Strata_resampling = Strata[[node[i,1]]][[node[i,2]]]; Strata_Validation = Validation
      Strata_split_num = dim(Strata_resampling)[2] - which(colnames(Strata_resampling) == 'Z') 
      
      # Find stratum boundaries
      variable_select = c(); variable_range = matrix(nrow=Strata_split_num,ncol = 2)
      for(ii in 1:Strata_split_num){
        variable_select[ii] = as.numeric(substring(colnames(Strata_resampling)[dim(Strata_resampling)[2]-Strata_split_num+ii],9,10))
        variable_range[ii,] = find_range(Strata_resampling,Strata_split_num,ii)
      }
      
      # Extract samples within the stratum
      for(ii in 1:Strata_split_num){
        Strata_Validation = Strata_Validation[!is.na(cut(Strata_Validation[,variable_select[ii]],breaks = variable_range[ii,])),]
      }
      ni = nrow(Strata_Validation)
      
      # Resampling 10 times and take average
      for(b in 1:10){
        if(ni <= 1){
          obj_vali_temp[b,i] = 0
        }else{
          boot_num = sample(1:ni,n_optimal[i],replace = TRUE)
          obj_vali_temp[b,i] = sd(Strata_Validation[boot_num,]$Z) * Tree_p[i]
        }
      }
    }
    
    obj_vali_resampling = rowSums(obj_vali_temp)
    obj_vali_save[S_alpha] = mean(obj_vali_resampling)
  }
  
  ## Choose the best substrata based on the validation result
  Best_strata_num = which.min(obj_vali_save)
  alpha = min(alpha_candi[Strata_alpha == Best_strata_num])
  #print(paste('alpha:',alpha))
  #print(paste('Full_strata size:',split_num))
  #print(paste('Best_strata:',Best_strata_num))
  
  ## Print result
  # Organize stratification from tree structure
  Best_strata = list()
  node = node_save[[Best_strata_num]]
  for(i in 1:nrow(node)){Best_strata[[i]] = Strata[[node[i,1]]][[node[i,2]]]}
  
  return(list('Best_strata'=Best_strata,'obj_value'=obj, 'Full_Strata' = Strata[1:max(node_save[[length(node_save)]][,1])], 'node' = node,'obj_vali_save'=obj_vali_save,
              'Best_strata_num'=Best_strata_num,'n_optimal' = n_optimal_save[[Best_strata_num]], 'p_strata' = Tree_p_save[[Best_strata_num]],'sd_strata' = Tree_p_save[[Best_strata_num]],'Validation' = Validation))
}


## Draw new samples from stratification
Draw_Sample_optimal_allocation = function(Best_strata,n_optimal,l){
  Tree = Best_strata$Best_strata; #n_optimal = Best_strata$n_optimal;
  I = length(n_optimal);  #node = Best_strata$node;
  n_allocation = sum(n_optimal); num_sample = rep(0,I)
  Sample = vector(mode='list', length=I)
  draw_num = 1
  while(sum(num_sample == n_optimal)<I){
    if(draw_num >= 10000){
      for(i in 1:I){
        if(num_sample[i] < n_optimal[i]){
          temp = Tree[[i]]
          split_num = dim(temp)[2] - which(colnames(temp) == 'Z') #node[i,1]-1
          variable_select = c(); variable_range = matrix(nrow=split_num,ncol = 2)
          for(ii in 1:split_num){
            variable_select[ii] = as.numeric(substring(colnames(temp)[dim(temp)[2]-split_num+ii],9,10))
            variable_range[ii,] = find_range(temp,split_num,ii)
          }
          
          ## for making samples which are very rare
          X = matrix(runif(5*n_allocation * p,min_x,max_x),ncol=p)
          
          for(ii in 1:split_num){
            X = X[!is.na(cut(X[,variable_select[ii]],breaks = variable_range[ii,])),]
            X = matrix(X,ncol = p)
          }
          Sample[[i]] = matrix(rbind(Sample[[i]],X),ncol = p)
          num_sample[i] = nrow(Sample[[i]])
          if(num_sample[i] > n_optimal[i]){
            Sample[[i]] = Sample[[i]][-((num_sample[i]-(num_sample[i] - n_optimal[i]-1)):num_sample[i]),]
            Sample[[i]] = matrix(Sample[[i]],ncol=p)
            num_sample[i] = nrow(Sample[[i]])
          }
        }
      }
    }else{
      for(i in 1:I){
        if(num_sample[i] < n_optimal[i]){
          temp = Tree[[i]]
          split_num = dim(temp)[2] - which(colnames(temp) == 'Z') #node[i,1]-1
          variable_select = c(); variable_range = matrix(nrow=split_num,ncol = 2)
          for(ii in 1:split_num){
            variable_select[ii] = as.numeric(substring(colnames(temp)[dim(temp)[2]-split_num+ii],9,10))
            variable_range[ii,] = find_range(temp,split_num,ii)
          }
          
          if(x_dist(0) ==  x_dist(0.1)){
            X = matrix(runif(5*n_allocation * p,min_x,max_x),ncol=p)
          }else if((case == 'Gamma_1')|(case == 'Gamma_2')|(case == 'Gamma_3')){
            X = matrix(rgamma(5*n_allocation * p,5,5),ncol=p)
          }else if((case == 't_1')|(case == 't_2')|(case == 't_3')){ # t
            X = matrix(rt(5*n_allocation * p,4),ncol=p)
          }else if((case == 't_trun_1')|(case == 't_trun_2')|(case == 't_trun_3')){ # trun_t
            X = matrix(rtt(5*n_allocation * p,df=1, left=min_x, right=max_x),ncol=p)
          }else{
            X = rmvnorm(n_allocation*5, mean = rep(0, p), sigma = diag(p)) # make enough samples
          }
          
          for(ii in 1:split_num){
            X = X[!is.na(cut(X[,variable_select[ii]],breaks = variable_range[ii,])),]
            X = matrix(X,ncol = p)
          }
          Sample[[i]] = matrix(rbind(Sample[[i]],X),ncol = p)
          num_sample[i] = nrow(Sample[[i]])
          if(num_sample[i] > n_optimal[i]){
            Sample[[i]] = Sample[[i]][-((num_sample[i]-(num_sample[i] - n_optimal[i]-1)):num_sample[i]),]
            Sample[[i]] = matrix(Sample[[i]],ncol=p)
            num_sample[i] = nrow(Sample[[i]])
          }
        }
      }
    }
    draw_num = draw_num + 1
  }
  
  for(i in 1:I){
    Sample[[i]] = data.frame(Sample[[i]])
    colnames(Sample[[i]]) = c(paste('X.',1:ncol(Sample[[i]]),sep=''))
    Sample[[i]]$Y = Make_Y(Sample[[i]]);
    Sample[[i]]$Z = Make_Z(Sample[[i]]$Y,l)
  }      
  return(Sample)
}


################################
##### Draw Tree structure ####
################################

## printing the split history
print_history <- function(branch){
  if(is.null(branch)){
    return(NULL)
  }else{
    start = which(colnames(branch) == 'Z')+1
    range = branch[1,][start:dim(branch)[2]]
    return(range)
  }
}

## collect split history same depth
collect_history <- function(Tree_depth){
  I = length(Tree_depth); collect = vector(mode = "list", length = I);
  for(i in 1:I){
    collect[[i]] = print_history(Tree_depth[[i]])
  }
  return(collect)
}

## find parent node number
find_parent <- function(branch){
  history = print_history(branch)
}


Tree_plot = function(Best_strata){
  ### adjust size 
  size = 8; font = 0.8; symbol = 16 #15; 
  background = 'grey'
  
  Tree = Best_strata$Full_Strata
  depth = max(Best_strata$node[,1])
  coordinate =  vector(mode = "list", length = depth)
  for(d in 0:(depth-2)){
    if(d==0){
      node = Best_strata$node;    parent_node_matrix = node
    }else{
      #node = cbind(rep(depth-d,length(Tree[[depth-d]])), 1:length(Tree[[depth-d]]))
      for(i in 1:nrow(draw_next)){
        sameparent = which(node_total[,3] == draw_next[i,][1] & node_total[,4] == draw_next[i,][2])
        node[sameparent,] = matrix(rep(draw_next[i,],length(sameparent)), nrow = length(sameparent),byrow=T)
      }
      node = unique(node)
      parent_node_matrix = node
    }
    
    ## save coordinate 
    if(d == 0){
      for(i in 1:nrow(node)){
        coordinate[[node[i,1]]][[node[i,2]]] = data.frame(x = i, y = node[i,1])
      }
    }
    
    
    ##find parent node number
    for(i in 1:nrow(node)){
      branch = Tree[[node[i,1]]][[node[i,2]]]
      split_num = length(branch) - which(colnames(branch) == 'Z') ## the number of split
      if(node[i,1]==2){
        parent_node_matrix[i,] = c(1,1)
      }
      else if(node[i,1] > 2){
        parent = collect_history(Tree[[split_num]])
        parent_history = print_history(branch)[-length(print_history(branch))]
        for(k in 1:length(parent)){
          if(length(parent[[k]]) == length(parent_history)){
            if(sum(parent[[k]] == parent_history)==length(parent_history)){
              parent_node = k
            }
          }
        }
        parent_node_matrix[i,] = c(split_num,parent_node)
      }
    }
    #draw depth line
    if(d == 0){
      
      terminal_size = Best_strata$n_optimal / sum(Best_strata$n_optimal) * 15 + 4
      plot(1:nrow(node),-node[,1],pch = symbol ,col='indianred3', cex = terminal_size, ylim = c(-depth-0.5,-0.5), xlim = c(0,nrow(node)+0.5),
           xlab = 'node',ylab = 'depth')
      for(i in 1:nrow(node)){
        temp = Tree[[node[i,1]]][[node[i,2]]]
        variable = colnames(temp)[dim(temp)[2]]
        v = substring(variable,9,11)
        interval = as.character(temp[,dim(temp)[2]][1])
        n = dim(temp)[1]
        text(i,-node[i,1],paste('X',v,'\n',interval,'\n',round(Best_strata$n_optimal[i])),cex=font)
      }
    }
    node_total = cbind(node,parent_node_matrix)
    draw_next = unique(parent_node_matrix[parent_node_matrix[,1]==depth-d-1,])
    for(i in 1:nrow(draw_next)){
      sameparent = which(node_total[,3] == draw_next[i,][1] & node_total[,4] == draw_next[i,][2])
      x_coordi = (coordinate[[depth-d]][[node[sameparent[1],2]]]$x + coordinate[[depth-d]][[node[sameparent[2],2]]]$x) / 2
      coordinate[[draw_next[i,1]]][[draw_next[i,2]]] = data.frame(x = x_coordi,y=depth-d-1)
      
      
      segments(coordinate[[draw_next[i,1]]][[draw_next[i,2]]]$x,-coordinate[[draw_next[i,1]]][[draw_next[i,2]]]$y
               ,coordinate[[depth-d]][[node[sameparent[1],2]]]$x, -coordinate[[depth-d]][[node[sameparent[1],2]]]$y,col='black')
      segments(coordinate[[draw_next[i,1]]][[draw_next[i,2]]]$x,-coordinate[[draw_next[i,1]]][[draw_next[i,2]]]$y
               ,coordinate[[depth-d]][[node[sameparent[2],2]]]$x, -coordinate[[depth-d]][[node[sameparent[2],2]]]$y,col='black')
      
      points(coordinate[[draw_next[i,1]]][[draw_next[i,2]]]$x,-coordinate[[draw_next[i,1]]][[draw_next[i,2]]]$y ,pch = symbol, cex = size,col=background)
      
      
      if(draw_next[1,1] == 1 & draw_next[1,2] == 1){
        temp = Tree[[1]]
        n = dim(temp)[1]
        text(coordinate[[draw_next[i,1]]][[draw_next[i,2]]]$x,-coordinate[[draw_next[i,1]]][[draw_next[i,2]]]$y,paste('Root'),cex=font)
      }else{
        temp = Tree[[draw_next[i,1]]][[draw_next[i,2]]]
        variable = colnames(temp)[dim(temp)[2]]
        v = substring(variable,9,11)
        interval = as.character(temp[,dim(temp)[2]][1])
        n = dim(temp)[1]
        text(coordinate[[draw_next[i,1]]][[draw_next[i,2]]]$x,-coordinate[[draw_next[i,1]]][[draw_next[i,2]]]$y, paste('X',v,'\n',interval),cex=font)
      }
    }
  }
}

########################################################################
## 3-1. Run simulation(MSS)
########################################################################
Simulation_MSS = function(Pilot_sample,Best_strata,l,n_simulation,nj_min,m,I,l_RR,r, Num_iter){
  Tree_save = vector(mode = "list", length = Num_iter)
  n_optimal = round(Best_strata$n_optimal)
  Tree = Best_strata$Best_strata
  p_MSS = Best_strata$p_strata
  Sample_merge = Pilot_sample
  Sample_merge$Z = Make_Z(Pilot_sample$Y,l)
  
  ########################################################################
  #Make estimator
  ########################################################################
  Ehat_MSS_inner = c(); Ehat_MSS_inner_obj = c()
  set.seed(12345 + r) # This r is experiment number that is from main file
  Sample = Draw_Sample_optimal_allocation(Best_strata,n_optimal,l)
  Ehat_MSS_stratum=c(); Ehat_MSS_stratum_obj=c()
  for(i in 1:length(n_optimal)){Ehat_MSS_stratum[i] = p_MSS[i] * mean(Sample[[i]]$Z)}
  for(i in 1:length(n_optimal)){Ehat_MSS_stratum_obj[i] = p_MSS[i]^2 * sd(Sample[[i]]$Z)^2 / n_optimal[i]}
  Ehat_MSS_inner[1] = sum(Ehat_MSS_stratum); Ehat_MSS_inner_obj[1] = sum(Ehat_MSS_stratum_obj)
  #print('iteration: 1')
  #print(paste0('Iteration 1 completed'))
  #print(round(rbind(n_optimal,Ehat_MSS_stratum),4))
  #print(paste('Ehat:',sum(Ehat_MSS_stratum)))
  #print("--------------------------------------------------")
  if(Num_iter == 1){
    Sample_frame = data.frame()
    for(i in 1:length(n_optimal)){Sample_frame = rbind(Sample_frame, Sample[[i]])}
    ### Merge data for sequential sampling
    Sample_merge = rbind(Sample_merge, Sample_frame)
    Tree_save[[1]] = Best_strata
    Ehat_MSS = mean(Ehat_MSS_inner); Ehat_MSS_obj = mean(Ehat_MSS_inner_obj)
    return(list('Best_strata' = Best_strata,'Ehat_MSS' = Ehat_MSS,'Tree'=Tree,'Tree_save'=Tree_save))
  }else{
    Sample_frame = data.frame()
    for(i in 1:length(n_optimal)){Sample_frame = rbind(Sample_frame, Sample[[i]])}
    ### Merge data for sequential sampling
    Sample_merge = rbind(Sample_merge, Sample_frame)
    h_hat <<- find_bandwidth_iteration(Sample_merge, h_hat) # save bandwidths as a global object
    #print(paste0('h_hat:',h_hat))
    
    ## Finding partition again
    Best_strata = Tree_stratification(Sample_merge,n_simulation,nj_min,l,I,m,l_RR)
    Tree_save[[1]] = Best_strata
    #Tree_plot(Best_strata)
    n_optimal = round(Best_strata$n_optimal)
    Tree = Best_strata$Best_strata
    p_MSS = Best_strata$p_strata 
    
    for(mm in 2:Num_iter){
      set.seed(12345 + mm*1000 + r)
      Sample = Draw_Sample_optimal_allocation(Best_strata,n_optimal,l)
      Ehat_MSS_stratum=c()
      for(i in 1:length(n_optimal)){Ehat_MSS_stratum[i] = p_MSS[i] * mean(Sample[[i]]$Z)}
      for(i in 1:length(n_optimal)){Ehat_MSS_stratum_obj[i] = p_MSS[i]^2 * sd(Sample[[i]]$Z)^2 / n_optimal[i]}
      Ehat_MSS_inner[mm] = sum(Ehat_MSS_stratum); Ehat_MSS_inner_obj[mm] = sum(Ehat_MSS_stratum_obj)
      #print(paste('iteration:',mm))
      #print(paste0('Iteration ', mm ,' completed'))
      #print(round(rbind(n_optimal,Ehat_MSS_stratum),4))
      #print(paste('Ehat:',sum(Ehat_MSS_stratum)))
      #print("--------------------------------------------------")
      ### Make Sample into one data.frame
      if(mm != Num_iter){
        Sample_frame = data.frame()
        for(i in 1:length(n_optimal)){Sample_frame = rbind(Sample_frame, Sample[[i]])}
        ### Merge data for sequential sampling
        Sample_merge = rbind(Sample_merge, Sample_frame)
        h_hat <<- find_bandwidth_iteration(Sample_merge, h_hat)
        #print(paste0('h_hat:',h_hat))
        ## Finding partition again
        Best_strata = Tree_stratification(Sample_merge,n_simulation,nj_min,l,I,m,l_RR)
        Tree_save[[mm]] = Best_strata
        #Tree_plot(Best_strata)
        n_optimal = round(Best_strata$n_optimal)
        Tree = Best_strata$Best_strata
        p_MSS = Best_strata$p_strata
      }else{
        Sample_frame = data.frame()
        for(i in 1:length(n_optimal)){Sample_frame = rbind(Sample_frame, Sample[[i]])}
        ### Merge data for sequential sampling
        Sample_merge = rbind(Sample_merge, Sample_frame)
        Tree_save[[mm]] = Sample_merge # To save samples
      }
    }
    Ehat_MSS = mean(Ehat_MSS_inner); Ehat_MSS_obj = mean(Ehat_MSS_inner_obj)
    return(list('Best_strata' = Best_strata,'Ehat_MSS' = Ehat_MSS,'Tree'=Tree,'Tree_save'=Tree_save,'Ehat_MSS_obj' = Ehat_MSS_obj, 'Ehat_MSS_inner' = Ehat_MSS_inner)) 
  }
}


########################################################################
## 3-2. Run simulation(Equi-distance)
########################################################################


Split_data_Equi = function(Pilot_sample,Z,I,p,v,min_x, max_x){
  a_equi = seq(min_x,max_x,length.out=I+1)
  Pilot_sample = cbind(Pilot_sample,cut(Pilot_sample[,v],a_equi,include.lowest = T))
  colnames(Pilot_sample)[dim(Pilot_sample)[2]] = paste('split_X.',v,sep='')
  return(split(Pilot_sample,Pilot_sample[,dim(Pilot_sample)[2]]))
}


Build_tree_Equi = function(Pilot_sample,I,p,min_x, max_x){
  a_equi = seq(min_x,max_x,length.out=I+1)
  tree = vector(mode = "list", length = p+1)
  tree[[1]] = Pilot_sample
  tree[[2]] = append(tree[[2]],Split_data_Equi(tree[[1]],tree[[1]]$Z,I,p,1,min_x, max_x))
  if(p >= 2){
    for(s in 3:(p+1)){
      for(ss in 1:length(tree[[s-1]])){
        tree[[s]] = append(tree[[s]], Split_data_Equi(tree[[s-1]][[ss]],tree[[s-1]][[ss]]$Z,I,p,s-1,min_x, max_x))
      }
    }
  }
  return(tree)
}

Draw_Sample_optimal_allocation_equi = function(Tree,n_optimal,l){
  depth = length(Tree); I = length(n_optimal);  
  n_allocation = sum(n_optimal); num_sample = rep(0,I)
  Sample = vector(mode='list', length=I)
  while(sum(num_sample == n_optimal)<I){
    for(i in 1:I){
      if(num_sample[i] < n_optimal[i]){
        temp = Tree[[depth]][[i]]
        split_num = dim(temp)[2]-(p+2)
        variable_select = c(); variable_range = matrix(nrow=split_num,ncol = 2)
        for(ii in 1:split_num){
          variable_select[ii] = as.numeric(substring(colnames(temp)[dim(temp)[2]-split_num+ii],9,9))
          variable_range[ii,] = find_range(temp,split_num,ii)
        }
        
        if(x_dist(0) ==  x_dist(0.1)){
          X <- matrix(runif(5*n_allocation * p,min_x,max_x),ncol=p)
        }else if(case == "4"){
          X = c()
          for(i in 1:(5 * n_allocation)){
          mu = c(0, 0, 0)
          Sigma = diag(3)
          x_i = rmvnorm(1, mu, Sigma)
          X1 = x_i[, 1]
          X2 = rnorm(1, X1, 1)
          X3 = rnorm(1, X1, 1)
          X4 = x_i[, 2]
          X5 = x_i[, 3]
          X_i = cbind(X1, X2, X3, X4, X5)
          X = rbind(X, X_i)
          }
        }
        else if((case == 'Gamma_1')|(case == 'Gamma_2')|(case == 'Gamma_3')){ # Gamma
          X <- matrix(rgamma(5*n_allocation * p,5,5),ncol=p)
        }else if((case == 't_1')|(case == 't_2')|(case == 't_3')){ # t
          X = matrix(rt(5*n_allocation * p,4),ncol=p)
        }else if((case == 't_trun_1')|(case == 't_trun_2')|(case == 't_trun_3')){ # trun_t
          X = matrix(rtt(5*n_allocation * p,df=1, left=min_x, right=max_x),ncol=p)
        }else{
          X = rmvnorm(n_allocation*5, mean = rep(0, p), sigma = diag(p)) # make enough samples
        }
        
        for(ii in 1:split_num){
          X = X[!is.na(cut(X[,variable_select[ii]],breaks = variable_range[ii,])),]    
          X = matrix(X,ncol = p)
        }
        Sample[[i]] = matrix(rbind(Sample[[i]],X),ncol = p)
        num_sample[i] = nrow(Sample[[i]])
        if(num_sample[i] > n_optimal[i]){
          Sample[[i]] = Sample[[i]][-((num_sample[i]-(num_sample[i] - n_optimal[i]-1)):num_sample[i]),]
          Sample[[i]] = matrix(Sample[[i]],ncol=p)
          num_sample[i] = nrow(Sample[[i]])
        }
      }
    }
  }
  for(i in 1:I){
    Sample[[i]] = data.frame(Sample[[i]])
    colnames(Sample[[i]]) = c(paste('X.',1:ncol(Sample[[i]]),sep=''))
    Sample[[i]]$Y = Make_Y(Sample[[i]]);
    Sample[[i]]$Z = Make_Z(Sample[[i]]$Y,l)
  }      
  return(Sample)
}

Simulation_MEqui = function(Pilot_sample,I,p,l,n_simulation,min_x, max_x,nj_min){
  a_equi = seq(min_x,max_x,length.out=I+1)
  mu_equi = c(); sd_equi = c();s_equi = c(); p_equi = c()
  depth = p+1
  Tree = Build_tree_Equi(Pilot_sample,I,p,min_x, max_x)
  range = vector(mode = "list", length = length(Tree[[depth]]))
  nj_min = round(nj_min)
  for(i in 1:length(Tree[[depth]])){
    temp = Tree[[depth]][[i]]
    p_temp = rep(1,p)
    temp_Pilot = Tree[[1]][paste('X.',substring(colnames(temp)[dim(temp)[2]],9,9),sep='')][,1] #select variable
    n_num <- nrow(temp)
    mu_equi[i] = mean(temp$Z)
    s_equi[i] = mean((temp$Z)^2)
    sd_equi[i] = sqrt(n_num / (n_num-1) * (s_equi[i] - mu_equi[i]^2))
    if(mu_equi[i] == 0 | mu_equi[i] == 1){sd_equi[i] = 0}
    range[[i]] = data.frame(lower = rep(min_x,p),upper = rep(max_x,p))
    ## finding the joint probability
    split_num <- dim(temp)[2]-dim(Tree[[1]])[2]
    if(split_num >= 1){
      for(ii in 1:split_num){  #number of split
        variable = colnames(temp)[dim(temp)[2]-split_num+ii]
        variable = as.numeric(substring(variable,unlist(gregexpr('X',variable))+2,unlist(gregexpr('X',variable))+2))
        interval = as.character(temp[, dim(temp)[2]-split_num+ii][1])
        lower = as.numeric(substring(interval,2,unlist(gregexpr(',',interval))-1))
        upper = as.numeric(substring(interval,unlist(gregexpr(',',interval))+1,unlist(gregexpr(']',interval))-1))
        range[[i]][variable,1] = lower; range[[i]][variable,2] = upper
        p_temp[variable] = integrate(x_dist,lower,upper)$value
      }
    }
    p_equi[i] = prod(p_temp)
  }
  
  ## Draw optimal allocation sample
  n_optimal_equi = round(n_simulation * (p_equi * sd_equi) / (sum(p_equi * sd_equi))) # set the optimal allocation
  #To assign minimum samples to calculate estimated mean and variance
  n_optimal_equi[n_optimal_equi < nj_min] = nj_min #replace optimal size less than minimal allocation size with minimal allocation size
  mu_equi = c(); Ehat_equi =c()
  
  
  ## Make Sample with optimal allocation
  Sample_equi = Draw_Sample_optimal_allocation_equi(Tree,n_optimal_equi,l)
  ########################################################################
  #Make estimator
  ########################################################################
  mu_equi_stratum = c(); Ehat_equi_stratum=c()
  for(i in 1:length(Sample_equi)){
    Ehat_equi_stratum[i] = p_equi[i] * mean(Sample_equi[[i]]$Z)
  }
  Ehat_equi = sum(Ehat_equi_stratum)
  return(Ehat_equi)
}


########################################################################
## 3-3. Run simulation(MC)
########################################################################
Simulation_MC = function(n_simulation,l,R,min_x,max_x){
  Ehat_MC = c()
  for(r in 1:R){
    if(x_dist(0) ==  x_dist(0.1)){
      Sample_MC <- matrix(runif(n_simulation * p,min_x,max_x),ncol=p)
    }else if(x_dist(0) == 0){
      Sample_MC <- matrix(rgamma(n_simulation * p,5,5),ncol=p)
    }else if((case == 't_1')|(case == 't_2')|(case == 't_3')){ # t
      Sample_MC <- matrix(rt(5*n_allocation * p,4),ncol=p)
    }else if((case == 't_trun_1')|(case == 't_trun_2')|(case == 't_trun_3')){ # trun_t
      Sample_MC = matrix(rtt(5*n_allocation * p,df=1, left=min_x, right=max_x),ncol=p)
    }else{
      Sample_MC = rmvnorm(n_simulation, mean = rep(0, p), sigma = diag(p)) # make enough samples
    }
    Ehat_MC[r] = mean(Make_Z(Make_Y(Sample_MC),l))
  }
  return(Ehat_MC)
}
