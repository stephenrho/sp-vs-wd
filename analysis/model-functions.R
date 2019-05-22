
# code for estimating k for features and bindings in
# Rhodes et al. Binding in Short-Term Visual Memory: Reassessing Whole Display Interference
# email: srhodes@research.baycrest.org

# this script contains all of the functions used for fitting the k models

calc_rates <- function(k, a, u, N, task, binding=T, informed=T){
  
  if (!(task %in% c("SP", "WD", "DP"))){
    stop("Error in calc_rates: task must be 'SP', 'WD', or 'DP'")
  }
  
  d <- min(1, k/N)
  c <- ifelse(k < (N - 1), 1 - (((N - k)*(N - k - 1))/(N*(N - 1))), 1)
  
  # SINGLE PROBE
  if (task == "SP"){
    if (informed){
      if (binding){
        g <- ((1 - c)*u)/((1 - c)*u + (1 - d)*(1 - u))
      } else {
        g <- u/(u + (1 - d)*(1 - u))
      }
    } else {
      g = u
    }
    
    if (binding){ # binding probe
      h <- a*(c + (1 - c)*g) + (1 - a)*u
      f <- a*(1 - d)*g + (1 - a)*u
    } else{ # feature probe
      h <- a*g + (1 - a)*u
      f <- a*(1 - d)*g + (1 - a)*u
    }
  }
  # WHOLE DISPLAY
  if (task == "WD"){
    if (informed){
      g <- ((1 - c)*u)/((1 - c)*u + (1 - u))
    } else{
      g = u
    }
    # same for feature and binding probes
    h <- a*(c + (1 - c)*g) + (1 - a)*u
    f <- a*g + (1 - a)*u
  }
  # DUAL PROBE
  if (task == "DP"){
    if (informed){
      g <- u/(u + (1 - c)*(1 - u))
    } else{
      g = u
    }
    # same for feature and binding probes
    f <- a*(1 - c)*g + (1 - a)*u
    h <- a*g + (1 - a)*u
  }
  
  return(c(f,h))
}

## models
# 1. same ks for sp and wd
# 2. different ks for sp and wd
# 3. different ks by task and set size
# (each model has separate guessing biases for sp and wd, and different ks for color shape and binding)

check_dat <- function(data, sp=TRUE){
  
  sp_str = ifelse(sp, "SingleProbe", "DualProbe")
  # makes sure its in the right order
  task = rep(c(sp_str, "WholeDisplay"), each=12)
  cond = rep(c("Binding", "Colour", "Shape"), each=4, times=2)
  ss = rep(c(4,6), each=2, times=6)
  tc = rep(c(0,1), times=12)
  
  if (!(all(data$ProbeType == task) & all(data$Condition == cond) & all(data$SetSize == ss) & all(data$TestChange == tc))){
    # make sure things are in right order
    data = data[with(data, order(ProbeType, Condition, SetSize, TestChange)),]
    
    warning(paste0(unique(data$ppt), " data set(s) needed to be reordered"))
  }
  
  return(data)
}

## initial round of models

ll_model1 <- function(par, data, informed, return_pred=F, sp=TRUE){
  # par = c(k_b, k_c, k_s, u_sp, u_wd, a)
  # different ks for features, bindings
  # but same ks across SP and WD
  
  k = c(b = par[1], c = par[2], s = par[3])
  #u = c(SP = par[4], WD = par[5])
  a = par[6]
  if (!sp){
    u = c(DP = par[4], WD = par[5])
    task = rep(c("DP", "WD"), each=6) 
  } else{
    u = c(SP = par[4], WD = par[5])
    task = rep(c("SP", "WD"), each=6) 
  }
  
  upper = c(6, 6, 6, 1, 1, 1)
  lower = rep(0, 6)
  
  # these match the order in data (enforced by check_dat)
  #task = rep(c("SP", "WD"), each=6)
  cond = rep(c("b", "c", "s"), each=2, times=2)
  ss = rep(c(4,6), times=6)
  
  pred = 1:24
  for (i in 1:12){
    pred[((i-1)*2+1):(i*2)] <- calc_rates(k = k[cond[i]], a = a, u = u[task[i]], N = ss[i], task = task[i], binding = cond[i] == 'b', informed = informed)
  }
  
  if (!return_pred){
    ll = dbinom(x = data$RespChange, size = data$Ntrials, prob = pred, log=T)
    ll = ifelse(ll == -Inf, 0, ll)
    
    if (any(par < lower) | any(par > upper)) {
      # penalize for non-sense par values
      ll = ll - 10000
    }
    return(-sum(ll))
  } else{
    return(pred)
  }
}


ll_model2 <- function(par, data, informed, return_pred=F, sp=TRUE){
  # par = c(k_SPb, k_SPc, k_SPs, k_WDb, k_WDc, k_WDs u_sp, u_wd, a)
  # separate ks for features, bindings
  # and separate for SP and WD
  
  if (!sp){
    u = c(DP = par[7], WD = par[8])
    task = rep(c("DP", "WD"), each=6) 
  } else{
    u = c(SP = par[7], WD = par[8])
    task = rep(c("SP", "WD"), each=6) 
  }
  
  #task = rep(c("SP", "WD"), each=6)
  cond = rep(c("b", "c", "s"), each=2, times=2)
  ss = rep(c(4,6), times=6)
  k_lab = paste0(task, cond)
  
  k = par[1:6]
  names(k) = unique(k_lab)
  #u = c(SP = par[7], WD = par[8])
  a = par[9]
  
  upper = c(rep(6, 6), 1, 1, 1)
  lower = rep(0, 9)
  
  pred = 1:24
  for (i in 1:12){
    pred[((i-1)*2+1):(i*2)] <- calc_rates(k = k[k_lab[i]], a = a, u = u[task[i]], N = ss[i], task = task[i], binding = cond[i] == 'b', informed = informed)
  }
  
  if (!return_pred){
    ll = dbinom(x = data$RespChange, size = data$Ntrials, prob = pred, log=T)
    ll = ifelse(ll == -Inf, 0, ll)
    
    if (any(par < lower) | any(par > upper)) {
      # penalize for non-sense par values
      ll = ll - 10000
    }
    return(-sum(ll))
  } else{
    return(pred)
  }
}

ll_model3 <- function(par, data, informed, return_pred=F, sp=TRUE){
  # par = c(k_SPb4, k_SPc4, k_SPs4, k_WDb4, k_WDc4, k_WDs4, 
  #         k_SPb6, k_SPc6, k_SPs6, k_WDb6, k_WDc6, k_WDs6,
  #         u_sp, u_wd, a)
  # freely estimate k by feature, probe type, and set size
  
  k = par[1:12]
  #names(k) = c('SPb4', 'SPc4', 'SPs4', 'WDb4', 'WDc4', 'WDs4', 'SPb6', 'SPc6', 'SPs6', 'WDb6', 'WDc6', 'WDs6')
  #u = c(SP = par[13], WD = par[14])
  a = par[15]
  
  if (!sp){
    u = c(DP = par[13], WD = par[14])
    task = rep(c("DP", "WD"), each=6)
    names(k) = c('DPb4', 'DPc4', 'DPs4', 'WDb4', 'WDc4', 'WDs4', 'DPb6', 'DPc6', 'DPs6', 'WDb6', 'WDc6', 'WDs6')
  } else{
    u = c(SP = par[13], WD = par[14])
    task = rep(c("SP", "WD"), each=6)
    names(k) = c('SPb4', 'SPc4', 'SPs4', 'WDb4', 'WDc4', 'WDs4', 'SPb6', 'SPc6', 'SPs6', 'WDb6', 'WDc6', 'WDs6')
  }
  
  #task = rep(c("SP", "WD"), each=6)
  cond = rep(c("b", "c", "s"), each=2, times=2)
  ss = rep(c(4,6), times=6)
  k_lab = paste0(task, cond, ss)

  upper = c(rep(6, 12), 1, 1, 1)
  lower = rep(0, 15)
  
  pred = 1:24
  for (i in 1:12){
    pred[((i-1)*2+1):(i*2)] <- calc_rates(k = k[k_lab[i]], a = a, u = u[task[i]], N = ss[i], task = task[i], binding = cond[i] == 'b', informed = informed)
  }
  
  if (!return_pred){
    ll = dbinom(x = data$RespChange, size = data$Ntrials, prob = pred, log=T)
    ll = ifelse(ll == -Inf, 0, ll)
    
    if (any(par < lower) | any(par > upper)) {
      # penalize for non-sense par values
      ll = ll - 10000
    }
    return(-sum(ll))
  } else{
    return(pred)
  }
}

## 

fit_model = function(dataset, model_no, informed=T, sp=TRUE, nstart=10){
  
  if (model_no==1){
    model_name = "ll_model1"
    npars = 6
    if (sp){
      par_names = c('k_b', 'k_c', 'k_s', 'u_sp', 'u_wd', 'a')
    } else{
      par_names = c('k_b', 'k_c', 'k_s', 'u_dp', 'u_wd', 'a')
    }
    
    start_upper = c(4,4,4,1,1,1)
    start_lower = rep(0, npars)
  }
  if (model_no==2){
    model_name = "ll_model2"
    npars = 9
    
    if (sp){
      par_names = c('k_SPb', 'k_SPc', 'k_SPs', 'k_WDb', 'k_WDc', 'k_WDs', 'u_sp', 'u_wd', 'a')
    } else{
      par_names = c('k_DPb', 'k_DPc', 'k_DPs', 'k_WDb', 'k_WDc', 'k_WDs', 'u_dp', 'u_wd', 'a')
    }
    start_upper = c(rep(4,6), 1, 1, 1)
    start_lower = rep(0, npars)
  }
  if (model_no==3){
    model_name = "ll_model3"
    npars = 15
    
    if (sp){
      par_names = c('k_SPb4', 'k_SPc4', 'k_SPs4', 'k_WDb4', 'k_WDc4', 'k_WDs4',
                    'k_SPb6', 'k_SPc6', 'k_SPs6', 'k_WDb6', 'k_WDc6', 'k_WDs6',
                    'u_sp', 'u_wd', 'a')
    } else{
      par_names = c('k_DPb4', 'k_DPc4', 'k_DPs4', 'k_WDb4', 'k_WDc4', 'k_WDs4',
                    'k_DPb6', 'k_DPc6', 'k_DPs6', 'k_WDb6', 'k_WDc6', 'k_WDs6',
                    'u_dp', 'u_wd', 'a')
    }
    start_upper = c(rep(4,12), 1, 1, 1)
    start_lower = rep(0, npars)
  }
  
  ids = unique(dataset$ppt)
  
  # create matrix to hold estimates
  res = matrix(NA, nrow = length(ids), ncol = npars+4)
  
  for (s in 1:length(ids)){
    this_id = ids[s]
    
    data = subset(dataset, ppt == this_id)
    data = check_dat(data, sp = sp)
    
    best_pars = rep(NA, npars)
    best_nll = Inf
    best_conv = NA
    for (r in 1:nstart){
      # start 10 times from different starting points
      par = runif(npars, start_lower, start_upper)
      out = optim(par = par, fn = get(model_name), data=data, informed=informed, sp=sp, control=list(maxit=10000))
      
      if(out$value < best_nll){
        best_nll = out$value
        best_pars = out$par
        best_conv = out$convergence
      }
    }
    
    # add best fitting pars to res matrix
    res[s,1:npars] <- best_pars
    res[s,npars+1] <- best_nll
    res[s,npars+2] <- best_conv
    
    n = nrow(data)
    res[s,npars+3] <- 2*best_nll + 2*npars # AIC
    res[s,npars+4] <- 2*best_nll + log(n)*npars # BIC
  }
  
  rownames(res) = ids
  colnames(res) = c(par_names, "nll", "conv", "AIC", "BIC")
  
  return(res)
}


extractFit = function(model_pat){
  # extracts AICs and BICs of particular models (named according to a particular pattern)
  model_names = ls(pattern = model_pat, envir = .GlobalEnv)
  
  aics = matrix(NA, nrow = nrow(get(model_names[1])), ncol = length(model_names))
  
  rownames(aics) = rownames(get(model_names[1]))
  colnames(aics) = model_names
  
  bics = aics
  conv = aics
  
  for (m in model_names){
    aics[,m] <- get(m)[,"AIC"]
    bics[,m] <- get(m)[,"BIC"]
    conv[,m] <- get(m)[,"conv"]
  }
  
  return(list(aic = aics, bic = bics, conv = conv))
}
