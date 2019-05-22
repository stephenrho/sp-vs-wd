
# code for estimating k for features and bindings in
# Rhodes et al. Binding in Short-Term Visual Memory: Reassessing Whole Display Interference
# email: srhodes@research.baycrest.org

# this script fits the models using the functions sourced below

source("analysis/model-functions.R")

RUN = F # T = run models, F = load models (if already run)

if (RUN){
  
  ### fit models
  # experiment 1
  e1 <- read.csv("data/Exp1.csv")
  
  m1_e1_inf = fit_model(dataset = e1, model_no = 1, informed = T)
  m2_e1_inf = fit_model(dataset = e1, model_no = 2, informed = T)
  m3_e1_inf = fit_model(dataset = e1, model_no = 3, informed = T)
  
  m1_e1 = fit_model(dataset = e1, model_no = 1, informed = F)
  m2_e1 = fit_model(dataset = e1, model_no = 2, informed = F)
  m3_e1 = fit_model(dataset = e1, model_no = 3, informed = F)
  
  # experiment 2
  e2 <- read.csv("data/Exp2.csv")
  
  m1_e2_inf = fit_model(dataset = e2, model_no = 1, informed = T)
  m2_e2_inf = fit_model(dataset = e2, model_no = 2, informed = T)
  m3_e2_inf = fit_model(dataset = e2, model_no = 3, informed = T)
  
  m1_e2 = fit_model(dataset = e2, model_no = 1, informed = F)
  m2_e2 = fit_model(dataset = e2, model_no = 2, informed = F)
  m3_e2 = fit_model(dataset = e2, model_no = 3, informed = F)
  
  # experiment 3
  e3 <- read.csv("data/Exp3.csv")
  
  m1_e3_inf = fit_model(dataset = e3, model_no = 1, informed = T, sp=FALSE)
  m2_e3_inf = fit_model(dataset = e3, model_no = 2, informed = T, sp=FALSE)
  m3_e3_inf = fit_model(dataset = e3, model_no = 3, informed = T, sp=FALSE)
  
  m1_e3 = fit_model(dataset = e3, model_no = 1, informed = F, sp=FALSE)
  m2_e3 = fit_model(dataset = e3, model_no = 2, informed = F, sp=FALSE)
  m3_e3 = fit_model(dataset = e3, model_no = 3, informed = F, sp=FALSE)
  
  save(list = ls(pattern = "m[0-9]_e"), file = "analysis/k-res.RData")
  
} else {
  load("analysis/k-res.RData")
}


## E1
e1_fit = extractFit(model_pat = "m[0-9]_e1")

# convergence
apply(e1_fit$conv, 2, function(x) sum(x==0))

# by individual
table(apply(e1_fit$aic, 1, function(x) colnames(e1_fit$aic)[which(x==min(x))]))
table(apply(e1_fit$bic, 1, function(x) colnames(e1_fit$bic)[which(x==min(x))]))

# overall
barplot(apply(e1_fit$aic, 2, mean), col=grey.colors(8))
abline(h = min(apply(e1_fit$aic, 2, mean)), col="red", lty=2)

barplot(apply(e1_fit$bic, 2, mean), col=grey.colors(8))
abline(h = min(apply(e1_fit$bic, 2, mean)), col="red", lty=2)


## E2
e2_fit = extractFit(model_pat = "m[0-9]_e2")
# convergence
apply(e2_fit$conv, 2, function(x) sum(x==0))

# by individual
table(apply(e2_fit$aic, 1, function(x) colnames(e2_fit$aic)[which(x==min(x))]))
table(apply(e2_fit$bic, 1, function(x) colnames(e2_fit$bic)[which(x==min(x))]))

# overall
barplot(apply(e2_fit$aic, 2, mean), col=grey.colors(8))
abline(h = min(apply(e2_fit$aic, 2, mean)), col="red", lty=2)

barplot(apply(e2_fit$bic, 2, mean), col=grey.colors(8))
abline(h = min(apply(e2_fit$bic, 2, mean)), col="red", lty=2)


## E3
e3_fit = extractFit(model_pat = "m[0-9]_e3")
# convergence
apply(e3_fit$conv, 2, function(x) sum(x==0))

# by individual
table(apply(e3_fit$aic, 1, function(x) colnames(e3_fit$aic)[which(x==min(x))]))
table(apply(e3_fit$bic, 1, function(x) colnames(e3_fit$bic)[which(x==min(x))]))

# overall
barplot(apply(e3_fit$aic, 2, mean), col=grey.colors(8))
abline(h = min(apply(e3_fit$aic, 2, mean)), col="red", lty=2)

barplot(apply(e3_fit$bic, 2, mean), col=grey.colors(8))
abline(h = min(apply(e3_fit$bic, 2, mean)), col="red", lty=2)


### anovas on estimate of k
# reshape to long

library(tidyr)

reshape_m2_k <- function(par_wide, sp=TRUE){
  
  if (!is.data.frame(par_wide)){
    par_wide = as.data.frame(par_wide)
  }
  
  par_wide$ppt = 1:nrow(par_wide)
  
  if (sp){
    cols = c("k_SPb", "k_SPc", "k_SPs", "k_WDb", "k_WDc", "k_WDs", "ppt")
    sp_str = "SP"
    par_long = gather(par_wide[, cols], "cond", "k", k_SPb:k_WDs, factor_key=TRUE) 
  } else {
    cols = c("k_DPb", "k_DPc", "k_DPs", "k_WDb", "k_WDc", "k_WDs", "ppt")
    sp_str = "DP"
    par_long = gather(par_wide[, cols], "cond", "k", k_DPb:k_WDs, factor_key=TRUE)
  }
  
  # add new columns
  par_long$probe = "WD"
  par_long$probe[grep(pattern = sp_str, par_long$cond)] = sp_str
  
  par_long$feature = "B"
  par_long$feature[grep(pattern = "c", par_long$cond)] = "C"
  par_long$feature[grep(pattern = "s", par_long$cond)] = "S"
  
  par_long = within(par_long, {
    ppt = as.factor(ppt)
    probe = as.factor(probe)
    feature = as.factor(feature)
  })
  
  return(par_long)
}

m2_e1_inf_l = reshape_m2_k(m2_e1_inf)
m2_e2_inf_l = reshape_m2_k(m2_e2_inf)
m2_e3_inf_l = reshape_m2_k(m2_e3_inf, sp = F)


e1_aov <- aov(k ~ probe*feature + Error(ppt/(probe*feature)), data = m2_e1_inf_l, contrasts = contr.sum)
summary(e1_aov)

e2_aov <- aov(k ~ probe*feature + Error(ppt/(probe*feature)), data = m2_e2_inf_l, contrasts = contr.sum)
summary(e2_aov)

e3_aov <- aov(k ~ probe*feature + Error(ppt/(probe*feature)), data = m2_e3_inf_l, contrasts = contr.sum)
summary(e3_aov)

