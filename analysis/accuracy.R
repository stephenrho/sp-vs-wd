
# code for analysis of change detection accuracy in
# Rhodes et al. Binding in Short-Term Visual Memory: Reassessing Whole Display Interference
# email: srhodes@research.baycrest.org

library(lme4)
library(car)

RUN = F # T = run models, F = load models (if already run)

if (RUN){
  ## Experiment 1
  e1 <- read.csv("data/Exp1.csv")
  
  e1 <- within(e1, {
    ppt = as.factor(ppt)
    SetSize = as.factor(SetSize)
  })
  
  e1$Ncorr <- with(e1, ifelse(TestChange == 0, Ntrials - RespChange, RespChange))
  
  contrasts(e1$Condition) <- cbind(FvB = c(-1, 1/2, 1/2), CvS = c(0, 1, -1))
  contrasts(e1$SetSize) <- c(-1, 1)
  contrasts(e1$ProbeType) <- c(-1, 1)
  
  e1_glmer = glmer(cbind(Ncorr, Ntrials - Ncorr) ~ ProbeType*Condition*SetSize + (1 + ProbeType+Condition+SetSize | ppt), data = e1, family = binomial(link = "logit"))
  #control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun=1e+05)))
  
  Anova(e1_glmer, type = 'III')
  
  ## Experiment 2
  e2 <- read.csv("data/Exp2.csv")
  
  e2 <- within(e2, {
    ppt = as.factor(ppt)
    SetSize = as.factor(SetSize)
  })
  
  e2$Ncorr <- with(e2, ifelse(TestChange == 0, Ntrials - RespChange, RespChange))
  
  contrasts(e2$Condition) <- cbind(FvB = c(-1, 1/2, 1/2), CvS = c(0, 1, -1))
  contrasts(e2$SetSize) <- c(-1, 1)
  contrasts(e2$ProbeType) <- c(-1, 1)
  
  e2_glmer = glmer(cbind(Ncorr, Ntrials - Ncorr) ~ ProbeType*Condition*SetSize + (1 + ProbeType+Condition+SetSize | ppt), data = e2, family = binomial(link = "logit"))
  
  Anova(e2_glmer, type = "III")
  
  ## Experiment 3
  e3 <- read.csv("data/Exp3.csv")
  
  e3 <- within(e3, {
    ppt = as.factor(ppt)
    SetSize = as.factor(SetSize)
  })
  
  e3$Ncorr <- with(e3, ifelse(TestChange == 0, Ntrials - RespChange, RespChange))
  
  contrasts(e3$Condition) <- cbind(FvB = c(-1, 1/2, 1/2), CvS = c(0, 1, -1))
  contrasts(e3$SetSize) <- c(-1, 1)
  contrasts(e3$ProbeType) <- c(-1, 1)
  
  e3_glmer = glmer(cbind(Ncorr, Ntrials - Ncorr) ~ ProbeType*Condition*SetSize + (1 + ProbeType+Condition+SetSize | ppt), data = e3, family = binomial(link = "logit"))
  
  Anova(e3_glmer, type = "III")
  
  ## bootstrap ses for figures
  
  newdat = data.frame(ProbeType = rep(levels(e1$ProbeType), each=6), Condition=rep(levels(e1$Condition), each = 2), SetSize=levels(e1$SetSize))
  
  e1_boot <- bootMer(e1_glmer, FUN=function(x) predict(x, newdat, re.form=NA),
                     nsim=250)
  
  e2_boot <- bootMer(e2_glmer, FUN=function(x) predict(x, newdat, re.form=NA),
                     nsim=250)
  
  newdat3 = data.frame(ProbeType = rep(levels(e3$ProbeType), each=6), Condition=rep(levels(e3$Condition), each = 2), SetSize=levels(e3$SetSize))
  
  e3_boot <- bootMer(e3_glmer, FUN=function(x) predict(x, newdat3, re.form=NA),
                     nsim=250)
  
  ## Joint analysis
  
  e1$exp = 1
  e2$exp = 2
  e3$exp = 3
  
  eall = rbind(e1, e2, e3)
  
  eall$id = as.factor(with(eall, paste0(exp, ppt)))
  eall$ProbeType[eall$ProbeType=="DualProbe"] = 'SingleProbe'
  eall = droplevels(eall)
  
  contrasts(eall$Condition) <- cbind(FvB = c(-1, 1/2, 1/2), CvS = c(0, 1, -1))
  contrasts(eall$SetSize) <- c(-1, 1)
  contrasts(eall$ProbeType) <- c(-1, 1)
  
  eall$exp = as.factor(eall$exp)
  contrasts(eall$exp) <- cbind(SvD = c(-1/2, -1/2, 1), OvT = c(-1, 1, 0))
  
  eall_glmer = glmer(cbind(Ncorr, Ntrials - Ncorr) ~ exp*ProbeType*Condition*SetSize + (1 + ProbeType+Condition+SetSize | id), data = eall, family = binomial(link = "logit"))
  
  Anova(eall_glmer, type = "III")
  
  # save
  save.image("analysis/acc-res.RData")
} else {
  load("analysis/acc-res.RData")
}
