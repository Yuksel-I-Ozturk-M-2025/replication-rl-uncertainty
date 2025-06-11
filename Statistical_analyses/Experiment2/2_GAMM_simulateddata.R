### This script reproduces the GAMM results reported in the section 
### "4.2. Results and discussion" within the "computational modelling" 
### section
###
### Script written by Dr Adnane Ez-zizi (last modified on 16/07/2022). 

##################################
# Set WD and load libraries
##################################

# Defining useful paths (change the path "TOP" as appropriate)
TOP = "C:/Users/aezzizi/Documents/Work/Lectureship_Suffolk/Research/My_work/RL_under_uncertainty/Package_to_share/"
WD = paste0(TOP, "Statistical_analyses/Experiment2")
# Simulated data from all models
SIM_DATA_EXP2 = paste0(TOP, "Computational_modelling/Experiment2/FitSimulations_exp2.mat")

# Set the working directory 
setwd(WD)

# Loading libraries
library(mgcv) # To build GAMMs
library(itsadug) # For visualising interactions effects in GAMs
library(data.table) 
library(R.matlab)
library(tictoc) 
library(ggplot2)

# Number of processor threads to use (change as appropriate)
NTHREADS = 3

options(show.signif.stars=FALSE)

####################
# Helper functions
####################

### To Add letter labels to figures 
### From: https://logfc.wordpress.com/2017/03/15/adding-figure-labels-a-b-c-in-the-top-left-corner-of-the-plotting-region/
fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {
  
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
  
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
    
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }
  
  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
  
  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100
  
  x1 <- switch(pos,
               topleft     =x[1] + sw, 
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)
  
  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)
  
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
  
  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}

### Make partial effect plots for a given model
partial_effect_plots <- function(gam_model) {
  
  par(mfrow = c(2,2))
  plot(gam_model, select = 1, shade = TRUE, scale = 0, rug = FALSE, cex.lab = 1.2, 
       main = 'SU condition: Different state', xlab = 'Trial', ylab = 'Partial effect')
  grid()
  fig_label("A", cex=1.7) 
  plot(gam_model, select = 3, shade = TRUE, scale = 0, rug = FALSE, cex.lab = 1.2, 
       main = 'SU condition: Same state', xlab = 'Trial', ylab = '')
  grid()
  fig_label("B", cex=1.7) 
  plot(gam_model, select = 2, shade = TRUE, scale = 0, rug = FALSE, cex.lab = 1.2, 
       main = 'RU condition: Different state', xlab = 'Trial', ylab = 'Partial effect')
  grid()
  fig_label("C", cex=1.7) 
  plot(gam_model, select = 4, shade = TRUE, scale = 0, rug = FALSE, cex.lab = 1.2, 
       main = 'RU condition: Same state', xlab = 'Trial', ylab = '')
  grid()
  fig_label("D", cex=1.7) 
  
}

### Make interaction plots for a given model
interaction_plots <- function(gam_model) {
  
  #par(mfrow=c(2,2))
  ### Difference SU:same - SU:different
  plot_diff(gam_model, view="TrialNo", 
            comp=list(CondByIsState=c("SU.1", "SU.0")),
            cond = list(RewardedOrNot_lag1 = '1'), # Doesn't matter as there is not interaction with reward
            rm.ranef = c("SubNo"),
            main = 'Difference SU:same - SU:different',
            xlab = '',
            ylab = 'Difference in fitted logit-transformed P(rewarded)')
  grid()
  fig_label("C", cex=1.7) 
  
  ### Difference RU:same - RU:different
  plot_diff(gam_model, view="TrialNo", 
            comp=list(CondByIsState=c("RU.1", "RU.0")),
            cond = list(RewardedOrNot_lag1 = '1'), # Doesn't matter as there is not interaction with reward
            rm.ranef = c("SubNo"),
            main = 'Difference RU:same - RU:different',
            xlab = '',
            ylab = ''
  )
  grid()
  fig_label("D", cex=1.7) 
  
  ### Difference SU:same - RU:same
  set.seed(1)
  plot_diff(gam_model, view="TrialNo", 
            comp=list(CondByIsState=c("SU.1", "RU.1")),
            cond = list(RewardedOrNot_lag1 = '1'), # Doesn't matter as there is not interaction with reward
            rm.ranef = c("SubNo"),
            main = 'Difference SU:same - RU:same',
            xlab = 'Trial',
            ylab = 'Difference in fitted logit-transformed P(rewarded)')
  grid()
  fig_label("E", cex=1.7) 
  
  ### Difference SU:different - RU:different
  plot_diff(gam_model, view="TrialNo", 
            comp=list(CondByIsState=c("SU.0", "RU.0")),
            cond = list(RewardedOrNot_lag1 = '1'), # Doesn't matter as there is not interaction with reward
            rm.ranef = c("SubNo"),
            main = 'Difference SU:different - RU:different',
            xlab = 'Trial',
            ylab = ''
  )
  grid()
  fig_label("F", cex=1.7) 
  
}

##########################
# Prepare the model data
##########################

### Here we will add simulated data from each model to  to the behavioural dataset

# Load behavioural data
data_all = read.csv(file = "./Data/Data_all.csv")

# Load model data
Rrel_exp2 = readMat(SIM_DATA_EXP2)

Rrel_SU_WTA = Rrel_exp2$R.rel.SU.WTA.rhofree
Rrel_SU_PWRL = Rrel_exp2$R.rel.SU.PWRL.rhofree
Rrel_SU_BISAW = Rrel_exp2$R.rel.SU.BISAW.rhofree
Rrel_RU_BISAW = Rrel_exp2$R.rel.RU.BISAW 
Rrel_RU_SARSA = Rrel_exp2$R.rel.RU.SARSA

# Add proportions of reward collected by each model
data_all$RewardProp_WTA = c(as.vector(Rrel_SU_WTA), as.vector(Rrel_RU_SARSA))
data_all$RewardProp_PWRL = c(as.vector(Rrel_SU_PWRL), as.vector(Rrel_RU_SARSA))
data_all$RewardProp_BISAW = c(as.vector(Rrel_SU_BISAW), as.vector(Rrel_RU_BISAW))

# Add number of reward successes for each model (out of 100 simulations per trial/model)
data_all$RewardSucc_WTA = data_all$RewardProp_WTA * 100
data_all$RewardSucc_PWRL = data_all$RewardProp_PWRL * 100
data_all$RewardSucc_BISAW = data_all$RewardProp_BISAW * 100

# Add number of reward successes for each model (out of 100 simulations per trial/model)
data_all$RewardFail_WTA = 100 - data_all$RewardSucc_WTA
data_all$RewardFail_PWRL = 100 - data_all$RewardSucc_PWRL
data_all$RewardFail_BISAW = 100 - data_all$RewardSucc_BISAW

# Export the final dataset
write.csv(data_all, file = "./Data/Data_all_withModels.csv")
saveRDS(data_all, file = "./Data/Data_all_withModels.rds")

############################################################
# Preprocessing steps before fitting GAMMs
############################################################

### Defining the reference level for the uncertainty condition 
### (SU is the ref level)
data_all$UncerCond = as.factor(as.character(data_all$UncerCond))
data_all = within(data_all, UncerCond <- relevel(UncerCond, ref = "SU")) 

### Number of subjects and trials
N_subj = length(unique(data_all$SubNo)); N_subj # 61
N_tr = length(unique(data_all$TrialNo)); N_tr # 200

### Add lags of the dependent variable up to the 5th order
# First add the necessary columns
N_obs = nrow(data_all)
data_all$RewardProp_PWRL_lag1 = NA
data_all$RewardProp_WTA_lag1 = NA
data_all$RewardProp_BISAW_lag1 = NA
data_all$IsStateIdentical = NA
# Transform IsStateIdentical into a numerical variable
data_all$IsStateIdentical = as.numeric(as.character(data_all$IsStateIdentical))
# Fill in the columns
for (j in 1:N_subj){
  data_all$RewardProp_PWRL_lag1[((j-1)*N_tr+2):(j*N_tr)] = data_all$RewardProp_PWRL[((j-1)*N_tr+1):(j*N_tr-1)]
  data_all$RewardProp_WTA_lag1[((j-1)*N_tr+2):(j*N_tr)] = data_all$RewardProp_WTA[((j-1)*N_tr+1):(j*N_tr-1)]
  data_all$RewardProp_BISAW_lag1[((j-1)*N_tr+2):(j*N_tr)] = data_all$RewardProp_BISAW[((j-1)*N_tr+1):(j*N_tr-1)]
  data_all$IsStateIdentical[((j-1)*N_tr+2):(j*N_tr)] = (data_all$State[((j-1)*N_tr+2):(j*N_tr)] == data_all$State[((j-1)*N_tr+1):(j*N_tr-1)])
 } 

# Transform into factor columns
data_all$IsStateIdentical = as.factor(as.character(data_all$IsStateIdentical))

### Create a new factor that combine Uncertainty cond and IsStateIdentical
data_all$CondByIsState = interaction(data_all$UncerCond,  data_all$IsStateIdentical)

###########################
# Model fit: GAMM for WTA
###########################

# Model summary (p-values reported in Table 5 - "Simple RL" column)
tic()
gamm_WTA = bam(cbind(RewardSucc_WTA, RewardFail_WTA) ~ CondByIsState
                                                     + RewardProp_WTA_lag1
                                                     + s(TrialNo, by = CondByIsState, bs="ad")
                                                     + s(TrialNo, SubNo, bs="fs"),
               data = data_all, family = "binomial", nthreads = NTHREADS)
toc() # 21 sec elapsed
summary(gamm_WTA)
# Formula:
# cbind(RewardSucc_WTA, RewardFail_WTA) ~ CondByIsState + RewardProp_WTA_lag1 + 
#   s(TrialNo, by = CondByIsState, bs = "ad") + s(TrialNo, 
#   SubNo, bs = "fs")
# 
# Parametric coefficients:
#                      Estimate Std. Error  z value Pr(>|z|)
# (Intercept)         -1.381355   0.011924 -115.843   <2e-16
# CondByIsStateRU.0   -0.003093   0.010969   -0.282    0.778
# CondByIsStateSU.1    0.187709   0.005454   34.417   <2e-16
# CondByIsStateRU.1    0.005069   0.010995    0.461    0.645
# RewardProp_WTA_lag1  2.891815   0.016918  170.934   <2e-16
# 
# Approximate significance of smooth terms:
#                                    edf    Ref.df  Chi.sq p-value
# s(TrialNo):CondByIsStateSU.0 1.514e+01 1.769e+01 511.379  <2e-16
# s(TrialNo):CondByIsStateRU.0 1.517e+00 1.858e+00   3.548   0.104
# s(TrialNo):CondByIsStateSU.1 1.367e+01 1.624e+01 203.929  <2e-16
# s(TrialNo):CondByIsStateRU.1 2.302e-04 3.881e-04   0.000   0.991
# s(TrialNo,SubNo)             2.365e+01 2.687e+01 181.942  <2e-16
# 
# Rank: 189/190
# R-sq.(adj) =  0.594   Deviance explained = 58.5%
# fREML =  31288  Scale est. = 1         n = 12139


### Effect plots (Figure 14)
png('./Results/GAMM_WTA_effectplots.png', he=9, wi=9, units='in', res=300)
par(mfrow = c(3,2))
plot(gamm_WTA, select = 1, shade = TRUE, scale = 0, rug = FALSE, cex.lab = 1.2, 
    main = 'SU condition: Different state', xlab = '', ylab = 'Partial effect')
grid()
fig_label("A", cex=1.7) 
plot(gamm_WTA, select = 3, shade = TRUE, scale = 0, rug = FALSE, cex.lab = 1.2, 
    main = 'SU condition: Same state', xlab = '', ylab = '')
grid()
fig_label("B", cex=1.7) 
interaction_plots(gamm_WTA)
dev.off()

###########################
# Model fit: GAMM for PWRL
###########################

# Model summary (p-values reported in Table 5 - "Bayesian RL" column) 
tic()
gamm_PWRL = bam(cbind(RewardSucc_PWRL, RewardFail_PWRL) ~ CondByIsState
                                                     + RewardProp_PWRL_lag1
                                                     + s(TrialNo, by = CondByIsState, bs="ad")
                                                     + s(TrialNo, SubNo, bs="fs"),
               data = data_all, family = "binomial", nthreads = NTHREADS)
toc() # 22 sec elapsed
summary(gamm_PWRL)
# Formula:
# cbind(RewardSucc_PWRL, RewardFail_PWRL) ~ CondByIsState + RewardProp_PWRL_lag1 + 
#   s(TrialNo, by = CondByIsState, bs = "ad") + s(TrialNo, 
#   SubNo, bs = "fs")
# 
# Parametric coefficients:
#                       Estimate Std. Error z value Pr(>|z|)
# (Intercept)          -1.209349   0.012687 -95.323  < 2e-16
# CondByIsStateRU.0    -0.067463   0.011148  -6.051 1.44e-09
# CondByIsStateSU.1     0.096954   0.005413  17.910  < 2e-16
# CondByIsStateRU.1    -0.059613   0.011174  -5.335 9.56e-08
# RewardProp_PWRL_lag1  2.677621   0.018352 145.905  < 2e-16
# 
# Approximate significance of smooth terms:
#                                    edf    Ref.df Chi.sq p-value
# s(TrialNo):CondByIsStateSU.0 1.536e+01 1.798e+01 564.32  <2e-16
# s(TrialNo):CondByIsStateRU.0 1.419e+00 1.718e+00   3.28   0.100
# s(TrialNo):CondByIsStateSU.1 1.227e+01 1.455e+01 341.97  <2e-16
# s(TrialNo):CondByIsStateRU.1 1.733e-04 3.007e-04   0.00   0.992
# s(TrialNo,SubNo)             2.498e+01 2.767e+01 369.13  <2e-16
# 
# Rank: 189/190
# R-sq.(adj) =  0.586   Deviance explained = 57.9%
# fREML =  29431  Scale est. = 1         n = 12139

### Effect plots
png('./Results/GAMM_PWRL_effectplots.png', he=9, wi=9, units='in', res=300)
par(mfrow = c(3,2))
plot(gamm_PWRL, select = 1, shade = TRUE, scale = 0, rug = FALSE, cex.lab = 1.2, 
    main = 'SU condition: Different state', xlab = '', ylab = 'Partial effect')
grid()
fig_label("A", cex=1.7) 
plot(gamm_PWRL, select = 3, shade = TRUE, scale = 0, rug = FALSE, cex.lab = 1.2, 
    main = 'SU condition: Same state', xlab = '', ylab = '')
grid()
fig_label("B", cex=1.7) 
interaction_plots(gamm_PWRL)
dev.off()

############################
# Model fit: GAMM for BISAW
############################

# Model summary (p-values reported in Table 5 - "Sampling" column) 
tic()
gamm_BISAW = bam(cbind(RewardSucc_BISAW, RewardFail_BISAW) ~ CondByIsState
                                                     + RewardProp_BISAW_lag1
                                                     + s(TrialNo, by = CondByIsState, bs="ad")
                                                     + s(TrialNo, SubNo, bs="fs"),
               data = data_all, family = "binomial", nthreads = NTHREADS)
toc() # 22 sec elapsed
summary(gamm_BISAW)
# Formula:
# cbind(RewardSucc_BISAW, RewardFail_BISAW) ~ CondByIsState + RewardProp_BISAW_lag1 + 
#     s(TrialNo, by = CondByIsState, bs = "ad") + s(TrialNo, SubNo, 
#     bs = "fs")

# Parametric coefficients:
#                        Estimate Std. Error z value Pr(>|z|)
# (Intercept)           -1.005614   0.011544 -87.112  < 2e-16
# CondByIsStateRU.0     -0.078781   0.011121  -7.084 1.40e-12
# CondByIsStateSU.1      0.193528   0.005416  35.730  < 2e-16
# CondByIsStateRU.1     -0.065737   0.011147  -5.897 3.69e-09
# RewardProp_BISAW_lag1  2.290118   0.015963 143.464  < 2e-16

# Approximate significance of smooth terms:
#                                 edf Ref.df Chi.sq p-value
# s(TrialNo):CondByIsStateSU.0 23.160 26.322 787.93 < 2e-16
# s(TrialNo):CondByIsStateRU.0  1.002  1.003  26.67 1.3e-06
# s(TrialNo):CondByIsStateSU.1 14.029 16.779 335.69 < 2e-16
# s(TrialNo):CondByIsStateRU.1  1.208  1.388  26.84 6.4e-06
# s(TrialNo,SubNo)             24.424 27.331 407.23 < 2e-16

# Rank: 189/190
# R-sq.(adj) =  0.469   Deviance explained = 46.3%
# fREML =  38693  Scale est. = 1         n = 12139

### Effect plots
png('./Results/GAMM_BISAW_effectplots.png', he=9, wi=9, units='in', res=300)
par(mfrow = c(3,2))
plot(gamm_BISAW, select = 1, shade = TRUE, scale = 0, rug = FALSE, cex.lab = 1.2, 
     main = 'SU condition: Different state', xlab = '', ylab = 'Partial effect')
grid()
fig_label("A", cex=1.7) 
plot(gamm_BISAW, select = 3, shade = TRUE, scale = 0, rug = FALSE, cex.lab = 1.2, 
     main = 'SU condition: Same state', xlab = '', ylab = '')
grid()
fig_label("B", cex=1.7) 
interaction_plots(gamm_BISAW)
dev.off()
