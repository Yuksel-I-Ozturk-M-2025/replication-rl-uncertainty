### This script reproduces the results reported in the section 
### "3.3.3. Analysis of learning performance using statistical modelling"
###
### Script written by Dr Adnane Ez-zizi (last modified on 16/07/2022). 

################################
# Set WD and load libraries
################################

# Defining the path to the working directory (change the path as appropriate)
WD = "./Statistical_analyses/Experiment2"
setwd(WD)

# Loading libraries
library(mgcv) # To build GAMMs
library(itsadug) #  For visualising interactions effects in GAMs
library(data.table)
library(tictoc)
library(ggplot2)

# Number of processor threads to use (change as appropriate)
NTHREADS = 3

options(show.signif.stars=FALSE)

####################
# A helper function
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

################################################################
# Loading data and preprocessing steps before running the GAMMs
################################################################

### Load data
data_all = as.data.table(readRDS("./Data/data_all.rds"))

### Defining the reference level for the uncertainty condition 
### (SU is the ref level)
data_all = within(data_all, UncerCond <- relevel(UncerCond, ref = "SU")) 

### Number of subjects and trials
N_subj = length(unique(data_all$SubNo)); N_subj # 61
N_tr = length(unique(data_all$TrialNo)); N_tr # 200

### Add lags of the dependent variable up to the 5th order
# First add the necessary columns
N_obs = nrow(data_all)
data_all$RewardedOrNot_lag1 = NA
data_all$IsStateIdentical = NA

# Transform RewardedOrNot into a numerical variable
data_all$RewardedOrNot = as.numeric(as.character(data_all$RewardedOrNot))
data_all$IsStateIdentical = as.numeric(as.character(data_all$IsStateIdentical))

for (j in 1:N_subj){
  data_all$RewardedOrNot_lag1[((j-1)*N_tr+2):(j*N_tr)] = data_all$RewardedOrNot[((j-1)*N_tr+1):(j*N_tr-1)]
  data_all$IsStateIdentical[((j-1)*N_tr+2):(j*N_tr)] = (data_all$State[((j-1)*N_tr+2):(j*N_tr)] == data_all$State[((j-1)*N_tr+1):(j*N_tr-1)])
 } 

# Transform into factor columns
data_all$RewardedOrNot = as.factor(as.character(data_all$RewardedOrNot))
data_all$RewardedOrNot_lag1 = as.factor(as.character(data_all$RewardedOrNot_lag1))
data_all$IsStateIdentical = as.factor(as.character(data_all$IsStateIdentical))

### Create a new factor that combine Uncertainty cond and IsStateIdentical
data_all$CondByIsState = interaction(data_all$UncerCond,  data_all$IsStateIdentical)

###################
# Model selection
###################

### Stage 1: Decide whether to include the by-participant smooth
# Model 1.1
tic()
RR1.1 = bam(RewardedOrNot ~ CondByIsState
                          + RewardedOrNot_lag1
                          + s(TrialNo, by = CondByIsState, bs="ad")
                          + s(TrialNo, SubNo, bs="fs"),
            data = data_all, family = "binomial", 
            discrete = TRUE, nthreads = NTHREADS)
toc() # 15.6 sec elapsed

# Model 1.2
tic()
RR1.2 = bam(RewardedOrNot ~ CondByIsState
                          + RewardedOrNot_lag1
                          + s(TrialNo, by = CondByIsState, bs="ad"),
            data = data_all, family = "binomial", 
            discrete = TRUE, nthreads = NTHREADS)
toc() # 1.4 sec elapsed

# Compare the models to decide which variables to keep
AIC(RR1.1, RR1.2)
#              df      AIC
# RR1.1 135.67927 15770.02
# RR1.2  23.78948 15957.59
#
# => We keep the by-participant smooths obviously (it was expected from seeing the individual learning curves)

### Stage 2: first round of removal of fixed effects 
### (whether to remove the non linear interaction of the two factor variables)
# Model 2.1 (full)
RR2.1 = RR1.1

# Model 2.2 (replacing s(TrialNo, by = CondByIsState, bs="ad") with  s(TrialNo, by = UncerCond, bs="ad"))
tic()
RR2.2 = bam(RewardedOrNot ~ UncerCond * IsStateIdentical
                         + RewardedOrNot_lag1
                         + s(TrialNo, by = UncerCond, bs="ad")
                         + s(TrialNo, SubNo, bs="fs"),
             data = data_all, family = "binomial", 
            discrete = TRUE, nthreads = NTHREADS)
toc() # 5.8 sec elapsed

# Model 2.3 (replacing s(TrialNo, by = CondByIsState, bs="ad") with  s(TrialNo, by = IsStateIdentical, bs="ad"))
tic()
RR2.3 = bam(RewardedOrNot ~ UncerCond * IsStateIdentical
                         + RewardedOrNot_lag1
                         + s(TrialNo, by = IsStateIdentical, bs="ad")
                         + s(TrialNo, SubNo, bs="fs"),
             data = data_all, family = "binomial", 
            discrete = TRUE, nthreads = NTHREADS)
toc() # 5.8 sec elapsed

# Compare the models to decide which variables to keep
AIC(RR2.1, RR2.2, RR2.3)
#             df      AIC
# RR2.1 135.6793 15770.02
# RR2.2 125.5357 15773.58
# RR2.3 129.2680 15786.96
#
# => We keep the full model

###############################
# Final model summary (Table 2)
###############################

# Final model (could take several minutes to run)
tic()
RR_final = bam(RewardedOrNot ~ CondByIsState
                             + RewardedOrNot_lag1
                             + s(TrialNo, by = CondByIsState, bs="ad")
                             + s(TrialNo, SubNo, bs="fs"),
               data = data_all, family = "binomial", nthreads = NTHREADS)
toc() # 536 sec elapsed

set.seed(1)
gam.check(RR_final)
# Basis dimension (k) checking results. Low p-value (k-index<1) may
# indicate that k is too low, especially if edf is close to k'.
# 
#                                  k'    edf k-index p-value
# s(TrialNo):CondByIsStateSU.0  39.00   8.59    1.01    0.87
# s(TrialNo):CondByIsStateRU.0  39.00   2.53    1.01    0.88
# s(TrialNo):CondByIsStateSU.1  39.00   1.46    1.01    0.83
# s(TrialNo):CondByIsStateRU.1  39.00   3.85    1.01    0.84
# s(TrialNo,SubNo)             610.00 107.49    1.01    0.91

# Checking autocorrelation between the residuals
acf(residuals(RR_final), lag.max = 10)

save(RR_final, file="./Results/RR_final.rda", compress='xz')
summary(RR_final)
# Formula:
#   RewardedOrNot ~ CondByIsState + RewardedOrNot_lag1 + s(TrialNo, 
#       by = CondByIsState, bs = "ad") + s(TrialNo, SubNo, 
#       bs = "fs")
# 
# Parametric coefficients:
#                      Estimate Std. Error z value Pr(>|z|)
# (Intercept)          0.143653   0.068121   2.109    0.035
# CondByIsStateRU.0    0.023890   0.090266   0.265    0.791
# CondByIsStateSU.1    0.605863   0.055323  10.951  < 2e-16
# CondByIsStateRU.1   -0.007955   0.090273  -0.088    0.930
# RewardedOrNot_lag11  0.228579   0.039149   5.839 5.26e-09
# 
# Approximate significance of smooth terms:
#                                  edf  Ref.df  Chi.sq p-value
# s(TrialNo):CondByIsStateSU.0   8.589  10.339  95.908  <2e-16
# s(TrialNo):CondByIsStateRU.0   2.530   3.131   3.785  0.2744
# s(TrialNo):CondByIsStateSU.1   1.456   1.728   1.278  0.5759
# s(TrialNo):CondByIsStateRU.1   3.853   4.705  13.173  0.0153
# s(TrialNo,SubNo)             107.490 606.000 294.872  <2e-16
# 
# R-sq.(adj) =  0.055   Deviance explained = 5.02%
# fREML =  17251  Scale est. = 1         n = 12139

### model validation: Model without extreme residuals
# You get the same model as we don't have extreme residuals 
which((abs(scale(resid(RR_final)))<2.5) == FALSE) 
# 0

#################################
# plots of the effects (Figure 9)
#################################

png('./Results/GAMM_behav_effectplots_exp2.png', he=9, wi=9, units='in', res=300)
par(mfrow = c(3,2))
plot(RR_final, select = 1, shade = TRUE, scale = 0, rug = FALSE, cex.lab = 1.2, 
     main = 'SU condition: Different state', xlab = '', ylab = 'Partial effect')
grid()
fig_label("A", cex=1.7) 
plot(RR_final, select = 4, shade = TRUE, scale = 0, rug = FALSE, cex.lab = 1.2, 
     main = 'RU condition: Same state', xlab = '', ylab = '')
grid()
fig_label("B", cex=1.7) 

### Difference SU:same - SU:different
plot_diff(RR_final, view="TrialNo", 
          comp=list(CondByIsState=c("SU.1", "SU.0")),
          cond = list(RewardedOrNot_lag1 = '1'), # Doesn't matter as there is not interaction with reward
          rm.ranef = c("SubNo"),
          main = 'Difference SU:same - SU:different',
          xlab = '',
          ylab = 'Difference in fitted logit-transformed P(rewarded)')
grid()
fig_label("C", cex=1.7) 
# TrialNo window(s) of significant difference(s):
#   2.000000 - 40.000000
#   66.000000 - 88.000000
#   100.000000 - 200.000000

### Difference RU:same - RU:different
plot_diff(RR_final, view="TrialNo", 
          comp=list(CondByIsState=c("RU.1", "RU.0")),
          cond = list(RewardedOrNot_lag1 = '1'), # Doesn't matter as there is not interaction with reward
          rm.ranef = c("SubNo"),
          main = 'Difference RU:same - RU:different',
          xlab = '',
          ylab = ''
          #ylab = 'Difference in fitted logit-transformed P(rewarded)'
)
grid()
fig_label("D", cex=1.7) 

plot_diff(RR_final, view="TrialNo", 
          comp=list(CondByIsState=c("SU.1", "RU.1")),
          cond = list(RewardedOrNot_lag1 = '1'), # Doesn't matter as there is not interaction with reward
          rm.ranef = c("SubNo"),
          main = 'Difference SU:same - RU:same',
          #xlab = '',
          xlab = 'Trial',
          ylab = 'Difference in fitted logit-transformed P(rewarded)')
grid()
fig_label("E", cex=1.7) 
# TrialNo window(s) of significant difference(s):
#   4.000000 - 200.000000

### Difference SU:different - RU:different
plot_diff(RR_final, view="TrialNo", 
          comp=list(CondByIsState=c("SU.0", "RU.0")),
          cond = list(RewardedOrNot_lag1 = '1'), # Doesn't matter as there is not interaction with reward
          rm.ranef = c("SubNo"),
          main = 'Difference SU:different - RU:different',
          xlab = 'Trial',
          ylab = ''
          #ylab = 'Difference in fitted logit-transformed P(rewarded)'
)
grid()
fig_label("F", cex=1.7) 
# TrialNo window(s) of significant difference(s):
#   38.000000 - 66.000000
#   94.000000 - 94.000000
#   100.000000 - 140.000000
dev.off()


