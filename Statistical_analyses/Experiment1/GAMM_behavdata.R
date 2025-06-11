### This script reproduces the results reported in the section 
### "2.4.3. Analysis of learning performance using statistical modelling"
###
### Script written by Dr Adnane Ez-zizi (last modified on 16/07/2022). 

################################
# Set WD and load libraries
################################

# Defining the path to the working directory (change the path as appropriate)
WD = "./Statistical_analyses/Experiment1"
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


##########################################################
# Loading data and preprocessing steps before fitting GAMMs
##########################################################

### Load data
data_all = as.data.table(readRDS("./Data/data_all.rds"))

### Number of subjects and trials
N_subj = length(unique(data_all$SubNo)); N_subj # 24
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

### Create a new factor that combine Condition and IsStateIdentical
data_all$CondByIsState = interaction(data_all$Condition,  data_all$IsStateIdentical)

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
toc() # 2.6 sec elapsed

# Model 1.2
tic()
RR1.2 = bam(RewardedOrNot ~ CondByIsState
                          + RewardedOrNot_lag1
                          + s(TrialNo, by = CondByIsState, bs="ad"),
            data = data_all, family = "binomial", 
            discrete = TRUE, nthreads = NTHREADS)
toc() # 0.4 sec elapsed

# Compare the models to decide which variables to keep
AIC(RR1.1, RR1.2)
#             df      AIC
# RR1.1 90.92533 5866.582
# RR1.2 26.26480 6093.004
#
# => We keep the by-participant smooths obviously (it was expected from seeing the individual learning curves)

### Stage 2: first round of removal of fixed effects 
### (whether to remove the non linear interaction of the two factor variables)
# Model 2.1 (full)
RR2.1 = RR1.1

# Model 2.2 (replacing s(TrialNo, by = CondByIsState, bs="ad") with  s(TrialNo, by = Condition, bs="ad"))
tic()
RR2.2 = bam(RewardedOrNot ~ Condition * IsStateIdentical
                         + RewardedOrNot_lag1
                         + s(TrialNo, by = Condition, bs="ad")
                         + s(TrialNo, SubNo, bs="fs"),
             data = data_all, family = "binomial", 
            discrete = TRUE, nthreads = NTHREADS)
toc() # 1.1 sec elapsed

# Model 2.3 (replacing s(TrialNo, by = CondByIsState, bs="ad") with  s(TrialNo, by = IsStateIdentical, bs="ad"))
tic()
RR2.3 = bam(RewardedOrNot ~ Condition * IsStateIdentical
                         + RewardedOrNot_lag1
                         + s(TrialNo, by = IsStateIdentical, bs="ad")
                         + s(TrialNo, SubNo, bs="fs"),
             data = data_all, family = "binomial", 
            discrete = TRUE, nthreads = NTHREADS)
toc() # 1.5 sec elapsed

# Compare the models to decide which variables to keep
AIC(RR2.1, RR2.2, RR2.3)
#             df      AIC
# RR2.1 90.92533 5866.582
# RR2.2 78.51707 5862.849
# RR2.3 98.01343 5903.794
#
# => We replace s(TrialNo, by = CondByIsState, bs="ad") with  s(TrialNo, by = Condition, bs="ad")


### Stage 3: Second round of removal of fixed effects 
# Model 3.1 (Initial)
RR3.1 = RR2.2

# Model 3.2 (replacing s(TrialNo, by = Condition, bs="ad") with s(TrialNo))
tic()
RR3.2 = bam(RewardedOrNot ~ Condition * IsStateIdentical
                         + RewardedOrNot_lag1
                         + s(TrialNo)
                         + s(TrialNo, SubNo, bs="fs"),
             data = data_all, family = "binomial", 
            discrete = TRUE, nthreads = NTHREADS)
toc() # 0.3 sec elapsed
# the warning message "model has repeated 1-d smooths of same variable"
# can be ignored since we need to have smooths over time both in the 
# fixed and random structure

# Model 3.3 (removing the interaction Condition:IsStateIdentical)
tic()
RR3.3 = bam(RewardedOrNot ~ Condition + IsStateIdentical
                         + RewardedOrNot_lag1
                         + s(TrialNo, by = Condition, bs="ad")
                         + s(TrialNo, SubNo, bs="fs"),
            data = data_all, family = "binomial", 
            discrete = TRUE, nthreads = NTHREADS)
toc() # 1.1 sec elapsed

# Compare the models to decide which variables to keep
AIC(RR3.1, RR3.2, RR3.3)
#    df      AIC
# RR3.1 78.51707 5862.849
# RR3.2 99.18056 5906.474
# RR3.3 77.72576 5861.043
#
# => We remove Condition:IsStateIdentical

### Stage 4: Third round of removal of fixed effects 
# Model 4.1 (Initial)
RR4.1 = RR3.3

# Model 4.2 (replacing s(TrialNo, by = Condition, bs="ad") with s(TrialNo))
tic()
RR4.2 = bam(RewardedOrNot ~ Condition + IsStateIdentical
                         + RewardedOrNot_lag1
                         + s(TrialNo)
                         + s(TrialNo, SubNo, bs="fs"),
             data = data_all, family = "binomial", 
            discrete = TRUE, nthreads = NTHREADS)
toc() # 0.3 sec elapsed
# the warning message "model has repeated 1-d smooths of same variable"
# can be ignored since we need to have smooths over time both in the 
# fixed and random structure

# Model 4.3 (removing IsStateIdentical)
tic()
RR4.3 = bam(RewardedOrNot ~ Condition
                         + RewardedOrNot_lag1
                         + s(TrialNo, by = Condition, bs="ad")
                         + s(TrialNo, SubNo, bs="fs"),
             data = data_all, family = "binomial", 
            discrete = TRUE, nthreads = NTHREADS)
toc() # 1.2 sec elapsed

# Compare the models to decide which variables to keep
AIC(RR4.1, RR4.2, RR4.3)
#             df      AIC
# RR4.1 77.72576 5861.043
# RR4.2 98.37830 5904.836
# RR4.3 77.86274 5985.605
#
# => No further reduction of the model

# Final model
tic()
RR_final = bam(RewardedOrNot ~ Condition + IsStateIdentical
                              + RewardedOrNot_lag1
                              + s(TrialNo, by = Condition, bs="ad")
                              + s(TrialNo, SubNo, bs="fs"),
               data = data_all, family = "binomial", nthreads = NTHREADS)
toc() # 1.1 sec elapsed

################################
# Final model summary (Table 1)
################################

set.seed(1)
gam.check(RR_final)
# Basis dimension (k) checking results. Low p-value (k-index<1) may
# indicate that k is too low, especially if edf is close to k'.
# 
#                                 k'    edf k-index p-value
# s(TrialNo):Conditioncontrol  39.00   1.00    1.01    0.89
# s(TrialNo):Conditiontest     39.00   6.03    1.01    0.84
# s(TrialNo,SubNo)            240.00  63.81    1.01    0.87

# Checking autocorrelation between the residuals
acf(residuals(RR_final), lag.max = 10)

### No extreme residuals
which((abs(scale(resid(RR_final)))<2.5) == FALSE) 
# 0

save(RR_final, file="./Results/RR_final.rda", compress='xz')

# Summary of the model (Table 1)
summary(RR_final)
# Formula:
#   RewardedOrNot ~ Condition + IsStateIdentical + RewardedOrNot_lag1 + 
#   s(TrialNo, by = Condition, bs = "ad") + s(TrialNo, 
#   SubNo, bs = "fs")
# 
# Parametric coefficients:
#                     Estimate Std. Error z value Pr(>|z|)
# (Intercept)          0.18686    0.13887   1.346   0.1784
# Conditiontest       -0.40458    0.18224  -2.220   0.0264
# IsStateIdentical1    0.70897    0.06390  11.095   <2e-16
# RewardedOrNot_lag11  0.16090    0.06695   2.403   0.0162
# 
# Approximate significance of smooth terms:
#                                edf  Ref.df  Chi.sq p-value
# s(TrialNo):Conditioncontrol  1.000   1.000   2.643   0.104
# s(TrialNo):Conditiontest     6.032   7.303 101.021  <2e-16
# s(TrialNo,SubNo)            63.814 236.000 271.942  <2e-16
# 
# R-sq.(adj) =  0.137   Deviance explained = 11.5%
# fREML = 6808.3  Scale est. = 1         n = 4776


#############################
# Plot the effects (Figure 6)
#############################

png('./Results/GAMM_behav_effectplots_exp1.png', he=6, wi=11, units='in', res=300)
par(mfrow = c(1,2))
plot(RR_final, select = 2, shade = TRUE, scale = 0, rug = FALSE, cex.lab = 1.2, 
     main = 'Test condition', xlab = 'Trial', ylab = 'Partial effect')
grid()
fig_label(LETTERS[1], cex=1.7) 
set.seed(1)
plot_diff(RR_final, view="TrialNo", 
          comp=list(Condition=c("control", "test")),
          cond = list(RewardedOrNot_lag1 = '1', IsStateIdentical = '1'), # Doesn't matter as there is not interaction with reward
          cex.lab = 1.2,
          main = 'Difference Control - Test',
          xlab = 'Trial',
          ylab = 'Difference in fitted logit-transformed P(rewarded)')
grid()
fig_label(LETTERS[2], cex=1.7) 
dev.off()
# TrialNo window(s) of significant difference(s):
#   102.000000 - 188.000000

#################################################################
# Check significance of the trial effect in the first 100 trials
#####################################################################

# Model run on the first 100 trials from both conditions 
tic()
RR_100tr = bam(RewardedOrNot ~ IsStateIdentical
               + RewardedOrNot_lag1
               + s(TrialNo, bs="ad")
               + s(TrialNo, SubNo, bs="fs"),
               data = data_all[data_all$TrialNo<101,], 
               family = "binomial", nthreads = NTHREADS)
toc() # 9 sec elapsed
# the warning message "model has repeated 1-d smooths of same variable"
# can be ignored since we need to have smooths over time both in the 
# fixed and random structure

summary(RR_100tr)
# Formula:
#   RewardedOrNot ~ IsStateIdentical + RewardedOrNot_lag1 + s(TrialNo, 
#     bs = "ad") + s(TrialNo, SubNo, bs = "fs")
# 
# Parametric coefficients:
#                     Estimate Std. Error z value Pr(>|z|)
# (Intercept)          0.09454    0.13206   0.716   0.4741
# IsStateIdentical1    0.76106    0.09014   8.443   <2e-16
# RewardedOrNot_lag11  0.17577    0.09384   1.873   0.0611
# 
# Approximate significance of smooth terms:
#                    edf Ref.df  Chi.sq p-value
# s(TrialNo)        1.00      1   8.882 0.00288
# s(TrialNo,SubNo) 32.32    238 126.470 < 2e-16
# 
# R-sq.(adj) =  0.0962   Deviance explained = 8.52%
# fREML = 3379.4  Scale est. = 1         n = 2376

