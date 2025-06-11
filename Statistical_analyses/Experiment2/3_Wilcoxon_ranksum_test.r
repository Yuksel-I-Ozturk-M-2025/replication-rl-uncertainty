### This script reproduces the results of the Wilcoxon rank-sum tests reported 
### in the first paragraph of the section "4.2. Results and discussion" 
###
### Script written by Dr Adnane Ez-zizi (last modified on 16/07/2022). 

###################
# Preliminary steps
###################

# Defining useful paths (change the path "TOP" as appropriate)
TOP = "./RL_under_uncertainty_package/"
WD = paste0(TOP, "Statistical_analyses/Experiment2")
FIT_DATA_EXP2 = paste0(TOP, "Computational_modelling/Experiment2/FitResults_exp2.mat")

### Load necessary packages
library(R.matlab)
library(car)

# Set the working directory
setwd(WD)

### Load necessary data
fit_exp2 = readMat(FIT_DATA_EXP2)

WTA_SU_rhofree = fit_exp2$FitResults.SU.WTA.rhofree[, c(1:4)]; colnames(WTA_SU_rhofree) = c("SubNo", "rho", "alpha", "T")
SARSA_RU = fit_exp2$FitResults.RU.SARSA[, c(1:3)]; colnames(SARSA_RU) = c("SubNo", "alpha", "T")

###############################
# Running the hypothesis tests
###############################

### Comparing exploration rates (T) ###

# Check normality assumption using Shapiro-Wilk normality test 
shapiro.test(WTA_SU_rhofree[, "T"]) # p = 0.06139
shapiro.test(SARSA_RU[, "T"]) # p = 1.927e-07
# there are a few strong ouliers, so we will use Wilcoxon test

# Running Wilcoxon test
wilcox_test_T = wilcox.test(WTA_SU_rhofree[, "T"], SARSA_RU[, "T"], alternative = "two.sided", exact = FALSE)
wilcox_test_T
# data:  WTA_SU_rhofree[, "T"] and SARSA_RU[, "T"]
# W = 296, p-value = 0.01465
# alternative hypothesis: true location shift is not equal to 0

# Estimate effect size
z = qnorm(wilcox_test_T$p.value/2)
r = z / sqrt(61) # -0.3125142

### Comparing learning rates (alpha) ###

# Check normality assumption using Shapiro-Wilk normality test 
shapiro.test(WTA_SU_rhofree[, "alpha"]) # p = 0.0007804
shapiro.test(SARSA_RU[, "alpha"]) # p = 0.0004234
# we will use Wilcoxon test

# Running Wilcoxon test
wilcox_test_alpha = wilcox.test(WTA_SU_rhofree[, "alpha"], SARSA_RU[, "alpha"], alternative = "two.sided", exact = FALSE)
wilcox_test_alpha
# data:  WTA_SU_rhofree[, "alpha"] and SARSA_RU[, "alpha"]
# W = 450, p-value = 0.8337
# alternative hypothesis: true location shift is not equal to 

# Estimate effect size
z = qnorm(wilcox_test_alpha$p.value/2)
r = z / sqrt(61) # -0.02688897
