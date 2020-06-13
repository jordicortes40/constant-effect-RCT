rm(list=ls())

#-----------------------------------------------------------------
#
# Read data
#
#-----------------------------------------------------------------
URL <- 'https://raw.githubusercontent.com/jordicortes40/constant-effect-RCT/master/'
source(paste0(URL,'code/read_data.R')) 

#-----------------------------------------------------------------
#
# Main Analysis --> All Models. This section produces all model outputs: 
#     - On complete datasets
#
#-----------------------------------------------------------------
##-- See file "../results_tables/MA_models.txt" for all the models
source(paste0(URL,'code/rma_models.R'),local=TRUE)

##-- Confidence intervals for heterogeneity
# Between arms
confint(rma.unadj)  # Heterogeneity between arms in unadjusted model
confint(rma.adj)    # Heterogeneity between arms in adjusted model
confint(rma.adjoff) # Heterogeneity between arms in model with offset (beta=1)

# Over time
confint(rma.unadj2) # Heterogeneity over time in unadjusted model
confint(rma.adj2)   # Heterogeneity over time in adjusted model
confint(rma.adjoff2)# Heterogeneity between arms in model with offset (beta=1)

##-- Test for assymmetries in the funnel plot
ranktest(rma.adj)
ranktest(rma.adj2)

#-----------------------------------------------------------------
#
# Subgroup analysis
# This section produces three plots with the subgroup analyses
#-----------------------------------------------------------------
source(paste0(URL,'code/subgroups.R'))

#-----------------------------------------------------------------
#
# Funnel plots
#
#-----------------------------------------------------------------
############################################################
# Funnel between arms
############################################################
##-- Labels for axes
xl <- bquote(bold(frac(S[OT]^2,S[OC]^2)))
yl <- 'Standard error'

##-- Colors according to mean treatment effect differences (n=208)
myFunnel(datos1,comparison='Between arms',model=rma.unadj,
         subgroup='effect', zoom=FALSE,xlab=xl,ylab=yl)


##-- Colors according to variance of treatment effect differences (n=208)
png('../results_figures/MA_funnel_BA.png',width=1530,height = 1190,res=144)
myFunnel(datos1,comparison='Between arms',model=rma.unadj,
         subgroup='variance',zoom=FALSE,xlab=xl,ylab=yl)
dev.off()

############################################################
# Funnel over time
############################################################
##-- Labels for axes
xl <- bquote(bold(frac(S[OT]^2,S[BT]^2)))
yl <- 'Standard error'

##-- Colors according to mean treatment effect differences (n=95)
myFunnel(datos1,comparison='Over time',model=rma.unadj2,
         subgroup='effect',zoom=FALSE,xlab=xl,ylab=yl)

##-- Colors according to variance of treatment effect differences (n=95)
png('../results_figures/MA_funnel_OT.png',width=1530,height = 1190,res=144)
myFunnel(datos1,comparison='Over time',model=rma.unadj2,
         subgroup='variance',zoom=FALSE,xlab=xl,ylab=yl)
dev.off()

############################################################
# Funnel between arms at baseline
############################################################
##-- Labels for axes
xl <- bquote(bold(frac(S[BT]^2,S[BC]^2)))
yl <- 'Standard error'

##-- Colors according to mean treatment effect differences
myFunnel(datos1,comparison='Baseline between arms',model=rma.unadjB,
         subgroup='effect',zoom=FALSE,xlab=xl,ylab=yl)

#-----------------------------------------------------------------
#
# Number of studies with different variances
#
#-----------------------------------------------------------------
##########################################################################
##-- Between arms
##########################################################################
##-- Based on model standard deviation
lower_var_BA <- with(datos1,yBetweenArmsRatio < (-qnorm(0.975)*seBetweenArmsRatio))
greater_var_BA <- with(datos1,yBetweenArmsRatio > (qnorm(0.975)*seBetweenArmsRatio))
significant_BA <- lower_var_BA | greater_var_BA
cat('Studies with significative lower variance in the treated arm:',sum(lower_var_BA),'\n')
cat('Studies with significative greater variance in the treated arm:',sum(greater_var_BA),'\n')


##########################################################################
##-- Over time
##########################################################################
##-- Based on model standard deviation
lower_var_OT <- with(datos1,yOverTimeRatioT < (-qnorm(0.975)*seOverTimeRatioT),na.rm=TRUE)
greater_var_OT <- with(datos1,yOverTimeRatioT > (qnorm(0.975)*seOverTimeRatioT),na.rm=TRUE)
significant_OT <- lower_var_OT | greater_var_OT
cat('Studies with significative lower variance at the end of the study:',sum(lower_var_OT,na.rm=TRUE),'\n')
cat('Studies with significative greater variance at the end of the study:',sum(greater_var_OT,na.rm=TRUE),'\n')







#-------------------------------------------------------------------
#
# Simulation: Compare the results from the delta method with 
# empirical SE obtained by simulation
#
#-------------------------------------------------------------------
source('SE_validation_BA.R')  # Between arms
source('SE_validation_OT.R')  # Over-time

#-----------------------------------------------------------------
#
# Combining two types of comparison
#
#-----------------------------------------------------------------
##-- Table with information available in comparison over time vs. the significance in comparsion between arms
t_2comparisons_1 <- table(is.na(significant_OT),significant_BA) 
oddsratio(t_2comparisons_1,rev="rows")

##-- Table with significance in comparison over time vs. the significance in comparison between arms
t_2comparisons_2 <- table(significant_OT,significant_BA) 
oddsratio(t_2comparisons_2)
