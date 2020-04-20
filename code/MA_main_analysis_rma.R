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
#     - On reduced datasets
#
#-----------------------------------------------------------------
source(paste0(URL,'code/rma_models.R'),local=TRUE)
source(paste0(URL,'code/rma_models_reduced_data.R'),local=TRUE)

#-----------------------------------------------------------------
#
# Subgroup analysis
# This section produces three plots with the subgroupa analyses
#-----------------------------------------------------------------
source(paste0(URL,'code/subgroups.R'))

#-----------------------------------------------------------------
#
# Pseudo-Funnel plots
#
#-----------------------------------------------------------------
############################################################
# Funnel between arms
############################################################
xl <- bquote(bold(frac(S[OT]^2,S[OC]^2)))
yl <- 'Standard error'

##-- Colors according to mean effect differences
myFunnel(datos1,comparison='Between arms',model=rma.unadj,
         subgroup='effect', zoom=FALSE,xlab=xl,ylab=yl)

##-- Colors according to variance differences
myFunnel(datos1,comparison='Between arms',model=rma.unadj,
         subgroup='variance',zoom=FALSE,xlab=xl,ylab=yl)

############################################################
# Funnel over time
############################################################
xl <- bquote(bold(frac(S[OT]^2,S[BT]^2)))
yl <- 'Standard error'

##-- Colors according to mean effect differences
myFunnel(datos1,comparison='Over time',model=rma.unadj2,
         subgroup='effect',zoom=FALSE,xlab=xl,ylab=yl)

##-- Colors according to mean effect differences
myFunnel(datos1,comparison='Over time',model=rma.unadj2,
         subgroup='variance',zoom=FALSE,xlab=xl,ylab=yl)

############################################################
# Funnel between arms basal
############################################################
xl <- bquote(bold(frac(S[BT]^2,S[BC]^2)))
yl <- 'Standard error'
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
lower.var.expM <- with(datos1,sum(yBetweenArmsRatio < (-qnorm(0.975)*seBetweenArmsRatio)))
greater.var.expM <- with(datos1,sum(yBetweenArmsRatio > (qnorm(0.975)*seBetweenArmsRatio)))
cat('Studies with significative lower variance in the treated arm:',lower.var.expM,'\n')
cat('Studies with significative greater variance in the treated arm:',greater.var.expM,'\n')


##########################################################################
##-- Over time
##########################################################################
##-- Based on model standard deviation
lower.var.outM <- with(datos1,sum(yOverTimeRatioT < (-qnorm(0.975)*seOverTimeRatioT),na.rm=TRUE))
greater.var.outM <- with(datos1,sum(yOverTimeRatioT > (qnorm(0.975)*seOverTimeRatioT),na.rm=TRUE))
cat('Studies with significative lower variance at the end of the study:',lower.var.outM,'\n')
cat('Studies with significative greater variance at the end of the study:',greater.var.outM,'\n')


plot(5,5)