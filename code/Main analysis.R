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
# Main Analysis --> All Models: 
#     On complete datasets
#     On reduced datasets
#-----------------------------------------------------------------
source(paste0(URL,'code/rma_models.R'))
source(paste0(URL,'code/rma_models_reduced_data.R'))

#-----------------------------------------------------------------
#
# Subgroup analysis --> Figures S2, S3 and S4
#
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
xl <- bquote(bold(frac(V[OT],V[OC])))
yl <- 'Standard error'

##-- Colors according mean effect differences
myFunnel(datos1,comparison='Between arms',model=rma.unadj,
         subgroup='effect', zoom=TRUE,xlab=xl,ylab=yl)

##-- Colors according variance differences
myFunnel(datos1,comparison='Between arms',model=rma.unadj,
         subgroup='variance',zoom=TRUE,xlab=xl,ylab=yl)

############################################################
# Funnel over time
############################################################
xl <- bquote(bold(frac(V[OT],V[BT])))
yl <- 'Standard error'

##-- Colors according mean effect differences
myFunnel(datos1,comparison='Over time',model=rma.unadj2,
         subgroup='effect',zoom=TRUE,xlab=xl,ylab=yl)

############################################################
# Funnel between arms basal
############################################################
xl <- bquote(bold(frac(V[BT],V[BC])))
yl <- 'Standard error'
myFunnel(datos1,comparison='Baseline between arms',model=rma.unadjB,
         subgroup='effect',zoom=FALSE,xlab=xl,ylab=yl)

#-----------------------------------------------------------------
#
# Differences --> Table 1
#
#-----------------------------------------------------------------
##########################################################################
##-- Between arms
##########################################################################
##-- Based on model standard deviation
greater.var.expM <- with(datos1,sum(yBetweenArmsRatio < (-qnorm(0.975)*seBetweenArmsRatio)))
lower.var.expM <- with(datos1,sum(yBetweenArmsRatio > (qnorm(0.975)*seBetweenArmsRatio)))

##-- Based on F test
Fest <- with(datos1,(final_sd_T1/final_sd_T2)^2)
LL <- with(datos1,qf(0.025,final_cases_T1,final_cases_T2))
UL <- with(datos1,qf(0.975,final_cases_T1,final_cases_T2))
greater.var.expF <- sum(Fest < LL)
lower.var.expF <- sum(Fest > UL)

##########################################################################
##-- Over time
##########################################################################
##-- Based on model standard deviation
greater.var.outM <- with(datos1,sum(yOverTimeRatioT < (-qnorm(0.975)*seOverTimeRatioT),na.rm=TRUE))
lower.var.outM <- with(datos1,sum(yOverTimeRatioT > (qnorm(0.975)*seOverTimeRatioT),na.rm=TRUE))

##-- Based on Q test
datos.paired <- datos1[!is.na(datos1$seOverTimeRatioT),]
sdx <- datos.paired$final_sd_T1
sdy <- datos.paired$base_sd_T1
n <- datos.paired$final_cases_T1
rho <- datos.paired$rho
Qest <- var.paired.test(sdx,sdy,n,rho)[[1]]
num <- var.paired.test(sdx,sdy,n,rho)[[2]]
den <- var.paired.test(sdx,sdy,n,rho)[[3]]
LL2 <- qt(0.025,n-2)
UL2 <- qt(0.975,n-2)
greater.var.outQ <- sum(Qest < LL2) 
lower.var.outQ <- sum(Qest > UL2)

##########################################################################
##-- Table 1
##########################################################################
table1 <- matrix(c(lower.var.expF,greater.var.expF,nrow(datos1)-lower.var.expF-greater.var.expF,
                   lower.var.outQ,greater.var.outQ,nrow(datos.paired)-lower.var.outQ-greater.var.outQ),
                 nrow=2,byrow=TRUE)
table1
prop.table(table1,1)
