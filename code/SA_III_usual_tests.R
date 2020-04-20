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
# Sensitivity analysis III
#
#-----------------------------------------------------------------
##########################################################################
##-- Between arms
##########################################################################
##-- Based on F test
Fest <- with(datos1,(final_sd_T1/final_sd_T2)^2)
LL <- with(datos1,qf(0.025,final_cases_T1,final_cases_T2))
UL <- with(datos1,qf(0.975,final_cases_T1,final_cases_T2))
greater.var.expF <- sum(Fest < LL)
lower.var.expF <- sum(Fest > UL)

##########################################################################
##-- Over time
##########################################################################
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