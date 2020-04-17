#-----------------------------------------------------------------
#
# All scripts are in 'http://www-eio.upc.es/teaching/best/variability_data/'
# You can see them in:
# http://www-eio.upc.es/teaching/best/variability_data/Main.R
# http://www-eio.upc.es/teaching/best/variability_data/functions.R
# http://www-eio.upc.es/teaching/best/variability_data/rma_models.R
# http://www-eio.upc.es/teaching/best/variability_data/rma_models_reduced_data.R
# http://www-eio.upc.es/teaching/best/variability_data/subgroups.R
#
#-----------------------------------------------------------------

##-- Remove objects in memory
rm(list=ls())

#-----------------------------------------------------------------
#
# Install and load packages
#
#-----------------------------------------------------------------
##-- Install packages
list.of.packages <- c('ggplot2','weights','catspec','alabama','metafor')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

##-- Load packages
library(ggplot2)
library(weights)
library(catspec)
library(alabama)
library(metafor)

##-- Penalize scientific notation
options(scipen=3)

#-----------------------------------------------------------------
#
# Load specific functions
#
#-----------------------------------------------------------------
URL <- 'http://www-eio.upc.es/teaching/best/variability_data/'
source(paste0(URL,'functions.R'))

#-----------------------------------------------------------------
#
# Read data
#
#-----------------------------------------------------------------

datos1 <- read.table(url(paste0(URL,'data.csv')),header=TRUE,sep=";",stringsAsFactors = TRUE,quote = "")
closeAllConnections()

#-----------------------------------------------------------------
#
# Create new variables: outcomes and their standard errors
#
#-----------------------------------------------------------------
##-- Between arms --> Baseline log variance ratio
datos1$yBaselineRatio <- with(datos1,2*log(base_sd_T1/base_sd_T2))                               # Outcome
datos1$seBaselineRatio <- with(datos1,sqrt(2*(1/(base_cases_T1-2)+1/(base_cases_T2-2))))         # Standard error (n-2 = df-1 --> best approximation)


##-- Between arms  --> Outcome log variance ratio
datos1$yBetweenArmsRatio <- with(datos1,2*log(final_sd_T1/final_sd_T2))                          # Outcome
datos1$seBetweenArmsRatio <- with(datos1,sqrt(2*(1/(final_cases_T1-2)+1/(final_cases_T2-2))))    # Standard error (n-2 = df-1 --> best approximation)


##-- Over-time --> Experimental log variance ratio
datos1$yOverTimeRatioT <- with(datos1,2*log(final_sd_T1/base_sd_T1))                             # Outcome Treateds

cor.var <- datos1$rho^2                                                                          # Correlation between variances are aprox. the squared correlation between measures (97 available)
cov.log <- with(datos1,log(1 + 2*cor.var/(final_cases_T1-2)))                                    # Covariance between log variances
datos1$seOverTimeRatioT <- with(datos1,sqrt(2*(1/(final_cases_T1-2) +                            # Standard error (adding covariance) for over-time ratio
                                                 1/(final_cases_T1-2) -
                                                 cov.log)))       


##-- Over-time --> Reference log variance ratio
datos1$yOverTimeRatioC <- with(datos1,2*log(final_sd_T2/base_sd_T2))                             # Outcome Controls
datos1$seOverTimeRatioC <- with(datos1,sqrt(2*(1/(final_cases_T2-2) +                            # Standard error in control group
                                                 1/(final_cases_T2-2) -
                                                 cov.log))) 
#-----------------------------------------------------------------
#
# Descriptive
#
#-----------------------------------------------------------------
##########################################################################
##-- Subgroups
##########################################################################
cat('Tables for subgroups----------------------------------------\n')
cat('\nNumber of medical areas-----\n');ctab(factor(datos1$N_Areas_WOS),dec.places=1)
cat('\nIntervention type-----------\n');ctab(datos1$Intervention_type,dec.places=1)
cat('\nOutcome type----------------\n');ctab(datos1$Outcome_type,dec.places=1)
cat('\nCondition type--------------\n');ctab(datos1$Condition_type,dec.places=1)
cat('\nMeasurement type------------\n');ctab(datos1$Measurement_type,dec.places=1)
cat('\nSignificant-----------------\n');ctab(datos1$significant,dec.places = 1)


##########################################################################
##-- Deviations - Table S2
##########################################################################
cat('\nDescriptive Variances with logs----------------------------------------\n')
TableS2 <- with(datos1,matrix(c(fd(2*log(base_sd_T2)),
                               fd(2*log(final_sd_T2)),
                               fd(2*log(base_sd_T1)),
                               fd(2*log(final_sd_T1)),
                               fd(2*log(final_sd_T1)-2*log(final_sd_T2)),
                               fd(2*log(final_sd_T1)-2*log(base_sd_T1)),
                               fd(2*log(final_sd_T1)-2*log(base_sd_T1)-
                                    (2*log(final_sd_T2)-2*log(base_sd_T2)))),nrow=7,byrow=TRUE))
rownames(TableS2) <- paste0('log',paste0('_Var_',c('BasalC','FinalC','BasalT','FinalT','DifT','DifO','DifDif')))
colnames(TableS2) <- c('n','mean','sd','P0%','P25%','P50%','P75%','P100%')
TableS2

##########################################################################
##-- Variances - Descriptive figure between arms
##########################################################################
##-- Previous calculations ------------------------------------------------
ord <- order(datos1$final_sd_T2)
ctrl <- datos1$final_sd_T2[ord]
trt <- datos1$final_sd_T1[ord]
d_figure1 <- data.frame(Treated=trt^2,Control=ctrl^2,Order=order(ctrl),ratio=(trt/ctrl)^2,
                        Group=ifelse(trt>ctrl,'Greater variance in Treated','Greater variance in Control'))


##-- Arrows Graphic ------------------------------------------------
ggplot(d_figure1,aes(x=Order,y=Control,colour=Group)) + geom_point(size=0) + 
  scale_y_log10() +
  geom_segment(aes(x = Order, y = Control, xend = Order, yend = Treated, colour=Group),
               arrow=arrow(length = unit(0.18,"cm")),size=1) +
  theme(legend.position="bottom",
        axis.title = element_text(face='bold'),
        legend.title = element_blank()) +
  xlab('Rank according to Control Outcome variability') + ylab('Variance')

##-- Histogram ------------------------------------------------
ggplot(d_figure1,aes(x=ratio)) + geom_histogram(color='white') + scale_x_log10(limits=c(1/25,25)) +
  ylab('n') + xlab(expression(bold(S[OT]^2/S[OC]^2))) +
  geom_vline(xintercept = 1,linetype=2,color='white',size=1.3) +
  annotate("text",x = 0.1, y = 40, label = "paste(bold(\"Higher \"),bold(S[OC]^2))", parse = TRUE, size=8)+
  annotate("text",x = 10, y = 40, label = "paste(bold(\"Higher \"),bold(S[OT]^2))", parse = TRUE, size=8)+
  theme(axis.title = element_text(face='bold',size=13),
        axis.text = element_text(face='bold',size=13)) 

##########################################################################
##-- Variances - Descriptive figures over-time
##########################################################################
##-- Previous calculations ------------------------------------------------
ord <- order(datos1$base_sd_T1)
baseline <- datos1$base_sd_T1[ord]
outcome <- datos1$final_sd_T1[ord]
d_figure2 <- data.frame(outcome=outcome^2,baseline=baseline^2,Order=order(baseline),ratio=(outcome/baseline)^2,
                        Group=ifelse(outcome>baseline,'Greater variance at the end of the study','Greater variance at baseline'))

##-- Arrows Graphic ------------------------------------------------
ggplot(d_figure2,aes(x=Order,y=baseline,colour=Group)) + geom_point(size=0) + 
  scale_y_log10() +
  geom_segment(aes(x = Order, y = baseline, xend = Order, yend = outcome, colour=Group),
               arrow=arrow(length = unit(0.18,"cm")),size=1) +
  theme(legend.position="bottom",
        axis.title = element_text(face='bold'),
        legend.title = element_blank()) +
  xlab('Rank according to Baseline variability') + ylab('Variance')

##-- Histogram ------------------------------------------------
ggplot(d_figure2,aes(x=ratio)) + geom_histogram(color='white') + scale_x_log10(limits=c(1/25,25)) +
  ylab('n') + xlab(expression(bold(S[OT]^2/S[BT]^2))) +
  geom_vline(xintercept = 1,linetype=2,color='white',size=1.3) +
  annotate("text",x = 0.1, y = 35, label = "paste(bold(\"Higher \"),bold(S[BT]^2))", parse = TRUE, size=8)+
  annotate("text",x = 10, y = 35, label = "paste(bold(\"Higher \"),bold(S[OT]^2))", parse = TRUE, size=8)+
  theme(axis.title = element_text(face='bold',size=13),
        axis.text = element_text(face='bold',size=13)) 

#-----------------------------------------------------------------
#
# Main Analysis --> Models table S1
#
#-----------------------------------------------------------------
source(paste0(URL,'rma_models.R'))
source(paste0(URL,'rma_models_reduced_data.R'))

#-----------------------------------------------------------------
#
# Subgroup analysis --> Figures S2, S3 and S4
#
#-----------------------------------------------------------------
source(paste0(URL,'subgroups.R'))

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
