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
# Relationship between baseline and outcome variances discrepancies
# 
#
#-----------------------------------------------------------------
##-- Size of the points proportional to the square root of the mean sample sizes 
datos1$size <- sqrt(apply(datos1[,c('final_cases_T1','final_cases_T2',
                                    'base_cases_T1','base_cases_T2')],1,mean))
ggplot(datos1,aes(x=yBaselineRatio,y=yBetweenArmsRatio,size=size)) + 
  geom_point(color='darkblue',alpha=0.5) +
  geom_vline(xintercept = 0,linetype=2) +
  geom_hline(yintercept = 0,linetype=2) +
  geom_abline(intercept = 0,slope=1,linetype=2) +
  xlab(expression(bold(log~bgroup("(",frac(S[BT]^2,S[BC]^2),")")))) +
  ylab(expression(bold(log~bgroup("(",frac(S[OT]^2,S[OC]^2),")")))) +
  theme(legend.position = 'none',
        axis.title = element_text(face = 'bold',size=13),
        axis.text = element_text(face = 'bold',size=13))

ggsave('../results_figures/SA_I_final_vs_baseline_discrepancies.jpeg')


#-----------------------------------------------------------------
#
# Forestplot of individual studies: at baseline baseline
#   - At baseline with simulated data under H0 (no differences in variances)
#   - At baseline with real data
#   - At the end of the sudy with real data
#-----------------------------------------------------------------

##-- Generate data under H0
# est0 is the value of the statistic log(S^2_{BT}/S^2_{BC}) assuming the randomization works
set.seed(12345)
est0 <- c()
for(j in 1:208){
  n1 <- datos1$final_cases_T1[j]
  n2 <- datos1$final_cases_T2[j]
  v1 <- datos1$final_sd_T1[j]^2
  v2 <- datos1$final_sd_T1[j]^2
  y01 <- rnorm(n1,0,sqrt(v1))                 # Generate normal data for one arm
  y02 <- rnorm(n2,0,sqrt(v2))                 # Generate normal data for the other arm
  est0[j] <- log(var(y01)/var(y02))           # Response
}

##-- Data for graphic with simulated data under H0
ord <- order(est0)
est <- est0[ord]
se <- with(datos1,sqrt(2/(final_cases_T1-2) + 2/(final_cases_T2-2)))[ord]
LI <- est-1.96*se
LS <- est+1.96*se
dd <- data.frame(y=1:208,est=exp(est),LI=exp(LI),LS=exp(LS),significance=ifelse(LI<0 & LS>0,'NS','S'))
dd$significance <- factor(dd$significance,levels=c('S','NS')) 

##-- Model for simulated data
rma.mod1 <- rma(est,sei=se)


##-- First plot
gg1 <- ggplot(dd,aes(x=est,y=y,col=significance,xmin=LI,xmax=LS)) + geom_point() +
  scale_x_log10(limits=c(0.01,100)) + scale_y_continuous(expand = c(0, 1)) +
  geom_vline(xintercept = 1,linetype=2) +
  xlab(expression(bold(widehat(S[BT]^2~"/"~S[BC]^2)~"[95%CI]"))) + ylab('Study') + 
  geom_errorbarh() + ggtitle(expression(bold('Simulated data under'~H[0]))) +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.text=element_text(face='bold',size=13),
        axis.title=element_text(face='bold',size=13),
        title = element_text(face='bold'))

#---------------------------------------------------------------------------
##-- Data for graphic with real data at baseline
ord <- order(datos1$yBaselineRatio)
p <- datos1$yBaselineRatio[ord]
se <- with(datos1,sqrt(2/(final_cases_T1-2) + 2/(final_cases_T2-2)))[ord]
LI <- p - 1.96*se
LS <- p + 1.96*se
dd <- data.frame(y=1:208,est=exp(p),LI=exp(LI),LS=exp(LS),significance=ifelse(LI<0 & LS>0,'NS','S'))
dd$significance <- factor(dd$significance,levels=c('S','NS')) 

##-- Model for real data at baseline
rma.mod2 <- rma(p,sei=se)

##-- Second plot
gg2 <- ggplot(dd,aes(x=est,y=y,col=significance,xmin=LI,xmax=LS)) + geom_point() +
  scale_x_log10(limits=c(0.01,100)) + scale_y_continuous(expand = c(0, 1)) +
  geom_vline(xintercept = 1,linetype=2) +
  xlab(expression(bold(widehat(S[BT]^2~"/"~S[BC]^2)~"[95%CI]"))) + ylab('Study') + 
  geom_errorbarh() + ggtitle('Baseline data') +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.text=element_text(face='bold',size=13),
        axis.title=element_text(face='bold',size=13),
        title = element_text(face='bold'))

#---------------------------------------------------------------------------
##-- Data for plot of the outcome
p0 <- with(datos1,2*log(final_sd_T1/final_sd_T2))
ord <- order(p0)
p <- p0[ord]
se <- with(datos1,sqrt(2/(final_cases_T1-2) + 2/(final_cases_T2-2)))[ord]
LI <- p-1.96*se
LS <- p+1.96*se
dd <- data.frame(y=1:208,est=exp(p),LI=exp(LI),LS=exp(LS),significance=ifelse(LI<0 & LS>0,'NS','S'))
dd$significance <- factor(dd$significance,levels=c('S','NS')) 

##-- Model for real data at the end of the study
rma.mod3 <- rma(p,sei=se)

##-- Third plot
gg3 <- ggplot(dd,aes(x=est,y=y,col=significance,xmin=LI,xmax=LS)) + geom_point() +
  scale_x_log10(limits=c(0.01,100)) + scale_y_continuous(expand = c(0, 1)) +
  geom_vline(xintercept = 1,linetype=2) +
  xlab(expression(bold(widehat(S[OT]^2~"/"~S[OC]^2)~"[95%CI]"))) + ylab('Study') + 
  geom_errorbarh() + ggtitle('Outcome data') +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.text=element_text(face='bold',size=13),
        axis.title=element_text(face='bold',size=13),
        title = element_text(face='bold'))
#---------------------------------------------------------------------------

##-- Arrange plots
ggarrange(gg1,gg2,gg3,nrow=1,common.legend = TRUE,legend = "bottom")
ggsave('../results_figures/SA_I_forest_single_simulation.jpeg',width = 8)

#-------------------------------------------------------------------
#
#
# Compare the model parameters with simulation under no variance                          
#
# 
#-------------------------------------------------------------------
####################################################################
# Parameters of the simulation
####################################################################
nsim <- 10000                    # Number of simulations
N1 <- datos1$final_cases_T1        # Sample size in one group
N2 <- datos1$final_cases_T2        # Sample size in other group
V1 <- datos1$final_sd_T1^2         # Variances in treated group
V2 <- V1                         # Variances in control group=treated group

####################################################################
# Simulation
####################################################################
##-- Generate data under no differences in variances (V2=V1)
set.seed(12345)
M <- matrix(ncol=4,nrow=nsim)
for (i in 1:nsim){
  est <- c()
  cat('iteration:',i,'/',nsim,'\n')                          # Print iteration
  for(j in 1:208){
    y01 <- rnorm(N1[j],0,sqrt(V1[j]))                # Generate normal data for one arm
    y02 <- rnorm(N2[j],0,sqrt(V1[j]))                # Generate normal data for the other arm
    est[j] <- log(var(y01)/var(y02))                 # Response
  }
  try(rma.model <- rma(yi=est,vi=2/(N1-2)+2/(N2-2))) # Random effect model
  LL <- est - 1.96*sqrt(2/(N1-2)+2/(N2-2))           # Lower limit for each single study  
  UL <- est + 1.96*sqrt(2/(N1-2)+2/(N2-2))           # Upper limit for each single study  
  M[i,] <- c(coef(rma.model)[1],sqrt(rma.model$tau2),rma.model$I2,sum(LL>0 | UL<0)/208)
}
colnames(M) <- c('mu','tau','I2','p') # p: proportion of studies in each iteration whose confidence intervals does not contain 0
summary(M)
dd <- as.data.frame(M)
write.table(x = dd,file='../results_tables/SA_I_simulated_data_under_H0.txt',
            row.names = FALSE,col.names = TRUE,sep='\t')

## Quantiles
round(t(apply(dd[,1:3],2,quantile,probs=seq(0,1,0.1))),2)


##-- Boxplots for the three parameters (mu, tau, I2) in all the simulations:
rma.unadjB <- rma(yBaselineRatio,sei=seBaselineRatio,data=datos1)  
rma.unadj  <-   rma(yBetweenArmsRatio,sei=seBetweenArmsRatio,data=datos1,method='REML')              # Adjusted by baseline
rma.adj    <-   rma(yBetweenArmsRatio,sei=seBetweenArmsRatio,data=datos1,mods=~yBaselineRatio,method='REML')              # Adjusted by baseline
REAL.BASELINE <- c(rma.unadjB$beta[1],sqrt(rma.unadjB$tau2),rma.unadjB$I2)
REAL.FINAL    <- c(rma.unadj$beta[1], sqrt(rma.unadj$tau2), rma.unadj$I2)

##-- Graphical parameters
cols <- c("BASELINE"="blue","OUTCOME"="red")         # Colors for the plot
sha <- c("BASELINE"=3,"OUTCOME"=4)                   # Shapes for the plot

# Common theme
common.theme <- theme(axis.text.x = element_blank(),
                      axis.ticks.x = element_blank(),
                      axis.title = element_text(size=15,face='bold'),
                      legend.text = element_text(face='bold'))
##-- First boxplot: mu
gg1 <- ggplot(dd,aes(y=mu)) + geom_boxplot() +
  geom_point(aes(x=0, y=REAL.BASELINE[1],colour="BASELINE",shape="BASELINE"),size=3,stroke =2) +
  geom_point(aes(x=0, y=REAL.FINAL[1],colour="OUTCOME",shape="OUTCOME"), size=3,stroke =2) +
  xlab('') + ylab(expression(bold(hat(mu)))) + common.theme +
  scale_colour_manual(name="",values=cols) + scale_shape_manual(name="",values=sha)

##-- Second boxplot: tau
gg2 <- ggplot(dd,aes(y=tau)) + geom_boxplot() +
  geom_point(aes(x=0, y=REAL.BASELINE[2]), colour="blue",size=3,shape=3,stroke =2) +
  geom_point(aes(x=0, y=REAL.FINAL[2]), colour="red",size=3,shape=4,stroke =2) +
  xlab('') + ylab(expression(bold(hat(tau)))) + common.theme     # "\u03c4"

##-- Third boxplot: I2
gg3 <- ggplot(dd,aes(y=I2)) + geom_boxplot() +
  geom_point(aes(x=0, y=REAL.BASELINE[3]), colour="blue",shape=3,size=3,stroke =2) +
  geom_point(aes(x=0, y=REAL.FINAL[3]), colour="red",shape=4,size=3,stroke =2) +
  xlab('') + ylab(expression(bold(I^2))) + common.theme

##-- Arrange boxplots
ggarrange(gg1,gg2,gg3,nrow=1, common.legend = TRUE,legend = 'bottom')
ggsave('../results_figures/SA_I_parameters_distribution_simulation.jpeg')

##-- Table to compare parameters
dd1 <- data.frame(mu= c(mean(dd$mu), REAL.BASELINE[1],REAL.FINAL[1]),
                  tau=c(mean(dd$tau),REAL.BASELINE[2],REAL.FINAL[2]),
                  I2= c(mean(dd$I2), REAL.BASELINE[3],REAL.FINAL[3]))
round(dd1,2)


#-----------------------------------------------------------------
#
# Exclusions based on discrepancies
#
#-----------------------------------------------------------------
##-- Tau_target is the (arbitrary) quantile 0.90 of the simulated distribution in a non-heterogeneity scenario.
tau_target <- quantile(dd$tau,0.90)

##-- Remove studies up to heterogeneity was close to tau target
limit_ratio_baseline <- uniroot(limit_absolute_ratio,interval=c(0.001,100),
                                data=datos1,
                                y=datos1$yBaselineRatio,
                                se=datos1$seBaselineRatio,comparison='Baseline',
                                quantile_tau=tau_target)$root
limit_ratio_BA <- uniroot(limit_absolute_ratio,interval=c(0.001,100),
                                data=datos1,
                                y=datos1$yBetweenArmsRatio,
                                se=datos1$seBetweenArmsRatio,comparison='BA',
                                quantile_tau=tau_target)$root
limit_ratio_OT <- uniroot(limit_absolute_ratio,interval=c(0.001,100),
                          data=datos1,
                          y=datos1$yOverTimeRatioT,
                          se=datos1$seOverTimeRatioT,comparison='OT',
                          quantile_tau=tau_target,extendInt='yes')$root
##-- Exclusion variable to achieve a heterogeneity parameter (tau) around 0.07-0.08
Exc_Basal       <- with(datos1,abs(yBaselineRatio/seBaselineRatio)>=limit_ratio_baseline)
Exc_BetweenArms <- with(datos1,abs(yBetweenArmsRatio/seBetweenArmsRatio)>=limit_ratio_BA)
Exc_Overtime    <- with(datos1,abs(yOverTimeRatioT/seOverTimeRatioT)>=limit_ratio_OT)

##-- Reduced datasets
dataB        <- datos1[!Exc_Basal,]                                           # Reduced data for baseline model
dataBetArm   <- datos1[!Exc_BetweenArms,]                                     # Reduced data for between arms models
dataOverTime <- datos1[!Exc_Overtime & !is.na(datos1$seOverTimeRatioT),]      # Reduced data for over-time models

cat('Reduced dataset at baseline (n=',nrow(dataB),')\n',sep = '')
cat('Reduced dataset between arms (n=',nrow(dataBetArm),')\n',sep = '')
cat('Reduced dataset over time (n=',nrow(dataOverTime),')\n',sep = '')

#-----------------------------------------------------------------
#
# Sensitive Analysis I --> All Models. This section produces all model outputs: 
#     - On reduced datasets
#
#-----------------------------------------------------------------
source(paste0(URL,'code/rma_models_reduced_data.R'),local=TRUE)

#-----------------------------------------------------------------
#
# Data for summary table
#
#-----------------------------------------------------------------
##-- Between arms
SAI_greater_BA <- sum(Exc_BetweenArms & datos1$yBetweenArmsRatio>0)
SAI_lower_BA <- sum(Exc_BetweenArms & datos1$yBetweenArmsRatio<0)
SAI_equal_BA <- nrow(datos1) - SAI_lower_BA - SAI_greater_BA
SAI_greater_BA; SAI_lower_BA; SAI_equal_BA

##-- Over time
SAI_greater_OT <- sum(Exc_Overtime & datos1$yOverTimeRatioT>0,na.rm=TRUE)
SAI_lower_OT <- sum(Exc_Overtime & datos1$yOverTimeRatioT<0,na.rm=TRUE)
SAI_equal_OT <- sum(!is.na(datos1$seOverTimeRatioT)) - SAI_lower_OT - SAI_greater_OT
SAI_greater_OT; SAI_lower_OT; SAI_equal_OT
