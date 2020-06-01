rm(list=ls())

#-----------------------------------------------------------------
#
# Read data
#
#-----------------------------------------------------------------
URL <- 'https://raw.githubusercontent.com/jordicortes40/constant-effect-RCT/master/'
source(paste0(URL,'code/read_data.R')) 
d <- datos1

####################################################################
# Information from data
####################################################################
N <- nrow(d)                   # Number of rows
n_treat <- d$final_cases_T1    # Sample size treated
n_ctrl <- d$final_cases_T2     # Sample size controls
n <- n_treat + n_ctrl          # Total sample size

#-----------------------------------------------------------------
#
# Simulation under main scenario:
# - All Additive effect with (100xpi_r)% Random and 100*(1-pi_r)% Fixed
#
#-----------------------------------------------------------------
##-- Parameters
THETA_MAX <- c(1,3,5,7)          # theta=sigma_C/sigma_T; theta_max = max(theta)
PROP_RANDOM <- seq(0,0.5,0.05)   # proportion of studies with random treatment effect
DELTA_C <- 0                     # Change in controls over time
DELTA_T <- 0                     # Change in treated over time 
EFFECT <- DELTA_T - DELTA_C      # Effect
nsim <- 100                      # number of simulations

##-- Martix to store the results
M <- matrix(nrow=nsim*length(THETA_MAX)*length(PROP_RANDOM),ncol=9)
colnames(M) <- c('theta','prop_random','iteration',
                 'mu','tau','I2',        # Between arms
                 'mu_2','tau_2','I2_2')  # Over time
#                 'p_greater','p_lower')  # proportion of studies with random effects and greater/lower variability in the 

##-- Simulations
set.seed(12345)
rowi <- 1                    # Row indicator
for(pa in PROP_RANDOM){
  for(theta in THETA_MAX){
    for(k in 1:nsim){
      y <- se <- yBaselineRatio <- c()
      y2 <- se2 <- yControlRatio <- c()
      
      ##-- Generate data for each single study
      for(j in 1:N){
        effect.type   <- sample(c('fixed','random'),1,prob=c(1-pa,pa))   # generate type of treatment effect: random or fixed
        sigma.control <- ifelse(effect.type=='fixed',0,runif(1,0,theta)) # generate standard deviation of the effect in control group. If fixed, sigma=0, otherwise sigma=runif(0,sigma_max)
        sigma.treated <- as.numeric(effect.type=='random')               # generate standard deviation of the effect in treated group. If fixed, sigma=0, otherwise sigma=1
        
        Y_B <- rnorm(n[j],0,1)                            # generate baseline values
        Y_OC <- Y_B + rnorm(n[j],DELTA_C,sigma.control)   # generate potential outcomes in control arm
        Y_OT <- Y_B + rnorm(n[j],DELTA_T,sigma.treated)   # generate potential outcomes in treatment arm
        
        sel.trt <- sample(1:n[j],n_treat[j],rep=FALSE)                    # indicator of treated patients 
        sel.ctrl <- (1:n[j])[-sel.trt]                                    # indicator of reference patients
        y_bc <- Y_B[sel.ctrl]                                                # observed outcome of treated patients
        y_bt <- Y_B[sel.trt]                                                 # observed outcome of treated patients
        y_oc <- Y_OC[sel.ctrl]                                               # observed outcome of treated patients
        y_ot <- Y_OT[sel.trt]                                                # observed outcome of treated patients
        
        # Between arms comparison
        y[j] <- log(var(y_ot)/var(y_oc))                                     # response of the model
        yBaselineRatio[j] <- log(var(y_bt)/var(y_bc))                        # response of the model at baseline
        se[j] <- sqrt(2/(n_treat[j]-2) + 2/(n_ctrl[j]-2))                    # within standard error of the response
        
        # Over time comparison
        y2[j] <- log(var(y_ot)/var(y_bt))                                    # response of the model
        yControlRatio[j] <- log(var(y_oc)/var(y_bc))                         # response of the model at baseline
        se2[j] <- sqrt(2/(n_treat[j]-2) + 
                         2/(n_treat[j]-2) - 
                         2*log(1 + 2*cor(y_ot,y_bt)^2/(n_treat[j]-1)))       # within standard error of the response
        
        
      }
      
      # Between arms
      data <- data.frame(y=y,se=se,yBaselineRatio=yBaselineRatio)                        # data to fit the model
      mod <- try(rma(y,sei=se,data=data,mods=~yBaselineRatio,method='REML'),silent=TRUE) # fit the model                                     
      
      # Over time
      data2 <- data.frame(y=y2,se=se2,yBaselineRatio=yControlRatio)                        # data to fit the model
      mod2 <- try(rma(y2,sei=se2,data=data2,mods=~yControlRatio,method='REML'),silent=TRUE) # fit the model                                     
      
      
      M[rowi,1] <- theta                     # THETA_MAX                            
      M[rowi,2] <- pa                        # proportion of random studies
      M[rowi,3] <- k                         # Iteration
      
      # Between arms
      if(class(mod)[1]!="try-error"){
        ##-- Store results
        M[rowi,4] <- as.numeric(mod$beta[1,1]) # coefficient that estimates log(mu)
        M[rowi,5] <- sqrt(mod$tau2)            # estimate of tau^2 
        M[rowi,6] <- mod$I2                    # estimate of I^2 
      }else{
        cat('Fisher scoring algorithm did not converge in between arms comparison.\n')
      }
      
      # Over time
      if(class(mod2)[1]!="try-error"){
        ##-- Store results
        M[rowi,7] <- as.numeric(mod2$beta[1,1]) # coefficient that estimates log(mu)
        M[rowi,8] <- sqrt(mod2$tau2)            # estimate of tau^2 
        M[rowi,9] <- mod2$I2                    # estimate of I^2 
      }else{
        cat('Fisher scoring algorithm did not converge in over time comparison.\n')
      }
      
      ##-- Update indicator
      rowi <- rowi+1
    }
    cat('Theta_max:',theta,'Prop. Random:',pa,'\n')
  }
}

##-- Summary of the results
summary(M)

##-- Store the data
df <- as.data.frame(M)
write.table(x = df,file='../results_tables/SA_II_simulated_data.txt',
            row.names = FALSE,col.names = TRUE,sep='\t')


##-- Arrange the data
df$I2 <- df$I2/100
df$I2_2 <- df$I2_2/100
df2 <- as.data.table(df)[,.(mu=mean(mu,na.rm=TRUE),sd.mu=sd(mu,na.rm=TRUE),
                            tau=mean(tau,na.rm=TRUE),sd.tau=sd(tau,na.rm=TRUE),
                            I2=mean(I2,na.rm=TRUE),sd.I2=sd(I2,na.rm=TRUE),
                            mu2=mean(mu_2,na.rm=TRUE),sd.mu2=sd(mu_2,na.rm=TRUE),
                            tau2=mean(tau_2,na.rm=TRUE),sd.tau2=sd(tau_2,na.rm=TRUE),
                            I22=mean(I2_2,na.rm=TRUE),sd.I22=sd(I2_2,na.rm=TRUE)),
                         by=.(theta,prop_random)]

## Data frame for plot between arms
df3 <- melt(df2,id.vars=1:2, measure.vars=c(3,5,7),variable.name = "stat")
df3[,theta:=factor(theta)]
df3[,theta_string:=paste0('theta[M]==',theta)]
df3$stat <- factor(df3$stat,levels=c('mu','tau','I2'))

## Data frame for plot over time (does not work because a specified correlation has to be fixed)
df4 <- melt(df2,id.vars=1:2, measure.vars=c(9,11,13),variable.name = "stat")
df4[,theta:=factor(theta)]
df4[,theta_string:=paste0('theta[M]==',theta)]

##-- Between arms estimation
rma.adj <-   rma(yBetweenArmsRatio,sei=seBetweenArmsRatio,data=datos1,mods=~yBaselineRatio,method='REML')              # Adjusted by baseline
REAL.FINAL <- c(rma.adj$beta[1],sqrt(rma.adj$tau2),rma.adj$I2/100)

##-- Between arms plot
ggplot(df3,aes(x=prop_random,y=value,color=stat)) + 
  stat_summary(fun.y=mean,geom='point',size=1.5) +
  stat_summary(fun.y=mean,geom='line',size=1) +
  facet_wrap(~theta_string,nrow = 2,labeller = label_parsed) +
  scale_x_continuous(limits=c(0,0.5)) +
  geom_hline(yintercept = REAL.FINAL,linetype=2,color=rep(gg_color_hue(3),length(THETA_MAX))) + 
  xlab(expression('Proportion ('~pi[R]~') of studies with random effect')) +
  labs(title = expression('All studies with additive treatment effect:'~pi[R]~'of them with random effect')) + # subtitle = expression("Fixed effect:"~Y[OT[i]]==Y[OC[i]]~~~~~"Random effect:"~Y[OT[i]]==Y[OC[i]] + T~~~~~~~~T~"~"~N(0,tau^2)~~~~~~~~tau^2~"~ ["~0~","~theta[M]~"]")
  theme(legend.title = element_blank(),legend.position = 'bottom') +
  ylab('Value') +
  scale_color_manual(labels = c(expression(~~mu),expression(~~tau),expression(~~I^2)),values=gg_color_hue(3))

ggsave(filename = '../results_figures/SA_II_simulation.jpg')

##-- Between arms plot --> Does not work
ggplot(df4,aes(x=prop_random,y=value,color=stat)) + 
  stat_summary(fun.y=mean,geom='point',size=1.5) +
  stat_summary(fun.y=mean,geom='line',size=1) +
  facet_wrap(~theta_string,nrow = 2,labeller = label_parsed) +
  scale_x_continuous(limits=c(0,0.5)) +
  geom_hline(yintercept = c(-0.1510,0.5897,0.9365),linetype=2,color=rep(gg_color_hue(3),length(THETA_MAX))) + 
  xlab(expression('Proportion ('~pi[R]~') of studies with random effect')) +
  labs(title = expression('All studies with additive treatment effect:'~pi[R]~'of them with random effect')) + # subtitle = expression("Fixed effect:"~Y[OT[i]]==Y[OC[i]]~~~~~"Random effect:"~Y[OT[i]]==Y[OC[i]] + T~~~~~~~~T~"~"~N(0,tau^2)~~~~~~~~tau^2~"~ ["~0~","~theta[M]~"]")
  theme(legend.title = element_blank(),legend.position = 'bottom') +
  ylab('Value') +
  scale_color_manual(labels = c(expression(~~mu),expression(~~tau),expression(~~I^2)),values=gg_color_hue(3))


#-----------------------------------------------------------------
#
# Data for summary table
#
#-----------------------------------------------------------------
##-- Between arms
PI_R <- 0.10
THETA_M <- 7
p_greater_BA <- PI_R*1/THETA_M
p_lower_BA <- PI_R*(THETA_M-1)/THETA_M
SAII_greater_BA <- round(p_greater_BA*nrow(datos1))
SAII_lower_BA <- round(p_lower_BA*nrow(datos1))
SAII_equal_BA <- nrow(datos1) - SAII_lower_BA - SAII_greater_BA
SAII_greater_BA; SAII_lower_BA; SAII_equal_BA

