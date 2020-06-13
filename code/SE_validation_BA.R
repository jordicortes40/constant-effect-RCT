#-------------------------------------------------------------------
#
#
#  Setting A                         
#
#
#-------------------------------------------------------------------

####################################################################
# Parameters of the simulation
####################################################################

nsim <- 100000                               # Number of simulations
N <- exp(seq(log(12),log(6144/4),length=8))  # Total sample size
ALLOCATION_RATIO <- c(1/2,1/3,1/4)           # Fraction in one group over the total sample size: allocation ratio 1:1 --> ALLOCATION_RATIO=1/2;allocation ratio 1:2 --> ALLOCATION_RATIO=1/3;allocation ratio 1:3 --> ALLOCATION_RATIO=1/4
V1 <- V2 <- 1                                # Actual variances

####################################################################
# Store information
####################################################################
M1 <- M2 <- M3 <- matrix(ncol=4,nrow=length(N)) # Store results each allocation ration
M <- list(M1,M2,M3)                             # All together ian a list

####################################################################
# Simulation
####################################################################
set.seed(12345)
for(n in N){
  for (allo_ratio in ALLOCATION_RATIO){
    
    # Print iteration
    cat('Total sample size:',n,'Allocation ratio',round(allo_ratio,2),'\n')
    
    n2 <- n*allo_ratio                       # Sample size in second group
    n1 <- n-n2                               # Sample size in first group   
    v1 <- v2 <- est1 <- c()                  # Arrays for storing information
    
    ##-- Estimate emprical Standard error E1 -----------------------------------
    for (i in 1:nsim){
      y1 <- rnorm(n1,0,sqrt(V1))             # Generate normal data for one arm
      y2 <- rnorm(n2,0,sqrt(V2))             # Generate normal data for the other arm
      v1[i] <- var(y1)                       # Sample variance for one arm
      v2[i] <- var(y2)                       # Sample variance for the other arm
      est1[i] <- log(v1[i]/v2[i])            # First estimation
    }
    
    E1 <- sd(est1)                           # Empirical standard error E1
    E_DM <- sqrt(2*(1/(n1-1)+1/(n2-1)))      # Estimation from delta method (DM): (n-1) in the denominator
    E_DMc <- sqrt(2*(1/(n1-2)+1/(n2-2)))     # Estimation from corrected delta method (DMc): (n-2) in the denominator
    
    ##-- Store information
    M[[which(ALLOCATION_RATIO==allo_ratio)]][which(N==n),] <- c(E_DM,E1,E_DM/E1,E_DMc/E1)
  }
}

#-------------------------------------------------------------------
#
#
#                         Graphic                         
#
#
#-------------------------------------------------------------------

####################################################################
# Graphic
####################################################################
##-- Summary of the real sample sizes
summary(with(datos1,final_cases_T1+final_cases_T2))
sum(with(datos1,final_cases_T1<10 | final_cases_T2<10))/208

##-- Data.frame for the plot
df0 <- as.data.table(rbind(M[[1]],M[[2]],M[[3]]))
colnames(df0) <- c('E_DM','E1','E_DM/E1','E_DMc/E1')
df0$SS <- rep(N,3)
df0$Ratio <- rep(paste0('Allocation ratio ',c('1:1','1:2','1:3')),each=8)
df1 <- melt.data.table(df0,id.vars = c('Ratio','SS'),measure.vars = c(3,4) )
mm <- matrix(unlist(strsplit(as.character(df1$variable),'/',fixed = TRUE)),ncol=2,byrow=TRUE)
df1$Estimator <- mm[,1]
df1$Comparator <- mm[,2]
write.table(x = df1,file = '../results_tables/MA_SE_validation_setting1_BA.txt',
            sep='\t',col.names = TRUE,row.names = FALSE, quote=FALSE)

ggplot(df1,aes(x=SS,y=value,colour=Estimator)) + 
  geom_hline(yintercept = 1) +
  geom_line(size=1.2) +
  scale_y_log10(limits=c(0.8,1.1),breaks=seq(0.80,1.25,0.05)) +
  scale_x_continuous(limits=c(0,769)) + # 769
  facet_wrap(.~Ratio) + # ,labeller = label_bquote(bold(Estimator: E[.(substr(Comparator,2,2))]~~~~~~~'Allocation'~.(Ratio)))
  xlab('Total sample size') +
  ylab('Ratio') +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size=13,face='bold'),
        strip.text = element_text(size=10,face='bold'),
        axis.title = element_text(size=13,face='bold')) +
  scale_color_discrete(labels = c(expression(bold(SE[DM])), expression(bold(SE[DM]^"C"))))
ggsave(file='../results_figures/MA_SE_validation_setting1_BA.jpeg',
       width=8,height = 3.5,units = 'in')
#-------------------------------------------------------------------
#
#
#                            Setting B                         
#
#
#-------------------------------------------------------------------
####################################################################
# Parameters of the simulation
####################################################################
nsim <- 10000                                # Number of simulations
N1 <- datos1$final_cases_T1                    # Sample size in experimental group
N2 <- datos1$final_cases_T2                    # Sample size in control group
V1 <- datos1$final_sd_T1^2                     # Actual variances in experimental group
V2 <- datos1$final_sd_T2^2                     # Actual variances in control group

####################################################################
# Store information
####################################################################
M <- matrix(ncol=4,nrow=208) # Matrix to store infiormation


####################################################################
# Simulation
####################################################################
set.seed(12345)
for(j in 1:208){
  n1 <- N1[j]
  n2 <- N2[j]
  v1 <- V1[j]
  v2 <- V2[j]
  cat('Iteration:',j,'of 208\n')                    # Print the process
  y01 <- rnorm(n1,0,sqrt(v1))                       # Generate normal data for one arm
  y02 <- rnorm(n2,0,sqrt(v2))                       # Generate normal data for the other arm
  
  var1 <- var2 <- var01 <- var02 <- est1 <- c()
  
  ##-- Estimate emprical Standard error E1 -----------------------------------
  for (i in 1:nsim){
    y1 <- rnorm(n1,0,sqrt(v1))                    # Generate normal data for one arm
    y2 <- rnorm(n2,0,sqrt(v2))                    # Generate normal data for the other arm
    var1[i] <- var(y1)                            # Sample variance for one arm
    var2[i] <- var(y2)                            # Sample variance for the other arm
    est1[i] <- log(var1[i]/var2[i])               # First estimation
  }
  
  E1 <- var(est1,na.rm=TRUE)                      # Empirical estimate
  E_DM <- 2/(n1-1)+2/(n2-1)                       # Estimation from delta method (DM): (n-1) in the denominator
  E_DMc <- 2/(n1-2)+2/(n2-2)                      # # Estimation from corrected delta method (DMc): (n-2) in the denominator

  ##-- Store information and update the counter
  M[j,] <- c(E_DM,E1,E_DM/E1,E_DMc/E1)
}
colnames(M) <- c('E_DM','E1','E_DMdivE1','E_DMcdivE1')

####################################################################
# Graphic with ggplot
####################################################################
df0 <- as.data.table(M)
df0$SS <- N1+N2
df2 <- melt.data.table(df0,id.vars = c('SS'),measure.vars = c(3,4))
mm <- matrix(unlist(strsplit(as.character(df2$variable),'div',fixed = TRUE)),ncol=2,byrow=TRUE)
df2$Estimator <- mm[,1]
df2$Comparator <- mm[,2]

write.table(x = df2,file = '../results_tables/MA_SE_validation_setting2_BA.txt',
            sep='\t',col.names = TRUE,row.names = FALSE, quote=FALSE)

##-- Descriptive categorized
df2$SS2 <- cut(df2$SS,c(seq(10,100,10),seq(200,900,100)),include.lowest = TRUE)
with(df2[df2$variable=='E_DMdivE1',],tapply(value,SS2,mean))


ggplot(df2,aes(x=SS,y=value,colour=Estimator)) + 
  geom_hline(yintercept = 1) +
  geom_smooth(size=1.2,se = FALSE) +
  scale_y_log10(limits=c(0.9,1.1),breaks=seq(0.9,1.1,0.05)) +
  xlim(c(0,900)) +
  xlab('Total sample size') +
  ylab('Ratio') +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size=11),
        axis.title = element_text(size=11,face='bold')) +
  scale_color_discrete(labels = c(expression(bold(SE[DM])), expression(bold(SE[DM]^C))))
ggsave(file='../results_figures/MA_SE_validation_setting2_BA.jpeg',
       width=8,height = 3.5,units = 'in')



