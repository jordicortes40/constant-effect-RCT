#-------------------------------------------------------------------
#
#
# Simulation for the standard error in over time comparison                        
#
#
#-------------------------------------------------------------------

########################################################################
##-- Parameters
########################################################################
n <- 100                                              # Study sample size
nsim <- 10000                                         # Number of simulations for the same correlation
SIGMA_EFFECT <- c(exp(seq(log(0.1),log(10),0.02)),10) # Effect deviation

########################################################################
##-- Store
########################################################################
cor.outcome <- cor.var <- c()
cov.outcome <- cov.var <- c()
sta <- c()
M <- matrix(ncol=5,nrow=length(SIGMA_EFFECT))

########################################################################
##-- Simulation
########################################################################
set.seed(12345)
for(sigma in SIGMA_EFFECT){
  vb <- c()                  # Estimated Baseline Variance
  vf <- c()                  # Estimated Final Variance
  cor.out <- cov.out <- c()  # Correlation and Covariance Baseline-Final
  for(j in 1:nsim){
    
    z <- rnorm(n)                    # Baseline values (Z)
    y <- z + rnorm(n,0,sigma)        # Final Values (Y)
    
    vb[j] <- var(z)                  # Baseline Variance
    vf[j] <- var(y)                  # Final Variance
    cor.out[j] <- cor(z,y)           # Correlation Baseline-Final
    cov.out[j] <- cov(z,y)           # Covariance Baseline-Final
  }
  i <- which(SIGMA_EFFECT==sigma)
  
  cor.outcome[i] <- mean(cor.out)                                                  # Correlation mean
  cov.outcome[i] <- mean(cov.out)                                                  # Covariance mean

  sta[i] <- log(1+2*cor.outcome[i]^2/(n-2))                                        # Deducted statistic for the covariance (right-side)
  cov.var[i] <- cov(vb,vf)                                                         # Covariance of the variance (step 1 of the proof)
  cor.var[i] <- cor(vb,vf)
  
  ##-- Store information -------------------------------
  M[i,1] <- cov(log(vf),log(vb))                                              # Covariance of the log-variance   (left-side)  --> Step 0
  M[i,2] <- log(1+cov(vf,vb)/(mean(vf)*mean(vb)))                             # After step 1
  M[i,3] <- log(1+cor.outcome[i]^2*sqrt(var(vf)*var(vb))/(mean(vf)*mean(vb))) # After step 2
  M[i,4] <- log(1+2*cor.outcome[i]^2/(n-1))                                   # After step 3
  M[i,5] <- log(mean(vf*vb))                                                  # Alternative approximation if Cov(x,y) = E(XY) - E(X)E(Y)

  ##-- Display sigma and correlation
  cat('sigma:',sigma,'correlation:',cor.outcome[i],'\n')
  
}

colnames(M) <- c('Real_covariance','Approx_covar_step_1',
                 'Approx_covar_step_2','Approx_covar_step_3',
                 'Approx_covar_alternative')
write.table(x = M,file = '../results_tables/MA_SE_validation_OT.txt',
            sep='\t',col.names = TRUE,row.names = FALSE, quote=FALSE)

##-- Descriptive of the differences
summary(M)
discrepancies <- data.frame(Real_covariance=1000*M[,'Real_covariance'],
                            abs_dif1=1000*abs(M[,'Approx_covar_step_1']-M[,'Real_covariance']),
                            abs_dif2=1000*abs(M[,'Approx_covar_step_2']-M[,'Real_covariance']),
                            abs_dif3=1000*abs(M[,'Approx_covar_step_3']-M[,'Real_covariance']))


##-- Absolute discrepancies (1)
apply(discrepancies,2,mean,na.rm=TRUE)

##-- Absolute discrepancies (2)
apply(discrepancies,2,sd,na.rm=TRUE)
##-- Relative discrepancies (1) --> There are a few negative values
mean(discrepancies[,2]/discrepancies[,1],na.rm=TRUE)
mean(discrepancies[,3]/discrepancies[,1],na.rm=TRUE)
mean(discrepancies[,4]/discrepancies[,1],na.rm=TRUE)

##- Relative discrepancies (2)
mean(discrepancies[,2],na.rm=TRUE)/mean(discrepancies[,1],na.rm=TRUE)
mean(discrepancies[,3],na.rm=TRUE)/mean(discrepancies[,1],na.rm=TRUE)
mean(discrepancies[,4],na.rm=TRUE)/mean(discrepancies[,1],na.rm=TRUE)


########################################################################
##-- Plots after each simulation
########################################################################

##-- Labels for axis
xl <- expression(Cov~bgroup("[",list(log(S[OT]^2),log(S[BT]^2)),"]"))
yl1 <- expression(Log~bgroup("[",1+frac(Corr~group("[",list(S[OT]^2,S[BT]^2),"]") %.% sqrt(V~group("[",S[OT]^2,"]") %.% V~group("[",S[BT]^2,"]")),S[OT]^2 %.% S[BT]^2),"]"))
yl2 <- expression(Log~bgroup("[",1+frac((Corr~group("[",list(Y[OT],Y[BT]),"]"))^2 %.% sqrt(V~group("[",S[OT]^2,"]") %.% V~group("[",S[BT]^2,"]")),S[OT]^2 %.% S[BT]^2),"]"))
yl3 <- expression(Log~bgroup("[",1+frac(2 %.% (Corr~group("[",list(Y[OT],Y[BT]),"]"))^2,df[OT]),"]"))  



dd0 <- as.data.table(M)
gg1 <- ggplot(dd0,aes(x=Real_covariance,y=Approx_covar_step_1)) + 
  geom_point(colour='blue',alpha=0.5) +
  xlab(xl) + ylab(yl1) +
  geom_abline(intercept=0,slope=1) +
  ggtitle('After Step 1') +
  theme(axis.title=element_text(size=8))
gg2 <- ggplot(dd0,aes(x=Real_covariance,y=Approx_covar_step_2)) + 
  geom_point(colour='blue',alpha=0.5) +
  xlab(xl) + ylab(yl2) +
  geom_abline(intercept=0,slope=1) +
  ggtitle('After Step 2') +
  theme(axis.title=element_text(size=8))
gg3 <- ggplot(dd0,aes(x=Real_covariance,y=Approx_covar_step_3)) + 
  geom_point(colour='blue',alpha=0.5) +
  xlab(xl) + ylab(yl3) +
  geom_abline(intercept=0,slope=1) +
  ggtitle('After Step 3') +
  theme(axis.title=element_text(size=8))
gg4 <- ggplot(dd0,aes(x=Real_covariance+Approx_covar_step_3,y=Real_covariance-Approx_covar_step_3)) + 
  geom_point(colour='blue',alpha=0.5) +
  xlab('Mean') + ylab('Difference') +
  geom_abline(intercept=0,slope=0) +
  ggtitle('Bland-Altman') +
  theme(axis.title=element_text(size=8))

gg_covariance <- ggarrange(gg1,gg2,gg3,gg4,nrow=2,ncol=2)
gg_covariance_def <- annotate_figure(gg_covariance,top = text_grob("Covariance approximation", face = "bold", size = 14))
ggsave(gg_covariance_def,file='../results_figures/MA_SE_validation_OT.jpeg',
       width=8,height = 7,units = 'in')



########################################################################
##-- Main plot cor.outcome versus cor.var
########################################################################
par(mfrow=c(1,1))
plot(cor.outcome,cor.var,pch=19,col=rgb(0,0,1,0.4))
lines(lowess(cor.outcome,cor.var,f=1/8),col=rgb(0,1,0,0.6),lwd=2)
curve(x^2,col=rgb(1,0,0,0.6),lwd=2,add=TRUE)
abline(0,1,lty=2)
abline(0,-1,lty=2)
abline(h=0,v=0,lty=2)

##-- Covariance plot
# plot(cov.outcome,cov.var,pch=19,col=rgb(0,0,1,0.4))
plot(sta,cov.var,pch=19,col=rgb(0,0,1,0.4))
lines(lowess(sta,cov.var),col=rgb(0,1,0,0.6),lwd=2)
abline(0,1,lty=2)
abline(0,-1,lty=2)
abline(h=0,v=0,lty=2)

