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
# - Independent samples (between arms): F test
# - Paired samples (over time): Q test
#
#-----------------------------------------------------------------
##########################################################################
##-- Between arms
##########################################################################
##-- Based on F test
Fest <- with(datos1,(final_sd_T1/final_sd_T2)^2)
LL <- with(datos1,qf(0.025,final_cases_T1-1,final_cases_T2-1))
UL <- with(datos1,qf(0.975,final_cases_T1-1,final_cases_T2-1))
p_values_BA <- data.frame(p_values_BA=pf(Fest,
                                         datos1$final_cases_T1-1,
                                         datos1$final_cases_T2-1))
write.table(x = p_values_BA,file='../results_tables/SA_III_pvalues_BA.txt',
            row.names = FALSE,col.names = TRUE,sep='\t')
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
p_values_OT <- data.frame(p_values_OT=pt(Qest,n-2))
write.table(x = p_values_OT,file='../results_tables/SA_III_pvalues_OT.txt',
            row.names = FALSE,col.names = TRUE,sep='\t')

#-----------------------------------------------------------------
#
# Funnel plots
#
#-----------------------------------------------------------------
############################################################
# Funnel plot between arms (vertical)
############################################################
##-- Previous calculations
y2 <- datos1$yBetweenArmsRatio                            # Between arms variance ratio
w <- with(datos1,sqrt((final_cases_T1+final_cases_T2)))   # Uncertainty
colbg <- 'grey80'                                         # color for background
co1 <- 1:2                                                # color for points
LL <- with(datos1,qf(0.025,final_cases_T1,final_cases_T2))# Lower limit
UL <- with(datos1,qf(0.975,final_cases_T1,final_cases_T2))# Upper limit
sign <- Fest<LL | Fest>UL                                 # Significant of the between arms analysis

##-- Plot
# Initialize
jpeg('../results_figures/SA_III_funnel_BA.jpeg',width=27,height=21,units = 'cm',res=144)
par(las=1,mfrow=c(1,1),mar=c(6,7,2,1),font.axis=4,font.lab=2,mgp=c(4.2,1,0))
yl <- bquote(bold(Uncertainty~~~bgroup("(",bold(frac(1,sqrt(n[OT]+n[CT]))),")"))) # 'Uncertainty'
xl <- bquote(bold(frac(S[OT]^2,S[CT]^2)))
# Base
plot(y2,1/w,pch=19,xlab=xl,ylab=yl,col=0,xlim=c(log(0.01),log(100)),ylim=c(0.42,0.03),xaxt='n') 
rect(-10,-10,100,300,col=colbg)
ticks <- c(0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100)
axis(1,at=log(ticks),lab=ticks)
ticks2 <- c(0.1,0.2,0.3,0.4)
abline(h=ticks2,lwd=2,col='white')

# Linear models for non-rejection region
xL <- xU <- 1/w
yL <- log(LL)
yU <- log(UL)
modL <- lm(xL~yL)  # Model lower
modU <- lm(xU~yU)  # Model upper

# Plot region
xpoly <- c(log(0.01),log(1),log(100))
ypoly <- as.numeric(c(predict(modL,data.frame(yL=xpoly[1:2])),predict(modU,data.frame(yU=xpoly[3]))))
polygon(x=xpoly,y=ypoly,col='white',border=1,lty=3)
xaxis <- 0.42+ 0.04*(0.42-0.03)
segments(log(0.01),xaxis,log(100),xaxis,col=1,lwd=1,lty=1)
abline(v=0,col=1,lwd=1,lty=1)
points(y2,1/w,pch=19,col=co1[sign+1])

# Labels
mtext("Greater Treated",1,adj=1,at=log(100),line=2.5,cex=.8,font=2,las=1)
mtext("Arm Variability",1,adj=1,at=log(100),line=3.3,cex=.8,font=2,las=1)
mtext("Greater Control",1,adj=0,at=log(0.01),line=2.5,cex=.8,font=2,las=1)
mtext("Arm Variability",1,adj=0,at=log(0.01),line=3.3,cex=.8,font=2,las=1)

dev.off()

############################################################
# Funnel plot over-time (vertical) --> Function of variance ratio
############################################################
##-- Previous calculations
sel_OT <- !is.na(datos1$rho) & datos1$rho<1 & datos1$rho>-1
datos2 <- datos1[sel_OT,]
n.datos2 <- nrow(datos2)
LS2 <- LI2 <- c()
for (i in 1:n.datos2){      # Calculate limits
  n <- datos2$final_cases_T1[i]
  Qy <- datos2$base_sd_T1[i]^2*(n-1)
  Qxy <- with(datos2,rho[i]*base_sd_T1[i]*final_sd_T1[i]*(n-1))
  
  roo <- uniroot(QLIM1,c(Qy/100,100*Qy),Qy=Qy,Qxy=Qxy,n=n)$root 
  LS2[i] <- sqrt(roo/(n-1))
  
  roo <- try(uniroot(QLIM2,c(Qy/1000,100*Qy),Qy=Qy,Qxy=Qxy,n=n)$root)
  if(class(roo)!='try-error') LI2[i] <- sqrt(roo/(n-1))
  
  print(i)
}

y2 <- datos2$yOverTimeRatioT
y_pos <- y2>=0
y_neg <- y2<0
x2_1 <- log(LS2/datos2$base_sd_T1)[y_pos] #1/sqrt(nest) # LS2/datos2$base_sd_T1
x2_2 <- -log(LI2/datos2$base_sd_T1)[y_neg]
x2 <- 1/sqrt(datos2$base_cases_T1)
 
# Significant studies
sign2 <- Qest<LL2 | Qest>UL2

##-- Plot
# Initialize
jpeg('../results_figures/SA_III_funnel__variance_ratio_OT.jpeg',width=27,height=21,units = 'cm',res=144)
yl <- bquote(bold(Uncertainty~~~bgroup("(",bold(frac(1,sqrt(n[BT]))),")"))) # 'Uncertainty'
xl <- bquote(bold(frac(S[OT]^2,S[BT]^2)))
par(las=1,mfrow=c(1,1),mar=c(6,7,2,1),font.axis=4,font.lab=2,mgp=c(4.2,1,0))

# Base
plot(y2,x2,pch=19,col=sign2+1,xlim=c(log(0.01),log(100)),ylim=c(0.42,0.03),xaxt='n',
     xlab=xl,ylab=yl) # bquote(bold(log~bgroup("(",frac(SD[O],SD[B]),")")))
ticks <- c(0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100)
axis(1,at=log(ticks),lab=ticks)
rect(-10,-10,100,300,col=colbg)
abline(h=seq(0.1,0.4,0.1),lwd=2,col='white')

# Plot region
a1 <- log(LS2/datos2$base_sd_T1)[y_pos]*2
a2 <- log(LI2/datos2$base_sd_T1)[y_neg]*2
xpoly <- c(log(0.01),log(1),log(100))
ypoly <- -as.numeric(c(predict(lm(x2_1~a1),data.frame(a1=xpoly[1:2])),predict(lm(x2_2~a2),data.frame(a2=xpoly[3]))))
polygon(x=xpoly,y=ypoly,col='white',border=1,lty=3)
xaxis <- 0.42+ 0.04*(0.42-0.03)
segments(log(0.01),xaxis,log(100),xaxis,col=1,lwd=1,lty=1)
abline(lm(x2_1~a1),lty=3,lwd=1,col=1)
abline(lm(x2_2~a2),lty=3,lwd=1,col=1)
abline(v=0,lty=1,lwd=1,col=1)
points(y2[y_pos],x2[y_pos],pch=19,col=co1[sign2+1][y_pos])
points(y2[y_neg],x2[y_neg],pch=19,col=co1[sign2+1][y_neg])

# Labels
mtext("Greater Outcome",1,adj=1,at=log(100),line=2.5,cex=.8,font=2,las=1)
mtext("Variability",1,adj=1,at=log(100),line=3.3,cex=.8,font=2,las=1)
mtext("Greater Baseline",1,adj=0,at=log(0.01),line=2.5,cex=.8,font=2,las=1)
mtext("Variability",1,adj=0,at=log(0.01),line=3.3,cex=.8,font=2,las=1)

dev.off()

############################################################
# Funnel plot over-time (vertical) --> Function of the statistic
############################################################
##-- Previos calculations
y2 <- 1/sqrt(datos2$base_cases_T1)
x2 <- -Qest

##-- Plot
jpeg('../results_figures/SA_III_funnel_Q_statistic_OT.jpeg',width=27,height=21,units = 'cm',res=144)
par(mgp=c(3,1,0))
plot(x2,y2,col=co1[sign2+1],pch=19,xlim=c(-10,10),ylim=c(0.35,0.05),ylab=bquote(bold("Uncertainty"~~~bgroup("(",frac(1,sqrt(n[BT])),")"))),xlab='Q')
rect(-15,0,15,1,col=colbg)
abline(h=seq(0.05,0.35,0.05),col='white',lwd=2)

# Plot region
ypoly <- c(0,.4,.4,0)
xpoly <- as.numeric(c(predict(lm(LL2~x2),data.frame(x2=ypoly[1:2])),predict(lm(UL2~x2),data.frame(x2=ypoly[3:4]))))
polygon(x=xpoly,y=ypoly,col='white',border=1,lty=3)
xaxis <- 0.35 + 0.04*(0.35-0.05)
segments(-10,xaxis,10,xaxis,col=1,lwd=1,lty=1)
abline(v=0,lty=1,col=1,lwd=1)
points(x2,y2,col=co1[sign2+1],pch=19)

# Labels
mtext("Greater Outcome",1,adj=1,at=10,line=2.5,cex=.8,font=2,las=0)
mtext("Variability",1,adj=1,at=10,line=3.3,cex=.8,font=2,las=0)
mtext("Greater Baseline",1,adj=0,at=-10,line=2.5,cex=.8,font=2,las=0)
mtext("Variability",1,adj=0,at=-10,line=3.3,cex=.8,font=2,las=0)

dev.off()



#-----------------------------------------------------------------
#
# Data for summary table
#
#-----------------------------------------------------------------
##-- Between arms
SAIII_greater_BA <- sum(Fest > UL)
SAIII_lower_BA <- sum(Fest < LL)
SAIII_equal_BA <- nrow(datos1) - SAIII_lower_BA - SAIII_greater_BA
SAIII_greater_BA; SAIII_lower_BA; SAIII_equal_BA

##-- Over time
SAIII_greater_OT <- sum(Qest > UL2) 
SAIII_lower_OT  <- sum(Qest < LL2)
SAIII_equal_OT <- sum(!is.na(datos1$seOverTimeRatioT)) - SAIII_lower_OT - SAIII_greater_OT
SAIII_greater_OT; SAIII_lower_OT; SAIII_equal_OT
