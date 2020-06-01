##########################################################################
# myForest. Perform a forest-plo for subgroups
# M: Matrix with data in four columns (Est, LL, UL, n)
# xl: X-axis lab
# lab1: categories/levels of subgroups (with pertinent NAs in gaps)
# lab2: names of variables (with pertinent NAs in gaps)
# tit: title
# laxis1: label for both sides of x-axis
##########################################################################
myForest <- function(M,xl='',lab1,lab2,tit='',laxis1=''){

  ##-- Introduce NA into the data for empty rows
  N <- rep(NA,4)
  M <- rbind(M[1,],N,
             N,M[2:3,],
             N,M[4:5,],
             N,M[6:7,],
             N,M[8:9,],
             N,M[10:11,])
  
  ##-- Plot region
  ymax <- length(lab1)
  xmax <- 2
  par(las=1,mar=c(5,12,5,10))
  plot(NA,ylim=c(0.5,ymax+.5),xlim=c(1/xmax,xmax),bty='n',yaxt='n',ylab='',xlab='',log='x',xaxt='n',main=tit,
       xaxs='i')
  par(xpd=NA)
  rect(0.001,seq(ymax-1.5,3.5,-6),1000,seq(ymax-4.5,.5,-6),border=NA,col='grey90')
  rect(0.001,seq(ymax-4.5,6.5,-6),1000,seq(ymax-7.5,3.5,-6),border=NA,col='grey80')
  rect(0.001,ymax+1,1000,ymax-1.5,border=NA,col='grey80')
  par(xpd=FALSE)
  tck <- ltck <- round(exp(seq(log(0.5),log(2),le=11)),1)
  abline(v=tck,col='white',lwd=2)
  abline(v=1,col='black',lwd=2)
  axis(1,at=tck,lab=ltck,font=4,cex.axis=0.9)
  
  ##-- Labels
  lab2[1] <- paste0(lab2[1],' (n = 208)')
  mtext(lab2,2,at=ymax:1,line=12,font=2,cex=1,adj=0)

  ##-- Text for confidence intervals
  for(i in ymax:1){
    row <- ymax-i+1
    if(!is.na(M[row,1])){
      mtext(paste0(formatC(as.numeric(M[row,1]),digits=2,format='f'),'  [',
                   formatC(as.numeric(M[row,2]),digits=2,format='f'),' , ',
                   formatC(as.numeric(M[row,3]),digits=2,format='f'),']'),4,at=i,line=5,font=3,cex=0.9,adj=0.5)
      if(i!=ymax) mtext(paste0(lab1[row],' (n = ',M[row,4],')'),2,at=i,line=11,font=3,cex=0.9,adj=0)
    }
  }
  mtext(expression(bold(hat(bold('\u03BC')))),4,at=ymax+2.25,cex=1,adj=1,line=3.6)
  mtext(' [CI95%]',4,at=ymax+2.25,cex=1,adj=0,line=4,font=2)

  ##-- Rectangle for overall estimation
  x0 <- c(M[1,2],M[1,1],M[1,3],M[1,1])
  y0 <- c(ymax,ymax-0.4,ymax,ymax+0.4)
  polygon(x0,y0,col='darkblue',border=NA)
  
  ##-- Confidence intervals
  segments(M[-1,2],(ymax-1):1,M[-1,3],(ymax-1):1,col='darkblue',lwd=2,lty=1)
  points(M[-1,1],(ymax-1):1,col='darkblue',pch=15,cex=1.4)
  
  ##-- Labels in x-axis
  mtext(c('Greater variability','Greater variability'),1,at=c(0.5,2),adj=c(0,1),line=2.2,font=2,cex=0.9)
  mtext(laxis1,1,at=c(0.5,2),adj=c(0,1),line=3,font=2,cex=0.9)
  mtext(xl,1,at=1,adj=0.5,line=3.7,font=2,cex=1)
  
}

##########################################################################
# fd. descriptive of a quantitative variable
# x: a continuous variables
##########################################################################
fd <- function(x){
  return(c(sum(!is.na(x)),mean(x),sd(x),as.numeric(quantile(x,seq(0,1,0.25)))))
}

##########################################################################
# figure1. descriptive of a quantitative variable
# ctrl: values for controls
# trt: values for treateds
# log: logical. If log transformation should be performed
##########################################################################

figure1 <- function(ctrl,trt,log){
  
  ##-- Number of patients
  n <- length(trt)
  
  ##-- Plot region
  graphics.off()
  windows(12,5)
  par(las=1,font.axis=4,font.lab=2,mar=c(5,5,1,1))
  
  ##-- Colors
  co0 <- c(rgb(1,0,0,0.5),rgb(0,0,1,0.8))
  colbg <- 'grey90'
  co <- ifelse(ctrl>trt,co0[1],co0[2])
  
  ##-- Log transformation
  if(!log){
    yl <- 'Standard Deviation'
    logy <- 'y'
    yminbg <- 0.01
    wlines <- c(0.1,1,10,100,1000)
  }else{
    yl <- 'log(SD)'
    logy <- ''
    yminbg <- -10
    wlines <- seq(-2,6,2)
  }
  
  ##-- Plot
  plot(ctrl,pch=1,col='grey',ylab=yl,main='',log=logy,
       xlab='Control outcome SD rank')
  rect(-10,yminbg,300,10000,col=colbg)
  abline(h=wlines,col='white',lwd=2)
  arrows(1:n,ctrl,1:n,trt,lwd=2,col=co,le=0.05)
  legend('topleft',c('SD greater in Controls','SD greater in Treateds'),
         col=co0,lwd=2,text.font = 2,bg=NA,pt.lwd = 2,cex=1.2)
}

##########################################################################
# var.paired.test. statistic for variance paired test
# sdx: standard deviation in one group
# sdy: standard deviation in the other group
# n: number of cases
# rho: correlation between variables
##########################################################################
var.paired.test <- function(sdx,sdy,n,rho){
  Qx <- sdx^2*(n-1)            
  Qy <- sdy^2*(n-1)
  Qxy <- rho*sdx*sdy*(n-1)
  num <- (Qx-Qy)*sqrt(n-2)          # Statistic numerator
  den <- 2*sqrt(Qx*Qy-Qxy^2)        # Statistic denominator
  Qest <- num/den
  return(list(Qest,num,den))
}

##########################################################################
# myFunnel. Pseudo funnel plot
# data: standard deviation in one group
# comparison: 'Between arms','Baseline between arms' or 'Over time'
# model: number of cases
# subgroup: point colors according differences in 'variance', 'effect' or no colors (NA)
# zoom: logical. TRUE if axis are fitted to data
# xlab: label for x-axis
# ...: other parameters for plot function
##########################################################################
myFunnel <- function(data,comparison,model,subgroup=NA,zoom=TRUE,xlab,...){
  
  ##-- Parameters according tipe of graphic
  if(comparison=='Between arms'){
    x=data$yBetweenArmsRatio
    y=data$seBetweenArmsRatio
    tit='Between arms'
    lab.11 <- 'Greater Treated'; lab.12 <- 'Greater Control'
    lab.21 <- 'Arm Variability'; lab.22 <- 'Arm Variability'
  }
  if(comparison=='Over time'){
    x=data$yOverTimeRatioT
    y=data$seOverTimeRatioT
    tit='Over time'
    lab.11 <- 'Greater Outcome';lab.12 <- 'Greater Baseline'
    lab.21 <- 'Variability';     lab.22 <- 'Variability'
  }
  if(comparison=='Baseline between arms'){
    x=data$yBaselineRatio
    y=data$seBaselineRatio
    tit='Baseline between arms'
    lab.11 <- 'Greater Treated'; lab.12 <- 'Greater Control'
    lab.21 <- 'Arm Variability'; lab.22 <- 'Arm Variability'
  }

  ##-- Plot region
  graphics.off()
  windows()
  if(zoom){xmax <- max(abs(exp(x)));ymax <- 1.02;ymin <- 0.05}else{xmax <- 100;ymax <- 1.1;ymin <- 0}
  par(las=1)
  plot(x,y,pch=19,xlab='',col=0,xaxt='n',xlim=log(c(1/xmax,xmax)),ylim=c(ymax,ymin),main=tit,...)
  mtext(xlab,1,at=0,line=4,adj=0.5)
  rect(-10,-10,300,100,col='grey85')
  ticksy <- seq(0,1,0.2)
  abline(h=ticksy,lwd=2,col='white')
  x0 <- c(4,0,-4,4) 
  y0 <- c(2,0,2,2)
  polygon(x0,y0,col='white',lty=3)
  abline(v=0)
  abline(v=coef(model)[1],lwd=2,lty=2,col=4)
  ticksx <- c(0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100)
  axis(1,at=log(ticksx),lab=ticksx)
  
  ##-- Point colors
  if(is.na(subgroup)){
    co <- 1
    pc <- 1
  }else{
    if(subgroup=='variance'){
      ##-- Significant studies regarding to variance
      sign1 <- x < (-2*y)
      sign2 <-  x > 2*y
      sign <- sign1 | sign2
      
      if(comparison=='Baseline between arms'){
        cat(sum(sign1,na.rm=TRUE),'studies with lower variance in experimental group\n')
        cat(sum(sign2,na.rm=TRUE),'studies with greater variance in experimental group\n')
        cat(sum(!is.na(x))-sum(sign1,na.rm=TRUE)-sum(sign2,na.rm=TRUE),'studies with no significant differences in variance\n')
      }else{
        cat(sum(sign1,na.rm=TRUE),'studies with lower variance at the end of the study\n')
        cat(sum(sign2,na.rm=TRUE),'studies with greater variance at the end of the study\n')
        cat(sum(!is.na(y))-sum(sign1,na.rm=TRUE)-sum(sign2,na.rm=TRUE),'studies with no significant differences in variance\n')
      }
      
      ##-- Point colors
      co1 <- rgb(1,0,0,0.5)
      co2 <- rgb(0,0,0,0.5)
      co <- ifelse(sign,co1,co2)
      pc <- 19
      legend('topright',c('Different variances','Non different variances'),
             pch=19,text.font = 2,pt.cex=1.3,
             co=c(co1,co2))
    }
    if(subgroup=='effect'){
      co <- with(data,ifelse(pvalue>=0.05 | significant=='No',rgb(0,0,0,0.8),
                             ifelse(pvalue>0.001,rgb(227,66,52,255/2,maxColorValue = 255),
                                    rgb(128,0,0,255/2,maxColorValue = 255))))#'#E34234','#800000'))))
      pc <- with(data,ifelse(significant=='No',1,19))
      co[is.na(co)] <- ifelse(data$significant[is.na(co)]=='No',rgb(0,0,0,0.8),rgb(128,0,0,255/2,maxColorValue = 255))
      legend('topright',c('p > 0.05','0.001 < p < 0.05','p < 0.001'),
             pch=c(1,19,19,19),pt.lwd = 2,text.font = 2,pt.cex=1.3,
             co=c(rgb(0,0,0,0.8),rgb(227,66,52,255/2,maxColorValue = 255),rgb(128,0,0,255/2,maxColorValue = 255)))
    }
  }

  ##-- plot points
  points(x,y,pch=pc,col=co,lwd=2,cex=1.1)

  ##-- Labels for both sides of x-axis
  if(zoom){
    mtext(lab.11,1,adj=1,at=log(xmax),line=2.7,cex=.8,font=2,las=0)
    mtext(lab.21,1,adj=1,at=log(xmax),line=3.5,cex=.8,font=2,las=0)
    mtext(lab.12,1,adj=0,at=log(1/xmax),line=2.7,cex=.8,font=2,las=0)
    mtext(lab.22,1,adj=0,at=log(1/xmax),line=3.5,cex=.8,font=2,las=0) 
  }else{
    mtext(lab.11,1,adj=1,at=log(100),line=2.7,cex=.8,font=2,las=0)
    mtext(lab.21,1,adj=1,at=log(100),line=3.5,cex=.8,font=2,las=0)
    mtext(lab.12,1,adj=0,at=log(0.01),line=2.7,cex=.8,font=2,las=0)
    mtext(lab.22,1,adj=0,at=log(0.01),line=3.5,cex=.8,font=2,las=0)
  }
}

##########################################################################
# limit_absolute_ratio. Given a limit for the ratio between th signal (y) and the noise (se), 
# this function selects the subset of studies with lower ratio than this limit and
# fit the rma model with this data. It returns the difference between the heterogeneity
# of the model (tau) and a desired target heterogeneity (quantile_tau). Thid function is used
# to find the optimal limit that minimizes this difference.
#
# limit: ratio between signal (y) and noise (se).
# data: complete dataset
# y: variance ratio
# quantile_tau: target tau; target heterogeneity.
##########################################################################
limit_absolute_ratio <- function(limit,data,y,se,comparison,quantile_tau){
  sel <- abs(y/se)<limit
  data_reduced <- data[sel & !is.na(se),]
  data_reduced$y <- y[sel & !is.na(se)]
  data_reduced$se<- se[sel & !is.na(se)]
  
  if(comparison=='Baseline'){
    mod <- rma(y,sei=se,data=data_reduced,method='REML')
  }else if(comparison=='BA'){
    mod <- rma(y,sei=se,data=data_reduced,mods=~yBaselineRatio,method='REML')
  }else if(comparison=='OT'){
    mod <- rma(y,sei=se,data=data_reduced,mods=~yOverTimeRatioC,method='REML')
  }
  return(sqrt(mod$tau2)-quantile_tau)
}


####################################################################
# gg_color_hue. Function to defin colors in a hue scale
#
# n: number of colors to define
####################################################################
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

####################################################################
# QLIM. Estimates the limits of the paired test. This function is optimized to finf
# a specific
#
# Qx: square sum (x-\bar{x})^2 at baseline
# Qx: square sum (y-\bar{y})^2 at the end of the study
# Qxy: (n-1)*Covariance
# n: sample size
####################################################################
QLIM <- function(Qx,Qy,Qxy,n){
  num <- (Qx-Qy)*sqrt(n-2)
  aux <- Qx*Qy-Qxy^2
  den <- 2*sqrt(max(-aux,aux))
  Qest <- num/den
  return(Qest - qt(c(0.975,0.025),n-2))
}
QLIM1 <- function(Qx,Qy,Qxy,n) QLIM(Qx,Qy,Qxy,n)[1]
QLIM2 <- function(Qx,Qy,Qxy,n) QLIM(Qx,Qy,Qxy,n)[2]