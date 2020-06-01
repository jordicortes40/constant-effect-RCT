rm(list=ls())

#-----------------------------------------------------------------
#
# Read data
#
#-----------------------------------------------------------------
URL <- 'https://raw.githubusercontent.com/jordicortes40/constant-effect-RCT/master/'
source(paste0(URL,'code/read_data.R')) 
source(paste0(URL,'code/pvalue_distribution_mixture_functions.R'))      # Functions for the likelihood

#-----------------------------------------------------------------
#
# Between arms
#
#-----------------------------------------------------------------
d <- read.table(paste0(URL,'results_tables/SA_III_pvalues_BA.txt'),
                header=TRUE,sep='\t')

##############################################################
# Descriptive
##############################################################
gg1 <- ggplot(d,aes(x=p_values,y=stat(density))) + geom_histogram(col='white',bins=10,breaks=seq(0,1,0.1),fill=rgb(0.1992188,0.1992188,0.6953125,maxColorValue = 1)) + 
  scale_x_continuous() + xlab('p-values') + ylab('n') + labs(title='Histogram') +
  geom_hline(yintercept = 1,col='grey10',linetype=2,size=1.2) +
  theme(axis.title = element_text(face='bold',size = 13),
        axis.text = element_text(face='bold',size = 10),
        title = element_text(face='bold',size = 15))
gg2 <- ggplot(d,aes(sample=p_values)) + stat_qq(distribution = qunif,size=2,alpha=0.5,color=rgb(0.1992188,0.1992188,0.6953125,maxColorValue = 1)) + 
  stat_qq_line(distribution = qunif) +
  scale_x_continuous() + xlab('Theoretical Uniform') + ylab('Sample') + labs(title='QQplot') +
  theme(axis.title = element_text(face='bold',size = 13),
        axis.text = element_text(face='bold',size = 10),
        title = element_text(face='bold',size = 15))
ggarrange(gg1,gg2,nrow=1)
ggsave(filename = '../results_figures/SA_IV_pvalue_descriptive_BA.jpg',width=8,height=4)

##############################################################
# Models estimation
##############################################################
# If some pvalues are almost 0 or 1, some bounds could be added --> pmin(pmax(datos,10^-5),1-10^-5)
var = list(datos=d[,1])     

##-- Uniform + 2 Triangulars
P0 = c(0.4, 0.3, 0.05, 0.95)
ans1 = auglag(par=P0, fn=fn1, heq=heq1,hin=hin1, hin.jac=hin.jac1, control.outer=list(trace=FALSE), var=var)

##-- Uniform + 2 Exponential
P0 = c(0.8, 0.1, 0.1, 0.1)
ans2 = auglag(par=P0, fn=fn2, heq=heq2,hin=hin2, hin.jac=hin.jac2, control.outer=list(trace=FALSE), var=var)

##-- Uniform + Beta 
P0 = c(0.4, 1, 1)
ans3 = auglag(par=P0, fn=fn3, heq=heq3,hin=hin3, hin.jac=hin.jac3, control.outer=list(trace=FALSE), var=var)

##-- Uniform + 2 Betas
P0 = c(0.4, 0.3, 1, 1, 1, 1)
ans4 = auglag(par=P0, fn=fn4, heq=heq4,hin=hin4, hin.jac=hin.jac4, control.outer=list(trace=FALSE), var=var)


##############################################################
# Check distributions using a QQplot
##############################################################
x <- seq(0,1,0.001)
fit <- data.frame(x=x,
                  fit1=F1(x,ans1$par),fit2=F2(x,ans2$par),
                  fit3=F3(x,ans3$par),fit4=F4(x,ans4$par))
common_theme <- theme(axis.text = element_text(face='bold',size=10),
                      axis.title = element_text(face='bold',size=13) )
gg1 <- ggplot(d, aes(x=p_values)) + stat_ecdf(geom = "step") + 
  geom_line(data = fit,mapping=aes(x=x,y=fit1),size=1,col='darkblue',alpha=0.5) +
  ggtitle('Uniform + 2 Triangular') + xlab('p values') + ylab(expression(F[n](x))) +
  common_theme

gg2 <- ggplot(d, aes(x=p_values)) + stat_ecdf(geom = "step") + 
  geom_line(data = fit,mapping=aes(x=x,y=fit2),size=1,col='darkblue',alpha=0.5) +
  ggtitle('Uniform + 2 Exponential') + xlab('p values') + ylab(expression(F[n](x))) +
  common_theme

gg3 <- ggplot(d, aes(x=p_values)) + stat_ecdf(geom = "step") + 
  geom_line(data = fit,mapping=aes(x=x,y=fit3),size=1,col='darkblue',alpha=0.5) +
  ggtitle('Uniform + Beta') + xlab('p values') + ylab(expression(F[n](x))) +
  common_theme

gg4 <- ggplot(d, aes(x=p_values)) + stat_ecdf(geom = "step") + 
  geom_line(data = fit,mapping=aes(x=x,y=fit4),size=1,col='darkblue',alpha=0.5) +
  ggtitle('Uniform + 2 Beta') + xlab('p values') + ylab(expression(F[n](x))) +
  common_theme

ggarrange(gg1,gg2,gg3,gg4,nrow=2,ncol=2)
ggsave(filename = '../results_figures/SA_IV_qqplots_BA.jpg',width=8,height=8)

##############################################################
# Check distributions using:
# 1. Histogram with mixture distribution density overlapped
# 2. Statistics, such as AIC and Komogorov- Smirnov
# results also provide the proportion (pu) of pvalues comming 
# from the uniform distribution
##############################################################
par(las=1,mfrow=c(2,2))
RES <- rbind(results(d[,1],ans1,nparam=4,F1,f1,'Uniform + 2 Triangular'),
             results(d[,1],ans2,nparam=4,F2,f2,'Uniform + 2 Exponentials'),
             results(d[,1],ans3,nparam=3,F3,f3,'Uniform + Beta'),
             results(d[,1],ans4,nparam=6,F4,f4,'Uniform + 2 Betas'))
colnames(RES) <- c('nparam','LL','AIC','KS1','KS2','pu','95%LL(pu)','95%UL(pu)','max_f(x)','iterations')
rownames(RES) <- c('Triangulars','Exponentials','Beta1','Beta2')

RES
write.table(x = as.data.frame(RES),
            file='../results_tables/SA_IV_results_BA.txt',
            row.names = TRUE,col.names = TRUE,sep='\t',quote = FALSE)


## Points to assess
x <- seq(0,1,0.001)

##-- Uniform + Beta
d1 <- data.frame(x=x,y=pmin(2,f1(x,ans1$par)))
gg1 <- ggplot(d,aes(x=p_values)) + geom_histogram(aes(y=..density..),col='white',bins=10,breaks=seq(-0.1,1,0.1),fill=rgb(0.1992188,0.1992188,0.6953125,maxColorValue = 1)) + 
  scale_x_continuous(limits=c(0,1)) + scale_y_continuous(limits=c(0,2)) + 
  xlab('pvalues') + ylab('n') + labs(title='Uniform + 2 Triangulars') +
  geom_hline(yintercept = 1,col='grey10',linetype=2,size=1.2) +
  geom_ribbon(data=d1,mapping=aes(x=x,ymax=y),ymin=0,fill='red',alpha=0.3) +
  theme(axis.title = element_text(face='bold',size = 13),
        axis.text = element_text(face='bold',size = 10),
        title = element_text(face='bold',size = 15))

##-- Uniform + 2 Betas
d2 <- data.frame(x=x,y=pmin(2,f2(x,ans2$par)))
gg2 <- ggplot(d,aes(x=p_values)) + geom_histogram(aes(y=..density..),col='white',bins=10,breaks=seq(-0.1,1,0.1),fill=rgb(0.1992188,0.1992188,0.6953125,maxColorValue = 1)) + 
  scale_x_continuous(limits=c(0,1)) + scale_y_continuous(limits=c(0,2)) + 
  xlab('pvalues') + ylab('n') + labs(title='Uniform + 2 Exponentials') +
  geom_hline(yintercept = 1,col='grey10',linetype=2,size=1.2) +
  geom_ribbon(data=d2,mapping=aes(x=x,ymax=y),ymin=0,fill='red',alpha=0.3) +
  theme(axis.title = element_text(face='bold',size = 13),
        axis.text = element_text(face='bold',size = 10),
        title = element_text(face='bold',size = 15))

##-- Uniform + 2 Triangular
d3 <- data.frame(x=x,y=pmin(2,f3(x,ans3$par)))
gg3 <- ggplot(d,aes(x=p_values)) + geom_histogram(aes(y=..density..),col='white',bins=10,breaks=seq(-0.1,1,0.1),fill=rgb(0.1992188,0.1992188,0.6953125,maxColorValue = 1)) + 
  scale_x_continuous(limits=c(0,1)) + scale_y_continuous(limits=c(0,2)) + 
  xlab('pvalues') + ylab('n') + labs(title='Uniform + Beta') +
  geom_hline(yintercept = 1,col='grey10',linetype=2,size=1.2) +
  geom_ribbon(data=d3,mapping=aes(x=x,ymax=y),ymin=0,fill='red',alpha=0.3) +
  theme(axis.title = element_text(face='bold',size = 13),
        axis.text = element_text(face='bold',size = 10),
        title = element_text(face='bold',size = 15))

##-- Uniform + 2 Exponentials
d4 <- data.frame(x=x,y=pmin(2,f4(x,ans4$par)))
gg4 <- ggplot(d,aes(x=p_values)) + geom_histogram(aes(y=..density..),col='white',bins=10,breaks=seq(-0.1,1,0.1),fill=rgb(0.1992188,0.1992188,0.6953125,maxColorValue = 1)) + 
  scale_x_continuous(limits=c(0,1)) + scale_y_continuous(limits=c(0,2)) + 
  xlab('pvalues') + ylab('n') + labs(title='Uniform + 2 Betas') +
  geom_hline(yintercept = 1,col='grey10',linetype=2,size=1.2) +
  geom_ribbon(data=d4,mapping=aes(x=x,ymax=y),ymin=0,fill='red',alpha=0.3) +
  theme(axis.title = element_text(face='bold',size = 13),
        axis.text = element_text(face='bold',size = 10),
        title = element_text(face='bold',size = 15))

ggarrange(gg1,gg2,gg3,gg4,nrow=2,ncol=2)
ggsave(filename = '../results_figures/SA_IV_histograms_BA.jpg',width=8,height=8)


#-----------------------------------------------------------------
#
# Data for summary table
#
#-----------------------------------------------------------------
##-- Between arms (1 beta)
p_equal_BA <- ans3$par[1]
p_greater_BA <- (1-p_equal_BA)*(ans3$par[2]/(ans3$par[2]+ans3$par[3])) # heuristic
p_lower_BA <- (1-p_equal_BA)*(ans3$par[3]/(ans3$par[2]+ans3$par[3])) # heuristic
SAIV_greater_BA <- round(nrow(d)*p_greater_BA)
SAIV_lower_BA <- round(nrow(datos1)*p_lower_BA)
SAIV_equal_BA <- round(nrow(datos1)*p_equal_BA)
SAIV_greater_BA; SAIV_lower_BA; SAIV_equal_BA


#-----------------------------------------------------------------
#
# Over time
#
#-----------------------------------------------------------------
d <- read.table(paste0(URL,'results_tables/SA_III_pvalues_OT.txt'),
                header=TRUE,sep='\t')

##############################################################
# Descriptive
##############################################################
gg1 <- ggplot(d,aes(x=p_values,y=stat(density))) + geom_histogram(col='white',bins=10,breaks=seq(0,1,0.1),fill=rgb(0.1992188,0.1992188,0.6953125,maxColorValue = 1)) + 
  scale_x_continuous() + xlab('p-values') + ylab('n') + labs(title='Histogram') +
  geom_hline(yintercept = 1,col='grey10',linetype=2,size=1.2) +
  theme(axis.title = element_text(face='bold',size = 13),
        axis.text = element_text(face='bold',size = 10),
        title = element_text(face='bold',size = 15))
gg2 <- ggplot(d,aes(sample=p_values)) + stat_qq(distribution = qunif,size=2,alpha=0.5,color=rgb(0.1992188,0.1992188,0.6953125,maxColorValue = 1)) + 
  stat_qq_line(distribution = qunif) +
  scale_x_continuous() + xlab('Theoretical Uniform') + ylab('Sample') + labs(title='QQplot') +
  theme(axis.title = element_text(face='bold',size = 13),
        axis.text = element_text(face='bold',size = 10),
        title = element_text(face='bold',size = 15))
ggarrange(gg1,gg2,nrow=1)
ggsave(filename = '../results_figures/SA_IV_pvalue_descriptive_OT.jpg',width=8,height=4)

##############################################################
# Models estimation
##############################################################
# If some pvalues are almost 0 or 1, some bounds could be added --> pmin(pmax(datos,10^-5),1-10^-5)
var = list(datos=d[,1])     

##-- Uniform + 2 Triangulars
P0 = c(0.4, 0.3, 0.05, 0.95)
ans1 = auglag(par=P0, fn=fn1, heq=heq1,hin=hin1, hin.jac=hin.jac1, control.outer=list(trace=FALSE), var=var)

##-- Uniform + 2 Exponential
P0 = c(0.8, 0.1, 0.1, 0.1)
ans2 = auglag(par=P0, fn=fn2, heq=heq2,hin=hin2, hin.jac=hin.jac2, control.outer=list(trace=FALSE), var=var)

##-- Uniform + Beta 
P0 = c(0.4, 1, 1)
ans3 = auglag(par=P0, fn=fn3, heq=heq3,hin=hin3, hin.jac=hin.jac3, control.outer=list(trace=FALSE), var=var)

##-- Uniform + 2 Betas
P0 = c(0.4, 0.3, 1, 1, 1, 1)
ans4 = auglag(par=P0, fn=fn4, heq=heq4,hin=hin4, hin.jac=hin.jac4, control.outer=list(trace=FALSE), var=var)


##############################################################
# Check distributions using a QQplot
##############################################################
x <- seq(0,1,0.001)
fit <- data.frame(x=x,
                  fit1=F1(x,ans1$par),fit2=F2(x,ans2$par),
                  fit3=F3(x,ans3$par),fit4=F4(x,ans4$par))
common_theme <- theme(axis.text = element_text(face='bold',size=10),
                      axis.title = element_text(face='bold',size=13) )
gg1 <- ggplot(d, aes(x=p_values)) + stat_ecdf(geom = "step") + 
  geom_line(data = fit,mapping=aes(x=x,y=fit1),size=1,col='darkblue',alpha=0.5) +
  ggtitle('Uniform + 2 Triangular') + xlab('p values') + ylab(expression(F[n](x))) +
  common_theme

gg2 <- ggplot(d, aes(x=p_values)) + stat_ecdf(geom = "step") + 
  geom_line(data = fit,mapping=aes(x=x,y=fit2),size=1,col='darkblue',alpha=0.5) +
  ggtitle('Uniform + 2 Exponential') + xlab('p values') + ylab(expression(F[n](x))) +
  common_theme

gg3 <- ggplot(d, aes(x=p_values)) + stat_ecdf(geom = "step") + 
  geom_line(data = fit,mapping=aes(x=x,y=fit3),size=1,col='darkblue',alpha=0.5) +
  ggtitle('Uniform + Beta') + xlab('p values') + ylab(expression(F[n](x))) +
  common_theme

gg4 <- ggplot(d, aes(x=p_values)) + stat_ecdf(geom = "step") + 
  geom_line(data = fit,mapping=aes(x=x,y=fit4),size=1,col='darkblue',alpha=0.5) +
  ggtitle('Uniform + 2 Beta') + xlab('p values') + ylab(expression(F[n](x))) +
  common_theme

ggarrange(gg1,gg2,gg3,gg4,nrow=2,ncol=2)
ggsave(filename = '../results_figures/SA_IV_qqplots_OT.jpg',width=8,height=8)

##############################################################
# Check distributions using:
# 1. Histogram with mixture distribution density overlapped
# 2. Statistics, such as AIC and Komogorov- Smirnov
# results also provide the proportion (pu) of pvalues comming 
# from the uniform distribution
##############################################################
par(las=1,mfrow=c(2,2))
RES <- rbind(results(d[,1],ans1,nparam=4,F1,f1,'Uniform + 2 Triangular'),
             results(d[,1],ans2,nparam=4,F2,f2,'Uniform + 2 Exponentials'),
             results(d[,1],ans3,nparam=3,F3,f3,'Uniform + Beta'),
             results(d[,1],ans4,nparam=6,F4,f4,'Uniform + 2 Betas'))
colnames(RES) <- c('nparam','LL','AIC','KS1','KS2','pu','95%LL(pu)','95%UL(pu)','max_f(x)','iterations')
rownames(RES) <- c('Triangulars','Exponentials','Beta1','Beta2')

RES
write.table(x = as.data.frame(RES),
            file='../results_tables/SA_IV_results_OT.txt',
            row.names = TRUE,col.names = TRUE,sep='\t',quote = FALSE)


## Points to assess
x <- seq(0,1,0.001)

##-- Uniform + Beta
d1 <- data.frame(x=x,y=pmin(2,f1(x,ans1$par)))
gg1 <- ggplot(d,aes(x=p_values)) + geom_histogram(aes(y=..density..),col='white',bins=10,breaks=seq(-0.1,1,0.1),fill=rgb(0.1992188,0.1992188,0.6953125,maxColorValue = 1)) + 
  scale_x_continuous(limits=c(0,1)) + scale_y_continuous(limits=c(0,2)) + 
  xlab('pvalues') + ylab('n') + labs(title='Uniform + 2 Triangulars') +
  geom_hline(yintercept = 1,col='grey10',linetype=2,size=1.2) +
  geom_ribbon(data=d1,mapping=aes(x=x,ymax=y),ymin=0,fill='red',alpha=0.3) +
  theme(axis.title = element_text(face='bold',size = 13),
        axis.text = element_text(face='bold',size = 10),
        title = element_text(face='bold',size = 15))

##-- Uniform + 2 Betas
d2 <- data.frame(x=x,y=pmin(2,f2(x,ans2$par)))
gg2 <- ggplot(d,aes(x=p_values)) + geom_histogram(aes(y=..density..),col='white',bins=10,breaks=seq(-0.1,1,0.1),fill=rgb(0.1992188,0.1992188,0.6953125,maxColorValue = 1)) + 
  scale_x_continuous(limits=c(0,1)) + scale_y_continuous(limits=c(0,2)) + 
  xlab('pvalues') + ylab('n') + labs(title='Uniform + 2 Exponentials') +
  geom_hline(yintercept = 1,col='grey10',linetype=2,size=1.2) +
  geom_ribbon(data=d2,mapping=aes(x=x,ymax=y),ymin=0,fill='red',alpha=0.3) +
  theme(axis.title = element_text(face='bold',size = 13),
        axis.text = element_text(face='bold',size = 10),
        title = element_text(face='bold',size = 15))

##-- Uniform + 2 Triangular
d3 <- data.frame(x=x,y=pmin(2,f3(x,ans3$par)))
gg3 <- ggplot(d,aes(x=p_values)) + geom_histogram(aes(y=..density..),col='white',bins=10,breaks=seq(-0.1,1,0.1),fill=rgb(0.1992188,0.1992188,0.6953125,maxColorValue = 1)) + 
  scale_x_continuous(limits=c(0,1)) + scale_y_continuous(limits=c(0,2)) + 
  xlab('pvalues') + ylab('n') + labs(title='Uniform + Beta') +
  geom_hline(yintercept = 1,col='grey10',linetype=2,size=1.2) +
  geom_ribbon(data=d3,mapping=aes(x=x,ymax=y),ymin=0,fill='red',alpha=0.3) +
  theme(axis.title = element_text(face='bold',size = 13),
        axis.text = element_text(face='bold',size = 10),
        title = element_text(face='bold',size = 15))

##-- Uniform + 2 Exponentials
d4 <- data.frame(x=x,y=pmin(2,f4(x,ans4$par)))
gg4 <- ggplot(d,aes(x=p_values)) + geom_histogram(aes(y=..density..),col='white',bins=10,breaks=seq(-0.1,1,0.1),fill=rgb(0.1992188,0.1992188,0.6953125,maxColorValue = 1)) + 
  scale_x_continuous(limits=c(0,1)) + scale_y_continuous(limits=c(0,2)) + 
  xlab('pvalues') + ylab('n') + labs(title='Uniform + 2 Betas') +
  geom_hline(yintercept = 1,col='grey10',linetype=2,size=1.2) +
  geom_ribbon(data=d4,mapping=aes(x=x,ymax=y),ymin=0,fill='red',alpha=0.3) +
  theme(axis.title = element_text(face='bold',size = 13),
        axis.text = element_text(face='bold',size = 10),
        title = element_text(face='bold',size = 15))

ggarrange(gg1,gg2,gg3,gg4,nrow=2,ncol=2)
ggsave(filename = '../results_figures/SA_IV_histograms_OT.jpg',width=8,height=8)


#-----------------------------------------------------------------
#
# Data for summary table
#
#-----------------------------------------------------------------
##-- Over time (1 beta)
p_equal_OT <- ans3$par[1]
p_greater_OT <- (1-p_equal_OT)*(ans3$par[2]/(ans3$par[2]+ans3$par[3])) # heuristic
p_lower_OT <- (1-p_equal_OT)*(ans3$par[3]/(ans3$par[2]+ans3$par[3])) # heuristic
SAIV_greater_OT <- round(nrow(d)*p_greater_OT)
SAIV_lower_OT <- round(nrow(d)*p_lower_OT)
SAIV_equal_OT <- round(nrow(d)*p_equal_OT)
SAIV_greater_OT; SAIV_lower_OT; SAIV_equal_OT