#-----------------------------------------------------------
#
# Subgroups analysis
#
#-----------------------------------------------------------
############################################################
# Estimate coefficients
############################################################

##-- Subgroup names
var.sg <- c(NA,'significant','Intervention_type','Outcome_type','Condition_type','Measurement_type')

##-- Objects for storing information
m <- matrix(nrow=11,ncol=4,dimnames = list(1:11,c('Ratio','LL_CI95','UL_CI95','n')))
M <- list(BetweenArms=m,Overtime=m,Full=m)

##-- Estimate data for forest-plots
for(v in var.sg){
  d0 <- datos1
  
  if(!is.na(v)){
    d1 <- d0[which(d0[,v]==levels(d0[,v])[1]),]
    d2 <- d0[which(d0[,v]==levels(d0[,v])[2]),]
  }else{
    d1 <- d0
  }
  
  
  rma.Arms1 <- rma(yBetweenArmsRatio,sei=seBetweenArmsRatio,data=d1,mods=~yBaselineRatio,method='REML') ##-- Between arms (adjusted)
  rma.Time1 <- rma(yOverTimeRatioT,sei=seOverTimeRatioT,data=d1,mods=~yOverTimeRatioC,method='REML')    ##-- Over Time (adjusted)
  rma.Full1 <- rma(yBetweenArmsRatio-yBaselineRatio,sei=seBetweenArmsRatio,data=d1,method='REML')       ##-- Full (adjusted)
  
  ##-- Store data
  row <- ifelse(is.na(v),1,2*which(var.sg==v)-2)
  
  M[[1]][row,] <- c(exp(c(coef(rma.Arms1)[1],rma.Arms1$ci.lb[1],rma.Arms1$ci.ub[1])),nrow(d1))
  M[[2]][row,] <- c(exp(c(coef(rma.Time1)[1],rma.Time1$ci.lb[1],rma.Time1$ci.ub[1])),nrow(d1))
  M[[3]][row,] <- c(exp(c(coef(rma.Full1)[1],rma.Full1$ci.lb[1],rma.Full1$ci.ub[1])),nrow(d1))
  
  if(!is.na(v)){for(i in 1:3) rownames(M[[i]])[row] <- levels(d1[,v])[1]}else{for(i in 1:3) rownames(M[[i]])[row] <- ''}
  
  if(!is.na(v)){
    rma.Arms2 <- rma(yBetweenArmsRatio,sei=seBetweenArmsRatio,data=d2,mods=~yBaselineRatio,method='REML') ##-- Between arms (adjusted)
    rma.Time2 <- rma(yOverTimeRatioT,sei=seOverTimeRatioT,data=d2,mods=~yOverTimeRatioC,method='REML')    ##-- Over Time (adjusted)
    rma.Full2 <- rma(yBetweenArmsRatio-yBaselineRatio,sei=seBetweenArmsRatio,data=d2,method='REML')       ##-- Full (adjusted)
    
    M[[1]][row+1,] <- c(exp(c(coef(rma.Arms2)[1],rma.Arms2$ci.lb[1],rma.Arms2$ci.ub[1])),nrow(d2))
    M[[2]][row+1,] <- c(exp(c(coef(rma.Time2)[1],rma.Time2$ci.lb[1],rma.Time2$ci.ub[1])),nrow(d2))
    M[[3]][row+1,] <- c(exp(c(coef(rma.Full2)[1],rma.Full2$ci.lb[1],rma.Full2$ci.ub[1])),nrow(d2))
    
    if(!is.na(v)){for(i in 1:3) rownames(M[[i]])[row+1] <- levels(d1[,v])[2]}else{for(i in 1:3) rownames(M[[i]])[row+1] <- ''}
  }
}

############################################################
# Forest-plot
############################################################
##-- Labels for forest-plot
# Levels
lab1 <- with(datos1,c(NA,NA,NA,levels(significant),NA,levels(Intervention_type),
                      NA,levels(Outcome_type),NA,levels(Condition_type)[1:2],
                      NA,levels(Measurement_type)))
# Variables
lab2 <- c('Global',NA,'Significant',NA,NA,'Intervention type',NA,NA,
          'Outcome type',NA,NA,'Condition type*',NA,NA,
          'Measurement type',NA,NA)

##-- Between arms
png('../results_figures/MA_subgroup_analysis_BA.png',width=960,height = 670,res=72)
myForest(M=M[[1]],xl=bquote(bold(frac(S[OT]^2,S[OC]^2))),lab1,lab2,
         tit='Between Arms - Subgroups',
         laxis1=c('in reference arm','in experimental arm'))
dev.off()

##-- Over-time
png('../results_figures/MA_subgroup_analysis_OT.png',width=960,height = 670,res=72)
myForest(M=M[[2]],xl=bquote(bold(frac(S[OT]^2,S[BT]^2))),lab1,lab2,
         tit='Over time - Subgroups',
         laxis1=c('at baseline','at the end of study'))
dev.off()

##-- Between arms over-time
png('../results_figures/MA_subgroup_analysis_BA_OT.png',width=960,height = 670,res=72)
myForest(M=M[[3]],xl=bquote(bold(frac(S[OT]^2/S[BT]^2,S[OC]^2/S[BC]^2))),lab1,lab2,
         tit='Change Over Time Between Arms - Subgroups',
         laxis1=c('in reference arm','in experimental arm'))
dev.off()

##-- Mean by subgroups
sink('../results_tables/MA_subgroup_means.txt',split = FALSE)
with(datos1,tapply(yBetweenArmsRatio,significant,mean,na.rm=TRUE))
with(datos1,tapply(yBetweenArmsRatio,Intervention_type,mean,na.rm=TRUE))
with(datos1,tapply(yBetweenArmsRatio,Outcome_type,mean,na.rm=TRUE))
with(datos1,tapply(yBetweenArmsRatio,Condition_type,mean,na.rm=TRUE))
with(datos1,tapply(yBetweenArmsRatio,Measurement_type,mean,na.rm=TRUE))
sink()