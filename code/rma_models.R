sink("../results_tables/MA_models.txt",split = TRUE)

data <- datos1

############################################################
# Model - Reference model
# log(Vbt/Vbc)~1
############################################################
rma.unadjB <- rma(yBaselineRatio,sei=seBaselineRatio,data=data)                                      # Unadjusted

cat('\n*************************************************************************\n')
cat('Model 1 (Reference): log(Vbt/Vbc)~1 -------------------------------------\n')
cat('*************************************************************************\n')
print(rma.unadjB)
cat('-------------------------------------------------------------------------\n\n')

############################################################
# Model - Comparison between arms
# 
############################################################
rma.unadj <- rma(yBetweenArmsRatio,sei=seBetweenArmsRatio,data=data)                                                 # Unadjusted
rma.adj <-   rma(yBetweenArmsRatio,sei=seBetweenArmsRatio,data=data,mods=~yBaselineRatio,method='REML')              # Adjusted by baseline
rma.adjoff <- rma(yBetweenArmsRatio-yBaselineRatio,sei=seBetweenArmsRatio,data=data,method='REML')                   # Baseline as offset

##-- Output
cat('\n*************************************************************************\n')
cat('Model 2 (Between arms unadjusted): log(Vot/Voc)~1 -----------------------\n')
cat('*************************************************************************\n')
print(rma.unadj)
cat('-------------------------------------------------------------------------\n\n')

cat('\n*************************************************************************\n')
cat('Model 3 (Between arms adjusted): log(Vot/Voc)~ B*log(Vbt/Vbc)------------\n')
cat('*************************************************************************\n')
print(rma.adj)
cat('-------------------------------------------------------------------------\n\n')

cat('\n*************************************************************************\n')
cat('Model 4 (Between arms with offset): log(Vot/Voc)~ 1*log(Vbt/Vbc)---------\n')
cat('*************************************************************************\n')
print(rma.adjoff)
cat('-------------------------------------------------------------------------\n\n')

############################################################
# Model - Comparison over time
# 
############################################################
data <- datos1[!is.na(datos1$seOverTimeRatioT),]

rma.unadj2 <- rma(yOverTimeRatioT,sei=seOverTimeRatioT,data=data)                                                  # Unadjusted
rma.adj2 <-   rma(yOverTimeRatioT,sei=seOverTimeRatioT,data=data,mods=~yOverTimeRatioC,method='REML')              # Adjusted by baseline
rma.adjoff2 <- rma(yOverTimeRatioT-yOverTimeRatioC,sei=seOverTimeRatioT,data=data,method='REML')                   # Baseline as offset

##-- Print
cat('\n*************************************************************************\n')
cat('Model 5 (Over time unadjusted): log(Voc/Vbt)~1 --------------------------\n')
cat('*************************************************************************\n')
print(rma.unadj2)
cat('-------------------------------------------------------------------------\n\n')

cat('\n*************************************************************************\n')
cat('Model 6 (Over time adjusted): log(Vot/Vbt)~ B*log(Voc/Vbc) --------------\n')
cat('*************************************************************************\n')
print(rma.adj2)
cat('-------------------------------------------------------------------------\n\n')

cat('\n*************************************************************************\n')
cat('Model 7 (Over time with offset): log(Vot/Vbt)~ 1*log(Voc/Vbc) -----------\n')
cat('*************************************************************************\n')
print(rma.adjoff2)
cat('-------------------------------------------------------------------------\n\n')

sink()
