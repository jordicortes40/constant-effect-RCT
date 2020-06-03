sink('../results_tables/SA_I_reduced_models.txt',split = FALSE)

############################################################
# Reduced data --> Model - Reference model
# log(Vbt/Vbc)~1
############################################################
rma.red.unadjB <- rma(yBaselineRatio,sei=seBaselineRatio,data=dataB)                                      # Unadjusted

cat('\n******************************************************************************************\n')
cat('Reduced data --> Model 1 (Reference): log(Vbt/Vbc)~1 -------------------------------------\n')
cat('******************************************************************************************\n')
print(rma.red.unadjB)
cat('-------------------------------------------------------------------------------------------\n\n')

############################################################
# Reduced data --> Model - Comparison between arms
# 
############################################################
rma.red.unadj <- rma(yBetweenArmsRatio,sei=seBetweenArmsRatio,data=dataBetArm)                                                 # Unadjusted
rma.red.adj <-   rma(yBetweenArmsRatio,sei=seBetweenArmsRatio,data=dataBetArm,mods=~yBaselineRatio,method='REML')              # Adjusted by baseline
rma.red.adjoff <- rma(yBetweenArmsRatio-yBaselineRatio,sei=seBetweenArmsRatio,data=dataBetArm,method='REML')                   # Baseline as offset

##-- Output
cat('\n******************************************************************************************\n')
cat('Reduced data --> Model 2 (Between arms unadjusted): log(Vot/Voc)~1 -----------------------\n')
cat('******************************************************************************************\n')
print(rma.red.unadj)
cat('-------------------------------------------------------------------------------------------\n\n')

cat('\n******************************************************************************************\n')
cat('Reduced data --> Model 3 (Between arms adjusted): log(Vot/Voc)~ B?log(Vbt/Vbc)------------\n')
cat('******************************************************************************************\n')
print(rma.red.adj)
cat('-------------------------------------------------------------------------------------------\n\n')

cat('\n******************************************************************************************\n')
cat('Reduced data --> Model 4 (Between arms with offset): log(Vot/Voc)~ 1?log(Vbt/Vbc)---------\n')
cat('******************************************************************************************\n')
print(rma.red.adjoff)
cat('-------------------------------------------------------------------------------------------\n\n')

############################################################
# Reduced data --> Model - Comparison over time
# 
############################################################
rma.red.unadj2 <- rma(yOverTimeRatioT,sei=seOverTimeRatioT,data=dataOverTime)                                                  # Unadjusted
rma.red.adj2 <-   rma(yOverTimeRatioT,sei=seOverTimeRatioT,data=dataOverTime,mods=~yOverTimeRatioC,method='REML')              # Adjusted by baseline
rma.red.adj2off <- rma(yOverTimeRatioT-yOverTimeRatioC,sei=seOverTimeRatioT,data=dataOverTime,method='REML')                   # Baseline as offset

##-- Output
cat('\n******************************************************************************************\n')
cat('Reduced data --> Model 5 (Over time unadjusted): log(Voc/Vbt)~1 --------------------------\n')
cat('******************************************************************************************\n')
print(rma.red.unadj2)
cat('-------------------------------------------------------------------------------------------\n\n')

cat('\n******************************************************************************************\n')
cat('Reduced data --> Model 6 (Over time adjusted): log(Vot/Vbt)~ B?log(Voc/Vbc) --------------\n')
cat('******************************************************************************************\n')
print(rma.red.adj2)
cat('-------------------------------------------------------------------------------------------\n\n')

cat('\n******************************************************************************************\n')
cat('Reduced data --> Model 7 (Over time with offset): log(Vot/Vbt)~ 1?log(Voc/Vbc) -----------\n')
cat('******************************************************************************************\n')
print(rma.red.adj2off)
cat('-------------------------------------------------------------------------------------------\n\n')

sink()