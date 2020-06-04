rm(list=ls())

#-----------------------------------------------------------------
#
# Read data
#
#-----------------------------------------------------------------
URL <- 'https://raw.githubusercontent.com/jordicortes40/constant-effect-RCT/master/'
source(paste0(URL,'code/MA_main_analysis_rma.R'))
source(paste0(URL,'code/SA_I_heuristic.R')) 
source(paste0(URL,'code/SA_II_simulation.R')) 
source(paste0(URL,'code/SA_III_usual_tests.R')) 
source(paste0(URL,'code/SA_IV_distribution_mixture.R')) 

##########################################################################
##-- Table 1
##########################################################################
table1 <- matrix(c(MA_greater_BA,MA_lower_BA,MA_equal_BA,
                   SAI_greater_BA,SAI_lower_BA,SAI_equal_BA,
                   SAII_greater_BA,SAII_lower_BA,SAII_equal_BA,
                   SAIII_greater_BA,SAIII_lower_BA,SAIII_equal_BA,
                   SAIV_greater_BA,SAIV_lower_BA,SAIV_equal_BA, # Lower and greater are heuristic in SA_IV
                   MA_greater_OT,MA_lower_OT,MA_equal_OT,
                   SAI_greater_OT,SAI_lower_OT,SAI_equal_OT,
                   SAIII_greater_OT,SAIII_lower_OT,SAIII_equal_OT,
                   SAIV_greater_OT,SAIV_lower_OT,SAIV_equal_OT), # Lower and greater are heuristic in SA_IV,
                 nrow=9,byrow=TRUE)
table1
prop.table(table1,1)