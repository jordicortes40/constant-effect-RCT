rm(list=ls())

quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 
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
# source(paste0(URL,'code/SA_IV_mixture_distribution.R')) 

##########################################################################
##-- Table 1
##########################################################################
table1 <- matrix(c(SAI_greater_BA,SAI_lower_BA,SAI_equal_BA,
                   SAII_greater_BA,SAII_lower_BA,SAII_equal_BA,
                   SAIII_greater_BA,SAIII_lower_BA,SAIII_equal_BA,
                   SAI_greater_OT,SAI_lower_OT,SAI_equal_OT,
                   SAIII_greater_OT,SAIII_lower_OT,SAIII_equal_OT),
                 nrow=5,byrow=TRUE)
table1
prop.table(table1,1)