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
# source(paste0(URL,'code/SA_I_heuristic.R')) 
# source(paste0(URL,'code/SA_II_simulation.R')) 
source(paste0(URL,'code/SA_III_usual_tests.R')) 
# source(paste0(URL,'code/SA_IV_mixture_distribution.R')) 

##########################################################################
##-- Table 1
##########################################################################
table1 <- matrix(c(lower.var.expF,greater.var.expF,nrow(datos1)-lower.var.expF-greater.var.expF,
                   lower.var.outQ,greater.var.outQ,nrow(datos.paired)-lower.var.outQ-greater.var.outQ),
                 nrow=2,byrow=TRUE)
table1
prop.table(table1,1)