#-----------------------------------------------------------------
#
# All scripts are in 'http://www-eio.upc.es/teaching/best/variability_data/'
# You can see them in:
# http://www-eio.upc.es/teaching/best/variability_data/Main.R
# http://www-eio.upc.es/teaching/best/variability_data/functions.R
# http://www-eio.upc.es/teaching/best/variability_data/rma_models.R
# http://www-eio.upc.es/teaching/best/variability_data/rma_models_reduced_data.R
# http://www-eio.upc.es/teaching/best/variability_data/subgroups.R
#
#-----------------------------------------------------------------

##-- Remove objects in memory
rm(list=ls())

#-----------------------------------------------------------------
#
# Install and load packages
#
#-----------------------------------------------------------------
##-- Install packages and load packages
list.of.packages <- c('data.table','weights','catspec',
                      'alabama','metafor','epitools',
                      'ggplot2','ggpubr','gridExtra')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
for (pkg in list.of.packages) library(pkg,character.only=TRUE)

##-- Penalize scientific notation
options(scipen=3)

#-----------------------------------------------------------------
#
# Load specific functions
#
#-----------------------------------------------------------------
URL <- 'https://raw.githubusercontent.com/jordicortes40/constant-effect-RCT/master/'
source(paste0(URL,'code/functions.R'))

#-----------------------------------------------------------------
#
# Read data
#
#-----------------------------------------------------------------
datos1 <- read.table(url(paste0(URL,'data/data.csv')),header=TRUE,sep=",",stringsAsFactors = TRUE)
closeAllConnections()

#-----------------------------------------------------------------
#
# Create new variables: outcomes and their standard errors
#
#-----------------------------------------------------------------
##-- Between arms --> Baseline log variance ratio
datos1$yBaselineRatio <- with(datos1,2*log(base_sd_T1/base_sd_T2))                               # Outcome
datos1$seBaselineRatio <- with(datos1,sqrt(2*(1/(base_cases_T1-2)+1/(base_cases_T2-2))))         # Standard error (n-2 = df-1 --> best approximation)


##-- Between arms  --> Outcome log variance ratio
datos1$yBetweenArmsRatio <- with(datos1,2*log(final_sd_T1/final_sd_T2))                          # Outcome
datos1$seBetweenArmsRatio <- with(datos1,sqrt(2*(1/(final_cases_T1-2)+1/(final_cases_T2-2))))    # Standard error (n-2 = df-1 --> best approximation)


##-- Over-time --> Experimental log variance ratio
datos1$yOverTimeRatioT <- with(datos1,2*log(final_sd_T1/base_sd_T1))                             # Outcome Treateds

cor.var <- datos1$rho^2                                                                          # Correlation between variances are aprox. the squared correlation between measures (97 available)
cov.log <- with(datos1,log(1 + 2*cor.var/(final_cases_T1-2)))                                    # Covariance between log variances
datos1$seOverTimeRatioT <- with(datos1,sqrt(2*(1/(final_cases_T1-2) +                            # Standard error (adding covariance) for over-time ratio
                                                 1/(final_cases_T1-2) -
                                                 cov.log)))       

##-- Over-time --> Reference log variance ratio
datos1$yOverTimeRatioC <- with(datos1,2*log(final_sd_T2/base_sd_T2))                             # Outcome Controls
datos1$seOverTimeRatioC <- with(datos1,sqrt(2*(1/(final_cases_T2-2) +                            # Standard error in control group
                                                 1/(final_cases_T2-2) -
                                                 cov.log)))