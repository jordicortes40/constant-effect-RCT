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
# Descriptive
#
#-----------------------------------------------------------------

##########################################################################
##-- Subgroups
##########################################################################
cat('Tables for subgroups----------------------------------------\n')
cat('\nNumber of medical areas-----\n');ctab(factor(datos1$N_Areas_WOS),dec.places=1)
cat('\nIntervention type-----------\n');ctab(datos1$Intervention_type,dec.places=1)
cat('\nOutcome type----------------\n');ctab(datos1$Outcome_type,dec.places=1)
cat('\nCondition type--------------\n');ctab(datos1$Condition_type,dec.places=1)
cat('\nMeasurement type------------\n');ctab(datos1$Measurement_type,dec.places=1)
cat('\nSignificant-----------------\n');ctab(datos1$significant,dec.places = 1)


##########################################################################
##-- Deviations - Table S2
##########################################################################
cat('\nDescriptive Variances with logs----------------------------------------\n')
TableS2 <- with(datos1,matrix(c(fd(2*log(base_sd_T2)),
                                fd(2*log(final_sd_T2)),
                                fd(2*log(base_sd_T1)),
                                fd(2*log(final_sd_T1)),
                                fd(2*log(final_sd_T1)-2*log(final_sd_T2)),
                                fd(2*log(final_sd_T1)-2*log(base_sd_T1)),
                                fd(2*log(final_sd_T1)-2*log(base_sd_T1)-
                                     (2*log(final_sd_T2)-2*log(base_sd_T2)))),nrow=7,byrow=TRUE))
rownames(TableS2) <- paste0('log',paste0('_Var_',c('BasalC','FinalC','BasalT','FinalT','DifT','DifO','DifDif')))
colnames(TableS2) <- c('n','mean','sd','P0%','P25%','P50%','P75%','P100%')
TableS2

##########################################################################
##-- Variances - Descriptive figure between arms
##########################################################################
##-- Previous calculations ------------------------------------------------
ord <- order(datos1$final_sd_T2)
ctrl <- datos1$final_sd_T2[ord]
trt <- datos1$final_sd_T1[ord]
d_figure1 <- data.frame(Treated=trt^2,Control=ctrl^2,Order=order(ctrl),ratio=(trt/ctrl)^2,
                        Group=ifelse(trt>ctrl,'Greater variance in Treated','Greater variance in Control'))


##-- Arrows Graphic ------------------------------------------------
ggplot(d_figure1,aes(x=Order,y=Control,colour=Group)) + geom_point(size=0) + 
  scale_y_log10() +
  geom_segment(aes(x = Order, y = Control, xend = Order, yend = Treated, colour=Group),
               arrow=arrow(length = unit(0.18,"cm")),size=1) +
  theme(legend.position="bottom",
        axis.title = element_text(face='bold'),
        legend.title = element_blank()) +
  xlab('Rank according to Control Outcome variability') + ylab('Variance')

##-- Histogram ------------------------------------------------
ggplot(d_figure1,aes(x=ratio)) + geom_histogram(color='white') + scale_x_log10(limits=c(1/25,25)) +
  ylab('n') + xlab(expression(bold(S[OT]^2/S[OC]^2))) +
  geom_vline(xintercept = 1,linetype=2,color='white',size=1.3) +
  annotate("text",x = 0.1, y = 40, label = "paste(bold(\"Higher \"),bold(S[OC]^2))", parse = TRUE, size=8)+
  annotate("text",x = 10, y = 40, label = "paste(bold(\"Higher \"),bold(S[OT]^2))", parse = TRUE, size=8)+
  theme(axis.title = element_text(face='bold',size=13),
        axis.text = element_text(face='bold',size=13)) 

##########################################################################
##-- Variances - Descriptive figures over-time
##########################################################################
##-- Previous calculations ------------------------------------------------
ord <- order(datos1$base_sd_T1)
baseline <- datos1$base_sd_T1[ord]
outcome <- datos1$final_sd_T1[ord]
d_figure2 <- data.frame(outcome=outcome^2,baseline=baseline^2,Order=order(baseline),ratio=(outcome/baseline)^2,
                        Group=ifelse(outcome>baseline,'Greater variance at the end of the study','Greater variance at baseline'))

##-- Arrows Graphic ------------------------------------------------
ggplot(d_figure2,aes(x=Order,y=baseline,colour=Group)) + geom_point(size=0) + 
  scale_y_log10() +
  geom_segment(aes(x = Order, y = baseline, xend = Order, yend = outcome, colour=Group),
               arrow=arrow(length = unit(0.18,"cm")),size=1) +
  theme(legend.position="bottom",
        axis.title = element_text(face='bold'),
        legend.title = element_blank()) +
  xlab('Rank according to Baseline variability') + ylab('Variance')

##-- Histogram ------------------------------------------------
ggplot(d_figure2,aes(x=ratio)) + geom_histogram(color='white') + scale_x_log10(limits=c(1/25,25)) +
  ylab('n') + xlab(expression(bold(S[OT]^2/S[BT]^2))) +
  geom_vline(xintercept = 1,linetype=2,color='white',size=1.3) +
  annotate("text",x = 0.1, y = 35, label = "paste(bold(\"Higher \"),bold(S[BT]^2))", parse = TRUE, size=8)+
  annotate("text",x = 10, y = 35, label = "paste(bold(\"Higher \"),bold(S[OT]^2))", parse = TRUE, size=8)+
  theme(axis.title = element_text(face='bold',size=13),
        axis.text = element_text(face='bold',size=13)) 
