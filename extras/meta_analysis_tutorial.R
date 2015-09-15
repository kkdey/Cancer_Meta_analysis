

###################   Meta Analysis Tutorial  #############################


########################   R packages   #################################

library(rmeta)
library(meta)
library(metafor)

#########################  Loading Data ################################

#############  BCG vaccine data

data(dat.bcg)
str(dat.bcg)

#############  Amlodipine data

data(amlodipine) # amlodipine meta-analytic dataset
str(amlodipine) # Describe meta-analysis structure

#####################   Log Odds Ratio  ##############################

Y <- with(dat.bcg, log(tpos * cneg/(tneg * cpos)))
V <- with(dat.bcg, 1/tpos + 1/cneg + 1/tneg + 1/cpos)
cbind(Y, V)

####################  Fixed effects model  ############################

result.or <- rma(yi = Y, vi = V, method = "FE") # Log Odds Ratio
summary(result.or)

###  contributions ####

contributions <- 1/result.or$vi/sum(1/result.or$vi) * 100
cbind(contributions)

par(mar = c(5, 10, 5, 5))
barplot(contributions, names = dat.bcg$trial,
xlim = c(0, 50), las = 2, horiz = T,
col = "royalblue", xlab="Contribution", ylab="Trials")

coef(result.or)

##############   Random effects model ##################################

###  Mantel random effects test

result.DL <- rma(yi = Y, vi = V, method = "DL") # Log Odds Ratio
summary(result.DL)

###  REML (Restricted Max Likelihood Estimation)  

result.REML <- rma(yi = Y, vi = V, method = "REML") # Log Odds Ratio
summary(result.REML)

result.REML$QEp  ########  p-value of the heterogeneity

result.REML$QE   #########   Q test statistic heterogeneity

confint(result.REML) ########  heterogeneity test CI 


#################   Visulaizing Heterogeneity ###########################

#################   The Forest Plot  ###############################


forest(result.REML)
dat.bcg$n <- with(dat.bcg, tpos + tneg + cpos + cneg)
forest(result.or, order = "prec",ilab = dat.bcg[, c("n", "year")],ilab.xpos = exp(result.or$b) - c(4, 2),transf = exp, refline = 1)


#####################  Sensitivity Analysis  ############################

leave1out(result.REML)

#####################  Funnel Plot : Publication Bias ###################3


funnel(result.REML,addtau2=TRUE)
trimfill(result.REML)
funnel(trimfill(result.REML))
value <- fsn(y = result.or$yi, v = result.or$vi)

