options(scipen=99999)

xfun::pkg_attach2(c("readxl", "dplyr","haven","stargazer","writexl","lmridge"))
library(mellonMisc)
load("cf_replication_files/Table01.RData")
source("scripts/functions.R")

reps <- 1000
set.seed(456123)
econ.all <- calcEverything(data.r = mydata.econ.r, 
                           data.b = mydata.econ.b, reps = reps)
set.seed(45613)
civ.all <- calcEverything(data.r = mydata.civ.r, 
                          data.b = mydata.civ.b, reps = reps)
set.seed(56653)
energy.all <- calcEverything(data.r = mydata.energy.r, 
                             data.b = mydata.energy.b, reps = reps)
set.seed(61213)
immig.all <- calcEverything(data.r = mydata.imm.r, 
                            data.b = mydata.imm.b, reps = reps)
set.seed(42613)
welfare.all <- calcEverything(data.r = mydata.wel.r, 
                              data.b = mydata.wel.b, reps = reps)
set.seed(15613)
foreign.all <- calcEverything(data.r = mydata.foreign.r,
                              data.b = mydata.foreign.b, reps = reps)

comb.all <- list(Economy = econ.all, "Civil rights"= civ.all, "Energy" = energy.all, 
                 Immigration = immig.all, Welfare = welfare.all, 
                 Foreign = foreign.all)



ols.comparison.resp <- lapply(comb.all, createOLSComparisonRowResp)
ols.comparison.resp <- do.call(rbind, ols.comparison.resp)
ols.comparison.both <- lapply(comb.all, createOLSComparisonRowResp, type = "both", 
                              var = "cooperation")
ols.comparison.both <- do.call(rbind, ols.comparison.both)

falsepos.comb <- dtf(Response = sapply(comb.all, 
                                       function(x) prop.table(table(x$resp.p.null<0.05))["TRUE"]), 
                     Cooperation = sapply(comb.all, 
                                          function(x) prop.table(table(x$coop.p.null<0.05))["TRUE"]))

load(file = "means.by.year.rda")
means.by.year <- means.by.year[sapply(means.by.year, function(x) all(c(min(mydata.econ.r$year), max(mydata.econ.r$year)) %in% x$year))]

cors.vars <- sapply(means.by.year, 
                    function(x) cor(as.matrix(x[, 1:2]), use = "pairwise.complete.obs")[2])
cors <- data.frame(Var = sapply(means.by.year, function(x) colnames(x)[2]), 
                   Correlation = cors.vars)
lm.vars <- sapply(means.by.year, 
                  function(x) summary(lm(x[,2 ]~ x[, 1]))$coefficients[2, "Pr(>|t|)"])

library(lmridge)

set.seed(3485)

issues.inputs.r <- list(Economy = list(model = mresponse.econFINAL, 
                                       data = mydata.econ.r), 
                        `Civil rights` = list(model = mresponse.civFINAL, 
                                              data = mydata.civ.r), 
                        Energy = list(model = mresponse.energyFINAL, 
                                      data = mydata.energy.r), 
                        `Foreign affairs` = list(model = mresponse.foreignFINAL, 
                                                 data = mydata.foreign.r),
                        Immigration = list(model = mresponse.immigrationFINAL, 
                                           data = mydata.imm.r),
                        welfare = list(model = mresponse.welfareFINAL, 
                                       data = mydata.wel.r))

issues.inputs.b <- list(Economy = list(model = mboth.econFINAL, 
                                       data = mydata.econ.b), 
                        `Civil rights` = list(model = mboth.civFINAL, 
                                              data = mydata.civ.b), 
                        Energy = list(model = mboth.energyFINAL, 
                                      data = mydata.energy.b), 
                        `Foreign affairs` = list(model = mboth.foreignFINAL, 
                                                 data = mydata.foreign.b),
                        Immigration = list(model = mboth.immigrationFINAL, 
                                           data = mydata.imm.b),
                        welfare = list(model = mboth.welfareFINAL, 
                                       data = mydata.wel.b))

set.seed(238)
multiple.choices <- c(1,10,100, 200)
power.r <- list()
power.b <- list()
for(ii in names(issues.inputs.r)) {
  print(ii)
  power.r[[ii]] <- replicate(1000, 
                             sapply(multiple.choices, 
                                    simFromResults, 
                                    model = issues.inputs.r[[ii]]$model, 
                                    data = issues.inputs.r[[ii]]$data,
                                    iv= "response"))
  power.b[[ii]] <- replicate(1000, 
                             sapply(multiple.choices, 
                                    simFromResults, 
                                    model = issues.inputs.b[[ii]]$model, 
                                    data = issues.inputs.b[[ii]]$data,
                                    iv= "cooperation"))
}

rr.power <- data.frame("Sample size multiple"=  multiple.choices, 
                       Variable = "Response",
                       sapply(power.r, function(x ) colMeans(t(x)<0.05)) * 100, 
                       check.names = F)

coop.power <- data.frame("Sample size multiple"=  multiple.choices, 
                         Variable = "Cooperation",
                         sapply(power.b, function(x ) colMeans(t(x)<0.05)) * 100, 
                         check.names = F)

comb.power <- rbind(rr.power, coop.power)


effect.size.units.r <- t(sapply(issues.inputs.r,
                                convertCoefToNaturalUnit, var = "response"))
effect.size.units.b <- t(sapply(issues.inputs.b,
                                convertCoefToNaturalUnit, var = "cooperation"))

relevant.ests <- c(effect.size.units.r[c("Immigration", "Economy", "Immigration"), "XYsd"], 
                   effect.size.units.b[c("Immigration", "Economy", "Immigration"), "XYsd"])


vif.all <- rbind(dtf(Issue = names(issues.inputs.r), Variable = "Response", 
                     vif = sapply(issues.inputs.r, vifProblems)),
                 dtf(Issue = names(issues.inputs.r), Variable = "Cooperation", 
                     vif = sapply(issues.inputs.b, vifProblems)))

relevant.vifs <- vif.all[vif.all$Issue %in% c("Economy", "Energy", "Immigration"), "vif"]

library(gridExtra)
library(mellonMisc)
library(lmridge)
library(reshape2)
library(ggplot2)
k.plot.b <- k.plot.r <- list()

for(ii  in names(issues.inputs.r)) {
  k.plot.r[[ii]] <- kPlot(x = issues.inputs.r[[ii]],
                          issue = ii)
  k.plot.b[[ii]] <- kPlot(x = issues.inputs.b[[ii]],
                          issue = ii)
}



save.image("intermediate data/recalculateddata.rdata")
# load("intermediate data/recalculateddata.rdata")




# KM2 
# summary(lmridge(abd_d ~ response + d_congress + year, data = mydata.econ.r, 
#                 scaling = "sc", K = 5.87355))
# KM4
# summary(lmridge(abd_d ~ response + d_congress + year, data = mydata.econ.r, 
#                 scaling = "sc", K = 3.31079))
# KM5
# summary(lmridge(abd_d ~ response + d_congress + year, data = mydata.econ.r, 
#                 scaling = "sc", K = 0.30204))
# KM6
# summary(lmridge(abd_d ~ response + d_congress + year, data = mydata.econ.r, 
#                 scaling = "sc", K = 2.56464))


# different scaling
# summary(lmridge(abd_d ~ response + d_congress + year, data = mydata.econ.r, 
#                 scaling = "sc", K = 3.31079))
# summary(lmridge(abd_d ~ response + d_congress + year, data = mydata.econ.r, 
# scaling = "scaled", K = 3.31079))
# summary(lmridge(abd_d ~ response + d_congress + year, data = mydata.econ.r, 
# scaling = "centered", K = 3.31079))


# Note on the models relating to foriegn affairs ----
# The survey uspew2004-10midpol, which has question relating to foreign affairs only,
# was excluded from the models presented in the paper because of a mistake in the 
# response rates provided by the pollsters in the original data.
# objects relating to foreign affairs with an addendum .M indicates that this survey has been excluded from the data.
# for robustness, we estimate the foreign affairs models a second time, including the survey with inaccurate response rates
# these objects have the same name, without the addendum .M

# Note on dataset objects ----
# objects starting with mydata are datasets used as input in the models.
# the addendum .r indicates the data include a single unit response rate predictor.
# the addendum .b indicates the data includes both contact and cooperation rates.
# the following abbreviations indicate the topic that the data relate to:
# econ = economy
# civ = civil rights
# energy = energy
# imm = immigration
# wel = welfare
# foreign = foreign affairs

# Estimating the correct value of K ----

# the following code creates txt files with various calculations of K to be used in the ridge regressions
# The most appropriate calculation for our data is Muniz et al. 2009 (KM4)
# These estimations are very slowly and memory-demanding
# we therefore have provided the txt files to allow users to view results without estimating 
# users who wish to estimate the results again should uncomment the following lines (mark all, shift+command+c)
# we recommend estimating on Macbook Pro with 32gb RAM. Computer with 16gb are likely to fail, especially the model on foreign data.

# mresponse refers to models with a single unit response rate predictor
# mboth refers to models where response rate is broken into two predictors: cooperation and contact.

# sink("keconomy20200806.txt")
# print("====================================================================================")
# mresponse.econK <- lmridge(abd_d ~ response + d_congress + year, data = mydata.econ.r, scaling = "sc", K = seq(0,10,0.001))
# print("abd_d ~ response + d_congress + year, data = mydata.econ.r")
# kest(mresponse.econK)
# print("====================================================================================")
# mboth.econK <- lmridge(abd_d ~ contact + cooperation + d_congress + year, data = mydata.econ.b, scaling = "sc", K = seq(0,10,0.001))
# print("abd_d ~ contact + cooperation + d_congress + year, data = mydata.econ.b")
# kest(mboth.econK)
# sink()
# 
# sink("kciv20200806.txt")
# print("====================================================================================")
# mresponse.civK <- lmridge(abd_d ~ response + d_congress + year, data = mydata.civ.r, scaling = "sc", K = seq(0,10,0.001))
# print("abd_d ~ response + d_congress + year, data = mydata.civ.r")
# kest(mresponse.civK)
# print("====================================================================================")
# mboth.civK <- lmridge(abd_d ~ contact + cooperation + d_congress + year, data = mydata.civ.b, scaling = "sc", K = seq(0,10,0.001))
# print("abd_d ~ contact + cooperation + d_congress + year, data = mydata.civ.b")
# kest(mboth.civK)
# sink()
# 
# sink("kenergy20200806.txt")
# print("====================================================================================")
# mresponse.energyK <- lmridge(abd_d ~ response + d_congress + year, data = mydata.energy.r, scaling = "sc", K = seq(0,10,0.001))
# print("abd_d ~ response + d_congress + year, data = mydata.energy.r")
# kest(mresponse.energyK)
# print("====================================================================================")
# mboth.energyK <- lmridge(abd_d ~ contact + cooperation + d_congress + year, data = mydata.energy.b, scaling = "sc", K = seq(0,10,0.001))
# print("abd_d ~ contact + cooperation + d_congress + year, data = mydata.energy.b")
# kest(mboth.energyK)
# sink()
# 
# sink("kimm20200806.txt")
# print("====================================================================================")
# mresponse.immK <- lmridge(abd_d ~ response + d_congress + year, data = mydata.imm.r, scaling = "sc", K = seq(0,10,0.001))
# print("abd_d ~ response + d_congress + year, data = mydata.imm.r")
# kest(mresponse.immK)
# print("====================================================================================")
# mboth.immK <- lmridge(abd_d ~ contact + cooperation + d_congress + year, data = mydata.imm.b, scaling = "sc", K = seq(0,10,0.001))
# print("abd_d ~ contact + cooperation + d_congress + year, data = mydata.imm.b")
# kest(mboth.immK)
# sink()
# 
# sink("kwel20200806.txt")
# print("====================================================================================")
# mresponse.welK <- lmridge(abd_d ~ response + d_congress + year, data = mydata.wel.r, scaling = "sc", K = seq(0,10,0.001))
# print("abd_d ~ response + d_congress + year, data = mydata.wel.r")
# kest(mresponse.welK)
# print("====================================================================================")
# mboth.welK <- lmridge(abd_d ~ contact + cooperation + d_congress + year, data = mydata.wel.b, scaling = "sc", K = seq(0,10,0.001))
# print("abd_d ~ contact + cooperation + d_congress + year, data = mydata.wel.b")
# kest(mboth.welK)
# sink()
# 
# sink("kforeign20200806.txt")
# print("====================================================================================")
# mresponse.foreignK <- lmridge(abd_d ~ response + d_congress + year, data = mydata.foreign.r, scaling = "sc", K = seq(0,10,0.001))
# print("abd_d ~ response + d_congress + year, data = mydata.foreign.r")
# kest(mresponse.foreignK)
# print("====================================================================================")
# mboth.foreignK <- lmridge(abd_d ~ contact + cooperation + d_congress + year, data = mydata.foreign.b, scaling = "sc", K = seq(0,10,0.001))
# print("abd_d ~ contact + cooperation + d_congress + year, data = mydata.foreign.b")
# kest(mboth.foreignK)
# sink()
# 
# sink("kforeignM20200806_response.txt")
# print("====================================================================================")
# mresponse.foreignK.m <- lmridge(abd_d ~ response + d_congress + year, data = mydata.foreign.r.m, scaling = "sc", K = seq(0,10,0.001))
# print("abd_d ~ response + d_congress + year, data = mydata.foreign.r.m")
# kest(mresponse.foreignK.m)
# sink()
# 
# sink("kforeignM20200806_both.txt")
# print("====================================================================================")
# mboth.foreignK.m <- lmridge(abd_d ~ contact + cooperation + d_congress + year, data = mydata.foreign.b.m, scaling = "sc", K = seq(0,10,0.001))
# print("abd_d ~ contact + cooperation + d_congress + year, data = mydata.foreign.b.m")
# kest(mboth.foreignK.m)
# sink()

# Estimating Ridge Regression ----

# the following code estimates the final models for each issue area using the correct value of K
# the models have already been estimated and appear as objects in the global environment
# users who wish to re-estiamte the models should uncomment the following lines

# unit response rate as single predictor:
# mresponse.econFINAL <- lmridge(abd_d ~ response + d_congress + year, data = mydata.econ.r, scaling = "sc", K = 3.31079)
# mresponse.civFINAL <- lmridge(abd_d ~ response + d_congress + year, data = mydata.civ.r, scaling = "sc", K = 2.32585)
# mresponse.energyFINAL <- lmridge(abd_d ~ response + d_congress + year, data = mydata.energy.r, scaling = "sc", K = 5.07074)
# mresponse.immigrationFINAL <- lmridge(abd_d ~ response + d_congress + year, data = mydata.imm.r, scaling = "sc", K = 4.71462)
# mresponse.welfareFINAL <- lmridge(abd_d ~ response + d_congress + year, data = mydata.wel.r, scaling = "sc", K = 3.04612)
# mresponse.foreignFINAL <- lmridge(abd_d ~ response + d_congress + year, data = mydata.foreign.r, scaling = "sc", K = 7.74354)
# mresponse.foreignFINAL.M <- lmridge(abd_d ~ response + d_congress + year, data = mydata.foreign.r.m, scaling = "sc", K = 8.00626)

# unit response rate broken into contact and cooperation:
# mboth.econFINAL <- lmridge(abd_d ~ contact + cooperation + d_congress + year, data = mydata.econ.b, scaling = "sc", K = 2.27576)
# mboth.civFINAL <- lmridge(abd_d ~ contact + cooperation + d_congress + year, data = mydata.civ.b, scaling = "sc", K = 2.82527)
# mboth.energyFINAL <- lmridge(abd_d ~ contact + cooperation + d_congress + year, data = mydata.energy.b, scaling = "sc", K = 3.92577)
# mboth.immigrationFINAL <- lmridge(abd_d ~ contact + cooperation + d_congress + year, data = mydata.imm.b, scaling = "sc", K = 4.41655)
# mboth.welfareFINAL <- lmridge(abd_d ~ contact + cooperation + d_congress + year, data = mydata.wel.b, scaling = "sc", K = 1.78701)
# mboth.foreignFINAL <- lmridge(abd_d ~ contact + cooperation + d_congress + year, data = mydata.foreign.b, scaling = "sc", K = 5.87731)
# mboth.foreignFINAL.M <- lmridge(abd_d ~ contact + cooperation + d_congress + year, data = mydata.foreign.b.m, scaling = "sc", K = 5.98371)

# Reviewing Results
summary(mresponse.econFINAL)
summary(mboth.econFINAL)

summary(mresponse.civFINAL)
summary(mboth.civFINAL)

summary(mresponse.energyFINAL)
summary(mboth.energyFINAL)

summary(mresponse.immigrationFINAL)
summary(mboth.immigrationFINAL)

summary(mresponse.welfareFINAL)
summary(mboth.welfareFINAL)

summary(mresponse.foreignFINAL.M)
summary(mboth.foreignFINAL.M)

# Robustness:
summary(mresponse.foreignFINAL)
summary(mboth.foreignFINAL)

