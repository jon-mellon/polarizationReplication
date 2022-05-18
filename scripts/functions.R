

kPlot <- function(x, issue) {
  res <- lmridge(formula = x$model$call$formula, K= seq(0, 10, 0.05),
                 data = x$data)
  
  ses <- sqrt(sapply(vcov(res), diag))
  ses <- ses / res$xscale
  ses <- data.frame(iv = rownames(ses), ses)
  ses <- melt(ses)
  colnames(ses)[colnames(ses)=="variable"] <- "K"
  colnames(ses)[colnames(ses)=="iv"] <- "variable"
  colnames(ses)[colnames(ses)=="value"] <- "se"
  ses$K <- as.character(ses$K)
  ses$K <- as.numeric(gsub("K\\.", "", ses$K))
  
  comb.ests <- data.frame(K = colnames(res$coef), est = t(res$coef / res$xscale))
  comb.ests <- melt(comb.ests, id.var= "K")
  levels(comb.ests$variable) <- gsub("est\\.", "", levels(comb.ests$variable))
  comb.ests$variable <- as.character(comb.ests$variable)
  comb.ests$K <- as.numeric(gsub("K=", "", comb.ests$K))
  
  comb.ests <- safemerge(comb.ests, ses, by = c("K", 'variable'), type = "1:1")
  
  comb.ests$uci <- comb.ests$value + comb.ests$se*1.96
  comb.ests$lci <- comb.ests$value - comb.ests$se*1.96
  comb.ests <- comb.ests[!comb.ests$variable %in% c("d_congress", "contact"), ]
  
  plot <- ggplot(comb.ests, 
                 aes(x = K, y = value, ymin =lci, ymax = uci, 
                     linetype = variable, colour = variable, 
                     shape = variable, fill = variable)) + 
    geom_line() + 
    geom_ribbon(alpha = 0.2, colour = NA) + 
    geom_hline(yintercept = 0, linetype = 3) + 
    theme_bes() + 
    geom_vline(xintercept = x$model$K, linetype = 2) + 
    ylab("Ridge regression coefficient") + ggtitle(issue)  +
    scale_color_manual(values = c("#1F85DE", "#00c317")) + 
    scale_fill_manual(values = c("#1F85DE", "#00c317")) + 
    xlab(expression(lambda))

  return(plot)
}



simpleRidgeTable<- function(ols.comp) {
  ridge.stars <- rep("", nrow(ols.comp))
  ridge.stars[ols.comp$ridge_p<0.05] <- "*"
  ridge.stars[ols.comp$ridge_p<0.01] <- "**"
  ridge.stars[ols.comp$ridge_p<0.001] <- "***"
  
  ols.stars <- rep("", nrow(ols.comp))
  ols.stars[ols.comp$ols_p<0.05] <- "*"
  ols.stars[ols.comp$ols_p<0.01] <- "**"
  ols.stars[ols.comp$ols_p<0.001] <- "***"
  
  out <- data.frame(Issue = rownames(ols.comp), 
                    Ridge = paste(round(ols.comp$ridge, 3), ridge.stars, sep = ""),
                    OLS = paste(round(ols.comp$ols, 3), ols.stars, sep = ""))
  return(out)
}

KM4estimate <- function(object, ...) {
  x <- object$xs
  K = object$K
  y <- object$y
  n <- nrow(x)
  p <- ncol(x)
  sigma2 <- sum(lm.fit(x, y)$residuals^2)/(n - p)
  
  P <- eigen(t(x) %*% x)$vectors
  
  EV <- eigen(t(x) %*% x)$values
  xstar <- x %*% P
  
  alphahat <- solve(diag(EV, p)) %*% t(xstar) %*% y
  mj <- sqrt(sigma2/alphahat^2)
  KM4 <- prod(1/mj)^(1/p)
  return(KM4)
}

summarizeVar <- function(var) {
  print(var)
  agg <- svyby(design = gss.b.d, 
               as.formula(paste("~", var)), 
               by = ~year,
               FUN= svymean, na.rm = T, 
               estimate.only	 = T)
  agg <- agg[agg$se!=0, ]
  return(agg)
}


onePowerSim <- function(data, formula, nullmodel, calc.k = TRUE, coef.measure) {
  rmse <- sqrt(mean(nullmodel$residuals^2)) 
  coef <- nullmodel$coefficients
  
  data$`(Intercept)` <- 1
  preds <- t(coef %*% t(as.matrix(data[, names(coef)])) + rnorm(n = nrow(data), mean = 0, sd = rmse))
  data$nullpred <- preds
  
  if(calc.k) {
    null.mresponse <- lmridge(formula, data = data, 
                              scaling = "sc", 
                              K = seq(0,10,0.001))
    # using KM4 in line with authors:
    k <- KM4estimate(null.mresponse)
  }
  
  null.mresponse <- lmridge(formula, data = data, 
                            scaling = "sc", K = k)
  test <- summary(null.mresponse)
  p <- test$summaries$`summary  1`$coefficients[coef.measure, "Pr(>|t|)"]
  return(p)
}

calcPForNull<- function(data.r, data.b, reps = 1000) {
  nullols <- lm(data = data.r, formula = abd_d ~ d_congress + year)
  p.mod.null <- replicate(reps, onePowerSim(data = data.r,
                                            formula =  nullpred ~ response + d_congress + year, 
                                            nullmodel = nullols,
                                            calc.k = TRUE, 
                                            coef.measure = "response"))
  p.mod.null.cc <- replicate(reps, onePowerSim(data = data.b,
                                               formula =  nullpred ~  contact + cooperation + d_congress + year, 
                                               nullmodel = nullols,
                                               calc.k = TRUE, 
                                               coef.measure = "cooperation"))
  output <- list(resp = p.mod.null, 
                 cc = p.mod.null.cc)
  return(output)
}

calcEverything <- function(data.r, data.b, reps = 1000) {
  p.null <- calcPForNull(data.r = data.r, data.b = data.b, reps = reps)
  
  # lm.response <- lm(abd_d ~ response + d_congress + year, data = data.r)
  # lm.both <- lm(abd_d ~ contact + cooperation + d_congress + year, data = data.b)  
  
  lm.response <- lmridge(abd_d ~ response + d_congress + year, data = data.r,
                         K=0)
  lm.both <- lmridge(abd_d ~ contact + cooperation + d_congress + year, data = data.b,
                     K=0)  
  
  lm.fe.response <- lmridge(abd_d ~ response + d_congress + factor(year), data = data.r,
                            K=0)
  lm.fe.both <- lmridge(abd_d ~ contact + cooperation + d_congress + factor(year), data = data.b,
                        K=0)  
  library(mgcv)
  
  gam.response <- gam(abd_d ~ response +
                        d_congress + s(year), data = data.r)
  gam.both <- gam(abd_d ~ contact + cooperation + 
                    d_congress + s(year), data = data.b)
  library(lme4)
  data.r$meanResponse <- tapply(data.r$response, 
                                data.r$year, mean)[as.character(data.r$year)]
  data.b$meanCooperation <- tapply(data.b$cooperation, 
                                   data.b$year, mean)[as.character(data.b$year)]
  data.b$meanContact <- tapply(data.b$contact, 
                               data.b$year, mean)[as.character(data.b$year)]
  lmer.response <- lmer(abd_d ~ response + meanResponse +
                          d_congress + (1|year), data = data.r)
  lmer.both <- lmer(abd_d ~ contact + cooperation + meanCooperation +
                      meanContact + 
                      d_congress + (1|year), data = data.b)
  
  ridge.k.resp <- lmridge(abd_d ~ response + d_congress + year, data = data.r, 
                          scaling = "sc", K = seq(0,10,0.001))
  
  ridge.final.response <- lmridge(abd_d ~ response + d_congress + year, data = data.r, 
                                  scaling = "sc", K = KM4estimate(ridge.k.resp))
  ridge.k.both <- lmridge(abd_d ~ contact + cooperation + d_congress + year, data = data.b, 
                          scaling = "sc", K = seq(0,10,0.001))
  ridge.final.both <- lmridge(abd_d ~ contact + cooperation + d_congress + year, data = data.b, 
                              scaling = "sc", K = KM4estimate(ridge.k.both))
  
  output <- list(lm.both = lm.both, 
                 lm.response = lm.response, 
                 lm.fe.response = lm.fe.response,
                 lm.fe.both = lm.fe.both,
                 lmer.response= lmer.response,
                 lmer.both = lmer.both,
                 gam.response = gam.response,
                 gam.both = gam.both,
                 ridge.final.response = ridge.final.response,
                 ridge.final.both = ridge.final.both,
                 resp.p.null = p.null$resp, 
                 coop.p.null = p.null$cc)
  return(output)
}

# createOLSComparisonRowResp <- function(all) {
#   ridge <- summary(all$ridge.final.resp)$summaries$`summary  1`$coefficient
#   ols <- summary(all$lm.response)$summaries$`summary  1`$coefficients
#   
#   row <- data.frame(ridge = ridge["response", "Estimate (Sc)"],
#                     ridge_se = ridge["response", "StdErr (Sc)"], 
#                     ridge_p = ridge["response", "Pr(>|t|)"], 
#                     ols = ols["response", "Estimate (Sc)"],
#                     ols_se = ols["response", "StdErr (Sc)"], 
#                     ols_p = ols["response", "Pr(>|t|)"])
#   return(row)  
# }
z.to.p <- function(z) {
  p <- pnorm(z, mean = 0, sd= 1, lower.tail = F) * 2  
  return(p)
}

createOLSComparisonRowBoth <- function(all) {
  s <- summary(all$ridge.final.both)$summaries$`summary  1`$coefficients
  ols <- summary(all$lm.both)$summaries$`summary  1`$coefficients
  
  row <- data.frame(ridge = s["cooperation", "Estimate (Sc)"],
                    ridge_se = s["cooperation", "StdErr (Sc)"], 
                    ridge_p = s["cooperation", "Pr(>|t|)"], 
                    ols = ols["cooperation", "Estimate (Sc)"], 
                    ols_se = ols["cooperation", "StdErr (Sc)"], 
                    ols_p = ols["cooperation", "Pr(>|t|)"])
  return(row)  
}


createOLSComparisonRowResp <- function(x, type = "response", var = "response") {
  
  extractLMCoefs <- function(z) {
    out <- dtf(summary(z)$summaries$`summary  1`$coefficients[-1, c("Estimate (Sc)", "StdErr (Sc)")] / z$xscale, 
               p = summary(z)$summaries$`summary  1`$coefficients[-1, "Pr(>|t|)"])
    colnames(out) <- c("est", "se", "p")
    return(out)
  }
  ols <- extractLMCoefs(x[[paste0("lm.", type)]])
  ols$Model <- "OLS (time trend)"
  fe <- extractLMCoefs(x[[paste0("lm.fe.", type)]])
  fe$Model <- "Year FEs"
  
  ridge <- extractLMCoefs(x[[paste0("ridge.final.", type)]])
  ridge$Model <- "Ridge"
  hlm <- dtf(summary(x[[paste0("lmer.", type)]])$coefficients[, c("Estimate", "Std. Error")],
             p = z.to.p(summary(x[[paste0("lmer.", type)]])$coefficients[, "t value"]))
  colnames(hlm) <- c("est", "se", "p")
  hlm$Model <- "Year REs"
  
  gam.temp <- x[[paste0("gam.", type)]]
  gam.temp2 <- summary(gam.temp)
  gam <- dtf(est = gam.temp$coefficients[names(gam.temp2$p.coeff)], 
             se = gam.temp2$se[names(gam.temp2$p.coeff)], 
             p = gam.temp2$p.coeff)
  
  gam$Model <- "GAM time trend"
  results <- rbind(ridge[var, ], 
                   ols[var, ], 
                   fe[var, ], 
                   hlm[var, ], 
                   gam[var, ])
  results$star <- ""
  results$star[results$p <0.05] <- "*"
  results$star[results$p <0.01] <- "**"
  results$star[results$p <0.001] <- "***"
  
  results$est.present <- paste0(round(results$est, 4), results$star)
  return(results)
}

simFromResults <- function(model, data, iv, multiples) {
  sigma <- sd(predict(model) - model$y)
  data <- data[rep(1:nrow(data), multiples), ]
  sim.y <- predict(model) + 
    rnorm(nrow(data),
          sd = sigma)
  data$abd_d <- sim.y
  p.out <- summary(lm(formula = model$call$formula, 
                      data = data))$coefficients[iv, "Pr(>|t|)"]
  return(p.out)
}

convertCoefToNaturalUnit <- function(x, var) {
  natural.units <- (x$model$coef / x$model$xscale)[var, ]
  Ysd.units <- (natural.units /   sd(x$data$abd_d))
  XYsd.units <- Ysd.units * sd(x$data[, var, drop = T])  
  out <- c(natural = natural.units,
           Ysd = Ysd.units, 
           XYsd = XYsd.units)
  return(out)
}


vifProblems <- function(x) {
  r2 <- summary(lm(x$model$xs[, 1] ~ x$model$xs[, -1]))$r.squared
  vif <- 1 / (1-r2)  
  return(vif)
}