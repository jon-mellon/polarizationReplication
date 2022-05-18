library(lmridge)
set.seed(1233)
x <- rnorm(150)
y <- x + rnorm(150)
z <- -x + rnorm(150)
dat <- data.frame(x,y,z)
ols.simple <- lm(data = dat, y~x+z)

ks <- seq(0, 10, 0.1)
ridge.vals <- lmridge(formula = y~x+z, data= dat, K = ks)

ridge.out <- data.frame(K = ks, t(ridge.vals$coef / ridge.vals$xscale))
ridge.sum <- summary(ridge.vals)
ridge.sum <- lapply(ridge.sum$summaries, 
                    function(x) x$coefficients[c("x", "z"),
                                               c("Estimate (Sc)", "StdErr (Sc)")])
ridge.sum <- do.call(rbind, ridge.sum)
ridge.sum <- data.frame(ridge.sum)
ridge.sum$K <- inverse.rle(list(values  = ks, lengths= rep(2, length(ks))))
ridge.sum$variable <- sapply(strsplit(rownames(ridge.sum), "\\."), function(x) x[1])
# ridge.out <- melt(ridge.out, id.vars = "K")
xl <- expression(lambda)

ridge.sum$est <- ridge.sum$Estimate..Sc. / ridge.vals$xscale
ridge.sum$se <- ridge.sum$StdErr..Sc. / ridge.vals$xscale
ridge.sum$uci <- ridge.sum$est + ridge.sum$se * 1.96
ridge.sum$lci <- ridge.sum$est - ridge.sum$se * 1.96

save.image("intermediate data/ridgesim.rdata")
