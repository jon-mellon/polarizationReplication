library(quantmod)
library(haven)
gss <-  read_stata("C:/Dropbox/cohort replacement/data/GSS_stata/GSS7218_R1.dta")

miss.vars <- aggregate(!is.na(gss), list(year = gss$year), sum)
gss.15 <- names(which(colSums(miss.vars>0) > 15))
gss.labs <- sapply(gss[,gss.15], 
                   function(x) all(na.omit(attributes(x)$labels) %in% 1:2 )& 
                     length(na.omit(attributes(x)$labels))>0 )
binary.vars.gss <- names(which(gss.labs))

library(readxl)
library(mellonMisc)
gss.cats <- read_excel("additional data/timeseries/gss_var_categorize.xlsx")
gss.cats <- dtf(gss.cats)
gss.cats <- gss.cats[!gss.cats$type=="exclude", ]

binary.vars.gss <- unique(c(gss.cats$var[gss.cats$type=="binary"], binary.vars.gss))
con.vars <- gss.cats[gss.cats$type %in% c("ord", "con"), "var"]

gss.binary <- gss[, binary.vars.gss]==1
gss.binary <- data.frame(gss.binary)
gss.binary <- sapply(gss.binary, as.numeric)
gss.binary <- dtf(gss.binary)
gss.binary$year <- gss$year
gss.binary$id <-1:nrow(gss.binary)
gss.binary$wtssall <- gss$wtssall
gss.binary <- dtf(gss.binary, gss[, con.vars])
gss.b.d <- makeDesign(data = gss.binary, weight.var = "wtssall", id.var = "id")

library(survey)

means.by.year <- lapply(c(con.vars, binary.vars.gss), summarizeVar)

means.by.year <- means.by.year[!sapply(means.by.year, is.null)] 

means.by.year <- means.by.year[sapply(means.by.year, nrow)>0]
means.by.year <- means.by.year[sapply(means.by.year, function(x) sd(x[, 2]))>0]

save(means.by.year, file = "intermediate data/means.by.year.rda")
