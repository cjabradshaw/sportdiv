## Global wildlife diversity in team sport organisations
## Corey Bradshaw & Ugo Arbieu

## libraries
library(lme4)
library(Hmisc)
library(boot)
library(ggpubr)
library(performance)
library(sjPlot)
library(ggplot2)

# source files
source("new_lmer_AIC_tables3.R") 
source("r.squared.R") 

# functions
AICc <- function(...) {
  models <- list(...)
  num.mod <- length(models)
  AICcs <- numeric(num.mod)
  ns <- numeric(num.mod)
  ks <- numeric(num.mod)
  AICc.vec <- rep(0,num.mod)
  for (i in 1:num.mod) {
    if (length(models[[i]]$df.residual) == 0) n <- models[[i]]$dims$N else n <- length(models[[i]]$residuals)
    if (length(models[[i]]$df.residual) == 0) k <- sum(models[[i]]$dims$ncol) else k <- (length(models[[i]]$coeff))+1
    AICcs[i] <- (-2*logLik(models[[i]])) + ((2*k*n)/(n-k-1))
    ns[i] <- n
    ks[i] <- k
    AICc.vec[i] <- AICcs[i]
  }
  return(AICc.vec)
}

## functions
delta.IC <- function(x) x - min(x) ## where x is a vector of an IC
weight.IC <- function(x) (exp(-0.5*x))/sum(exp(-0.5*x)) ## Where x is a vector of dIC
ch.dev <- function(x) ((( as.numeric(x$null.deviance) - as.numeric(x$deviance) )/ as.numeric(x$null.deviance))*100) ## % change in deviance, where x is glm object

linreg.ER <- function(x,y) { # where x and y are vectors of the same length; calls AICc, delta.AIC, weight.AIC functions
  fit.full <- lm(y ~ x); fit.null <- lm(y ~ 1)
  AIC.vec <- c(AICc(fit.full),AICc(fit.null))
  dAIC.vec <- delta.AIC(AIC.vec); wAIC.vec <- weight.AIC(dAIC.vec)
  ER <- wAIC.vec[1]/wAIC.vec[2]
  r.sq.adj <- as.numeric(summary(fit.full)[9])
  return(c(ER,r.sq.adj))
}

## Professional leagues in each country. To compare women’s and men’s differences in 'Percent2'
# TotTeams = total number of teams in the league
# AnimTeams = number of teams with animal (domesticated + wild) emblems in the league
# Percent = TotTeams / AnimTeams
# naf_all = total number of animal organisms found in the league
# naf_wild = total number of wildlife organisms found in the league
# AnimTeamsWild = number of teams with wildlife emblems in the league
# Percent2 = AnimTeamsWild / TotTeams

# load the data
league.dat <- read.csv("leagueDat.csv")
head(league.dat)

# scale & transform
hist(league.dat$Percent2, breaks=20, col="grey", xlab="% of teams with wildlife emblems", main="")
twmwldembllgtcor1 <- ifelse(league.dat$Percent2==0, 0.05, league.dat$Percent2/100)
twmwldembllgtcor <- ifelse(twmwldembllgtcor1==1, 0.95, twmwldembllgtcor1)
league.dat$tmwldembl.sc <- scale(logit(twmwldembllgtcor), scale=T, center=T)
hist(league.dat$tmwldembl.sc, col="grey", xlab="scaled logit prop of teams with wildlife emblems", main="")

table(league.dat$Gender)
table(league.dat$Continent)
xtabs(league.dat$tmwldembl.sc ~ league.dat$Continent) / table(league.dat$Continent)

ggdensity(league.dat$tmwldembl.sc, xlab = "scaled logit prop of teams with wildlife emblems")
ggqqplot(league.dat$tmwldembl.sc)
shapiro.test(league.dat$tmwldembl.sc)


## generalised linear mixed-effects models
# model set
m1 <- "tmwldembl.sc ~ Gender + Sport + Gender*Sport + (1|Continent)"
m2 <- "tmwldembl.sc ~ Gender + Sport + (1|Continent)"
m3 <- "tmwldembl.sc ~ Gender + (1|Continent)"
m4 <- "tmwldembl.sc ~ Sport + (1|Continent)"
m5 <- "tmwldembl.sc ~ 1 + (1|Continent)"

## make model vector
mod.vec <- c(m1,m2,m3,m4,m5)

## define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- lmer(as.formula(mod.vec[i]),data=league.dat, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- as.numeric(logLik(fit))
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- r.squared(fit)$AIC
  BIC.vec[i] <- BIC(fit)
  
  Rm[i] <- 100*r.squared(fit)$Marginal # marginal R-squared
  Rc[i] <- 100*r.squared(fit)$Conditional # conditional R-squared
  
  print(i)
}
dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)
dBIC <- delta.IC(BIC.vec)
wBIC <- weight.IC(dBIC)

sumtable <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,4),round(Rm,4),round(Rc,4))
colnames(sumtable) <- c("model","k","LL","AICc","dAICc","wAICc","Rm","Rc")
row.names(sumtable) <- as.character(mod.vec)
summary.table <- sumtable[order(sumtable[,6],decreasing=T),]
summary.table

## different structure (no reliance on proportions)
hist(league.dat$AnimTeamsWild)
hist(log10(league.dat$AnimTeamsWild+0.5))
hist(league.dat$TotTeams)
hist(log10(league.dat$TotTeams))


## generalised linear mixed-effects models
# model set
m1 <- "AnimTeamsWild ~ TotTeams + Gender + TotTeams*Gender + Sport + Gender*Sport + (1|Continent)"
m2 <- "AnimTeamsWild ~ TotTeams + Gender + Sport + Gender*Sport + (1|Continent)"
m3 <- "AnimTeamsWild ~ TotTeams + Gender + Sport + (1|Continent)"
m4 <- "AnimTeamsWild ~ TotTeams + Gender + TotTeams*Gender + (1|Continent)"
m5 <- "AnimTeamsWild ~ TotTeams + Gender + (1|Continent)"
m6 <- "AnimTeamsWild ~ TotTeams + Sport + (1|Continent)"
m7 <- "AnimTeamsWild ~ TotTeams + (1|Continent)"
m8 <- "AnimTeamsWild ~ 1 + (1|Continent)"

## make model vector
mod.vec <- as.character(mget(paste("m",rep(1:8,1),sep="")))

## define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- lmer(as.formula(mod.vec[i]),data=league.dat, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- as.numeric(logLik(fit))
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- r.squared(fit)$AIC
  BIC.vec[i] <- BIC(fit)
  
  Rm[i] <- 100*r.squared(fit)$Marginal # marginal R-squared
  Rc[i] <- 100*r.squared(fit)$Conditional # conditional R-squared
  
  print(i)
}
dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)
dBIC <- delta.IC(BIC.vec)
wBIC <- weight.IC(dBIC)

sumtable <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,4),round(Rm,4),round(Rc,4))
colnames(sumtable) <- c("model","k","LL","AICc","dAICc","wAICc","Rm","Rc")
row.names(sumtable) <- as.character(mod.vec)
summary.table <- sumtable[order(sumtable[,6],decreasing=T),]
summary.table


## log10
# define models
m1 <- "log10(AnimTeamsWild+0.5) ~ log10(TotTeams) + Gender + log10(TotTeams)*Gender + Sport + Gender*Sport + (1|Continent)"
m2 <- "log10(AnimTeamsWild+0.5) ~ log10(TotTeams) + Gender + Sport + Gender*Sport + (1|Continent)"
m3 <- "log10(AnimTeamsWild+0.5) ~ log10(TotTeams) + Gender + Sport + (1|Continent)"
m4 <- "log10(AnimTeamsWild+0.5) ~ log10(TotTeams) + Gender + log10(TotTeams)*Gender + (1|Continent)"
m5 <- "log10(AnimTeamsWild+0.5) ~ log10(TotTeams) + Gender + (1|Continent)"
m6 <- "log10(AnimTeamsWild+0.5) ~ log10(TotTeams) + Sport + (1|Continent)"
m7 <- "log10(AnimTeamsWild+0.5) ~ log10(TotTeams) + (1|Continent)"
m8 <- "log10(AnimTeamsWild+0.5) ~ 1 + (1|Continent)"

## make model vector
mod.vec <- as.character(mget(paste("m",rep(1:8,1),sep="")))

## define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- lmer(as.formula(mod.vec[i]),data=league.dat, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- as.numeric(logLik(fit))
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- r.squared(fit)$AIC
  BIC.vec[i] <- BIC(fit)
  
  Rm[i] <- 100*r.squared(fit)$Marginal # marginal R-squared
  Rc[i] <- 100*r.squared(fit)$Conditional # conditional R-squared
  
  print(i)
}
dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)
dBIC <- delta.IC(BIC.vec)
wBIC <- weight.IC(dBIC)

sumtable <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,4),round(Rm,4),round(Rc,4))
colnames(sumtable) <- c("model","k","LL","AICc","dAICc","wAICc","Rm","Rc")
row.names(sumtable) <- as.character(mod.vec)
summary.table <- sumtable[order(sumtable[,6],decreasing=T),]
summary.table
check_model(fit)
plot_model(fit, show.values=T, type="re", vline.color = "purple")


# general linear models
# define models
m1 <- "log10(AnimTeamsWild+0.5) ~ log10(TotTeams) + Gender + log10(TotTeams)*Gender + Sport + Gender*Sport"
m2 <- "log10(AnimTeamsWild+0.5) ~ log10(TotTeams) + Gender + Sport + Gender*Sport"
m3 <- "log10(AnimTeamsWild+0.5) ~ log10(TotTeams) + Gender + Sport"
m4 <- "log10(AnimTeamsWild+0.5) ~ log10(TotTeams) + Gender + log10(TotTeams)*Gender"
m5 <- "log10(AnimTeamsWild+0.5) ~ log10(TotTeams) + Gender"
m6 <- "log10(AnimTeamsWild+0.5) ~ log10(TotTeams) + Sport"
m7 <- "log10(AnimTeamsWild+0.5) ~ log10(TotTeams)"
m8 <- "log10(AnimTeamsWild+0.5) ~ 1"

## make model vector
mod.vec <- as.character(mget(paste("m",rep(1:8,1),sep="")))

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
mod.list <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=league.dat, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,5],decreasing=T),]
summary.table

i <- 1
fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=league.dat, na.action=na.omit)
plot(fit)
qqnorm(fit$residuals)
check_model(fit)
plot_model(fit, show.values=T, vline.color = "purple")

## without 'Sport'
# define models
m1 <- "log10(AnimTeamsWild+0.5) ~ log10(TotTeams) + Gender + log10(TotTeams)*Gender"
m2 <- "log10(AnimTeamsWild+0.5) ~ log10(TotTeams) + Gender"
m3 <- "log10(AnimTeamsWild+0.5) ~ log10(TotTeams)"
m4 <- "log10(AnimTeamsWild+0.5) ~ Gender"
m5 <- "log10(AnimTeamsWild+0.5) ~ 1"

## make model vector
mod.vec <- as.character(mget(paste("m",rep(1:5,1),sep="")))

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
mod.list <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=league.dat, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,5],decreasing=T),]
summary.table

i <- 1
fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=league.dat, na.action=na.omit)
plot(fit)
qqnorm(fit$residuals)
check_model(fit)
plot_model(fit, show.values=T, vline.color = "purple")

## no evidence for a gender effect on AnimTeamsWild/TotTeams


########################################
## Database focusing on the taxa that could be identified at the species level
## “Classification” columns = finest possible level of identification
# “N_Species_In_Country” =  number of corresponding species that occur in the respective country,
#   according to the IUCN database
# “EX” … “NE” = IUCN threat levels &  corresponding number of species in that country
## “Decreasing” … “Unknown” = population trends according to the IUCN database

## “N_Species_In_Country”  used for analysis of presence/absence (0=absent, 1=present) in the country

# load data
spp.dat <- read.csv("sppDat.csv")
head(spp.dat)


## differences in class representation
table(spp.dat$Continent)
table(spp.dat$Country)
table(spp.dat$Classification_Class)
table(spp.dat$N_Species_In_Country)
table(spp.dat$N_Species_In_Country, spp.dat$Continent, spp.dat$Classification_Class)
spp.dat$ClassComb <- ifelse(spp.dat$Classification_Class == "Actinopterygii" | 
                              spp.dat$Classification_Class == "Arachnida" |
                              spp.dat$Classification_Class == "Chondrichthyes" |
                              spp.dat$Classification_Class == "Reptilia" |
                              spp.dat$Classification_Class == "Insecta", "other", 
                            spp.dat$Classification_Class)
table(spp.dat$ClassComb)
table(spp.dat$Continent, spp.dat$ClassComb)

## simulation test of independence
class.tab <- table(spp.dat$Continent, spp.dat$ClassComb)
class.tab
contSums <- rowSums(class.tab)
classSums <- colSums(class.tab)

# exact binomial test
binom.test(classSums[2], sum(classSums), 1/3, alternative="greater")


exp.tab <- class.tab
for (i in 1:length(contSums)) {
  exp.tab[i,1] <- classSums[1] * (contSums[i]/sum(class.tab))
  exp.tab[i,2] <- classSums[2] * (contSums[i]/sum(class.tab))
  exp.tab[i,3] <- classSums[3] * (contSums[i]/sum(class.tab))
}
exp.tab
chi.tab <- (class.tab - exp.tab)^2/exp.tab
chi.obs <- sum(chi.tab)
chisq.test(class.tab, simulate.p.value=T, B=10000)


iter <- 100000
itdiv <- iter/10
chi.rnd <- rep(0,iter)

for (c in 1:iter) {
  rnd.tab <- class.tab
  chi.tab.adj <- chi.tab
  
  # AFR
  rnd.tab[1,1] <- rbinom(1,contSums[1],classSums[1]/sum(classSums))
  rnd.tab[1,2] <- rbinom(1,contSums[1],classSums[2]/sum(classSums))
  rnd.tab[1,3] <- contSums[1] - sum(rnd.tab[1,1], rnd.tab[1,2])
  mam.less <- ifelse(rnd.tab[1,3] < 0, rbinom(1, abs(rnd.tab[1,3]), classSums[2]/sum(classSums[1:2])), 0)
  ave.less <- ifelse(mam.less < abs(rnd.tab[1,3]), abs(rnd.tab[1,3]) - mam.less, 0)
  rnd.tab[1,2] <- rnd.tab[1,2] - mam.less
  rnd.tab[1,1] <- ifelse(rnd.tab[1,3] > 0, rnd.tab[1,1], rnd.tab[1,1] - ave.less)
  rnd.tab[1,3] <- ifelse(rnd.tab[1,3] < 0, 0, rnd.tab[1,3])
  
  # AMR
  rnd.tab[2,1] <- rbinom(1,contSums[2],classSums[1]/sum(classSums))
  rnd.tab[2,2] <- rbinom(1,contSums[2],classSums[2]/sum(classSums))
  rnd.tab[2,3] <- contSums[2] - sum(rnd.tab[2,1], rnd.tab[2,2])
  mam.less <- ifelse(rnd.tab[2,3] < 0, rbinom(1, abs(rnd.tab[2,3]), classSums[2]/sum(classSums[1:2])), 0)
  ave.less <- ifelse(mam.less < abs(rnd.tab[2,3]), abs(rnd.tab[2,3]) - mam.less, 0)
  rnd.tab[2,2] <- rnd.tab[2,2] - mam.less
  rnd.tab[2,1] <- ifelse(rnd.tab[2,3] > 0, rnd.tab[2,1], rnd.tab[2,1] - ave.less)
  rnd.tab[2,3] <- ifelse(rnd.tab[2,3] < 0, 0, rnd.tab[2,3])
  
  # ASI
  rnd.tab[3,1] <- rbinom(1,contSums[3],classSums[1]/sum(classSums))
  rnd.tab[3,2] <- rbinom(1,contSums[3],classSums[2]/sum(classSums))
  rnd.tab[3,3] <- contSums[3] - sum(rnd.tab[3,1], rnd.tab[3,2])
  mam.less <- ifelse(rnd.tab[3,3] < 0, rbinom(1, abs(rnd.tab[3,3]), classSums[2]/sum(classSums[1:2])), 0)
  ave.less <- ifelse(mam.less < abs(rnd.tab[3,3]), abs(rnd.tab[3,3]) - mam.less, 0)
  rnd.tab[3,2] <- rnd.tab[3,2] - mam.less
  rnd.tab[3,1] <- ifelse(rnd.tab[3,3] > 0, rnd.tab[3,1], rnd.tab[3,1] - ave.less)
  rnd.tab[3,3] <- ifelse(rnd.tab[3,3] < 0, 0, rnd.tab[3,3])
  
  # EUR
  rnd.tab[4,1] <- rbinom(1,contSums[4],classSums[1]/sum(classSums))
  rnd.tab[4,2] <- rbinom(1,contSums[4],classSums[2]/sum(classSums))
  rnd.tab[4,3] <- contSums[4] - sum(rnd.tab[4,1], rnd.tab[4,2])
  mam.less <- ifelse(rnd.tab[4,3] < 0, rbinom(1, abs(rnd.tab[4,3]), classSums[2]/sum(classSums[1:2])), 0)
  ave.less <- ifelse(mam.less < abs(rnd.tab[5,3]), abs(rnd.tab[4,3]) - mam.less, 0)
  rnd.tab[4,2] <- rnd.tab[4,2] - mam.less
  rnd.tab[4,1] <- ifelse(rnd.tab[4,3] > 0, rnd.tab[4,1], rnd.tab[4,1] - ave.less)
  rnd.tab[4,3] <- ifelse(rnd.tab[4,3] < 0, 0, rnd.tab[4,3])
  
  # OCE
  rnd.tab[5,1] <- rbinom(1,contSums[5],classSums[1]/sum(classSums))
  rnd.tab[5,2] <- rbinom(1,contSums[5],classSums[2]/sum(classSums))
  rnd.tab[5,3] <- contSums[5] - sum(rnd.tab[5,1], rnd.tab[5,2])
  mam.less <- ifelse(rnd.tab[5,3] < 0, rbinom(1, abs(rnd.tab[5,3]), classSums[2]/sum(classSums[1:2])), 0)
  ave.less <- ifelse(mam.less < abs(rnd.tab[5,3]), abs(rnd.tab[5,3]) - mam.less, 0)
  rnd.tab[5,2] <- rnd.tab[5,2] - mam.less
  rnd.tab[5,1] <- ifelse(rnd.tab[5,3] > 0, rnd.tab[5,1], rnd.tab[5,1] - ave.less)
  rnd.tab[5,3] <- ifelse(rnd.tab[5,3] < 0, 0, rnd.tab[5,3])
  
  chi.tab.rnd <- (rnd.tab - exp.tab)^2/exp.tab
  inf.sub <- which(is.infinite(chi.tab.rnd)==T, arr.ind=T)
  if (length(inf.sub) > 0) {
    chi.tab.rnd[inf.sub] <- 0
    chi.tab.adj[inf.sub] <- 0
  } # end if  
    
  chi.rnd[c] <- sum(chi.tab.rnd, na.rm=T)
  comp.vec <- ifelse(chi.rnd[c] >= sum(chi.tab.adj, na.rm=T), 1, 0)
  
  if (c %% itdiv==0) print(c)
  
} # end c
pr <- sum(comp.vec)/iter
pr # low type I error for hypothesis that distribution of species by major class differs by continent


## conservation status
## simulation test of independence
head(spp.dat)
threat.sub <- spp.dat[,which(colnames(spp.dat) %in% c("EX","EW","CR","EN","VU","NT","LC","DD","NE"))]
dim(threat.sub)
threat.stat <- which(threat.sub == 1, arr.ind=T)
nothreat.sub <- which(rowSums(threat.sub, na.rm=T) == 0)
threat.stat.ord <- threat.stat[order(threat.stat[,1],decreasing=F),]
threat <- ifelse(threat.stat.ord[,2] == 1, "EX", 
                          ifelse(threat.stat.ord[,2] == 2, "EW", 
                                 ifelse(threat.stat.ord[,2] == 3, "CR", 
                                        ifelse(threat.stat.ord[,2] == 4, "EN", 
                                               ifelse(threat.stat.ord[,2] == 5, "VU", 
                                                      ifelse(threat.stat.ord[,2] == 6, "NT", 
                                                             ifelse(threat.stat.ord[,2] == 7, "LC", 
                                                                    ifelse(threat.stat.ord[,2] == 8, "DD",
                                                                           ifelse(threat.stat.ord[,2] == 9, "NE", NA)))))))))
dim(spp.dat)
length(threat)
spp.red <- spp.dat[-nothreat.sub,]
spp.red$IUCN <- threat
spp.red$threat <- ifelse(spp.red$IUCN == "CR" | spp.red$IUCN == "EN" | spp.red$IUCN == "VU", 1, 0)
threat.tab <- table(spp.red$Continent, spp.red$threat)
threat.tab
threat.tab[,1]/rowSums(threat.tab)
contSums <- rowSums(threat.tab)
threatSums <- colSums(threat.tab)

# exact binomial test
binom.test(threatSums[2], sum(threatSums), 1/2, alternative="two.sided")

exp.tab <- threat.tab
for (i in 1:length(contSums)) {
  exp.tab[i,1] <- threatSums[1] * (contSums[i]/sum(threat.tab))
  exp.tab[i,2] <- threatSums[2] * (contSums[i]/sum(threat.tab))
}
exp.tab
rowSums(exp.tab)
chi.tab <- (threat.tab - exp.tab)^2/exp.tab
chi.obs <- sum(chi.tab)
chisq.test(threat.tab, simulate.p.value=T, B=10000)

iter <- 100000
itdiv <- iter/10
chi.rnd <- rep(NA,iter)

for (c in 1:iter) {
  rnd.tab <- threat.tab
  
  # AFR
  rnd.tab[1,1] <- rbinom(1,contSums[1],threatSums[1]/sum(threatSums))
  rnd.tab[1,2] <- contSums[1] - rnd.tab[1,1]

  # AMR
  rnd.tab[2,1] <- rbinom(1,contSums[2],threatSums[1]/sum(threatSums))
  rnd.tab[2,2] <- contSums[2] - rnd.tab[2,1]
  
  # ASI
  rnd.tab[3,1] <- rbinom(1,contSums[3],threatSums[1]/sum(threatSums))
  rnd.tab[3,2] <- contSums[3] - rnd.tab[3,1]

  # EUR
  rnd.tab[4,1] <- rbinom(1,contSums[4],threatSums[1]/sum(threatSums))
  rnd.tab[4,2] <- contSums[4] - rnd.tab[4,1]

  # OCE
  rnd.tab[5,1] <- rbinom(1,contSums[5],threatSums[1]/sum(threatSums))
  rnd.tab[5,2] <- contSums[5] - rnd.tab[5,1]

  
  chi.tab.rnd <- (rnd.tab - exp.tab)^2/exp.tab
  chi.rnd[c] <- sum(chi.tab.rnd, na.rm=T)

  if (c %% itdiv==0) print(c)
  
} # end c
pr <- length(which(chi.rnd >= chi.obs))/iter
pr # low type I error for hypothesis that distribution of species by major threat category (threatened or not) differs by continent


## population trend
head(spp.dat)
trend.sub <- spp.dat[,which(colnames(spp.dat) %in% c("Decreasing","Increasing","Stable","Unknown"))]
head(trend.sub)
dim(trend.sub)
notrend.sub <- which(rowSums(trend.sub, na.rm=T) == 0)
trend.stat <- which(trend.sub == 1, arr.ind=T)
trend.stat.ord <- trend.stat[order(trend.stat[,1],decreasing=F),]
trend <- ifelse(trend.stat.ord[,2] == 1, "decr", 
                 ifelse(trend.stat.ord[,2] == 2, "incr", 
                        ifelse(trend.stat.ord[,2] == 3, "stab", 
                               ifelse(trend.stat.ord[,2] == 4, "unkn", NA))))
trend
dim(spp.dat)
length(trend)
spp.red <- spp.dat[-notrend.sub,]
spp.red$trend <- trend
spp.red2 <- spp.red[spp.red$trend != "unkn", ]
dim(spp.red2)
trend.tab <- table(spp.red2$Continent, spp.red2$trend)
trend.tab
contSums <- rowSums(trend.tab)
trendSums <- colSums(trend.tab)

# exact binomial test
binom.test(trendSums[1], sum(threatSums), 1/3, alternative="greater")

exp.tab <- trend.tab
for (i in 1:length(contSums)) {
  exp.tab[i,1] <- trendSums[1] * (contSums[i]/sum(trend.tab))
  exp.tab[i,2] <- trendSums[2] * (contSums[i]/sum(trend.tab))
  exp.tab[i,3] <- trendSums[3] * (contSums[i]/sum(trend.tab))
}
exp.tab
trend.tab
chi.tab <- (trend.tab - exp.tab)^2/exp.tab
chi.obs <- sum(chi.tab)
chi.obs
chisq.test(trend.tab, simulate.p.value=T, B=10000)


iter <- 100000
itdiv <- iter/10
chi.rnd <- rep(0,iter)

for (c in 1:iter) {
  rnd.tab <- trend.tab
  chi.tab.adj <- chi.tab
  
  # AFR
  rnd.tab[1,1] <- rbinom(1,contSums[1],trendSums[1]/sum(trendSums))
  rnd.tab[1,2] <- rbinom(1,contSums[1],trendSums[2]/sum(trendSums))
  rnd.tab[1,3] <- contSums[1] - sum(rnd.tab[1,1], rnd.tab[1,2])
  mam.less <- ifelse(rnd.tab[1,3] < 0, rbinom(1, abs(rnd.tab[1,3]), trendSums[2]/sum(trendSums[1:2])), 0)
  ave.less <- ifelse(mam.less < abs(rnd.tab[1,3]), abs(rnd.tab[1,3]) - mam.less, 0)
  rnd.tab[1,2] <- rnd.tab[1,2] - mam.less
  rnd.tab[1,1] <- ifelse(rnd.tab[1,3] > 0, rnd.tab[1,1], rnd.tab[1,1] - ave.less)
  rnd.tab[1,3] <- ifelse(rnd.tab[1,3] < 0, 0, rnd.tab[1,3])
  
  # AMR
  rnd.tab[2,1] <- rbinom(1,contSums[2],trendSums[1]/sum(trendSums))
  rnd.tab[2,2] <- rbinom(1,contSums[2],trendSums[2]/sum(trendSums))
  rnd.tab[2,3] <- contSums[2] - sum(rnd.tab[2,1], rnd.tab[2,2])
  mam.less <- ifelse(rnd.tab[2,3] < 0, rbinom(1, abs(rnd.tab[2,3]), trendSums[2]/sum(trendSums[1:2])), 0)
  ave.less <- ifelse(mam.less < abs(rnd.tab[2,3]), abs(rnd.tab[2,3]) - mam.less, 0)
  rnd.tab[2,2] <- rnd.tab[2,2] - mam.less
  rnd.tab[2,1] <- ifelse(rnd.tab[2,3] > 0, rnd.tab[2,1], rnd.tab[2,1] - ave.less)
  rnd.tab[2,3] <- ifelse(rnd.tab[2,3] < 0, 0, rnd.tab[2,3])
  
  # ASI
  rnd.tab[3,1] <- rbinom(1,contSums[3],trendSums[1]/sum(trendSums))
  rnd.tab[3,2] <- rbinom(1,contSums[3],trendSums[2]/sum(trendSums))
  rnd.tab[3,3] <- contSums[3] - sum(rnd.tab[3,1], rnd.tab[3,2])
  mam.less <- ifelse(rnd.tab[3,3] < 0, rbinom(1, abs(rnd.tab[3,3]), trendSums[2]/sum(trendSums[1:2])), 0)
  ave.less <- ifelse(mam.less < abs(rnd.tab[3,3]), abs(rnd.tab[3,3]) - mam.less, 0)
  rnd.tab[3,2] <- rnd.tab[3,2] - mam.less
  rnd.tab[3,1] <- ifelse(rnd.tab[3,3] > 0, rnd.tab[3,1], rnd.tab[3,1] - ave.less)
  rnd.tab[3,3] <- ifelse(rnd.tab[3,3] < 0, 0, rnd.tab[3,3])
  
  # EUR
  rnd.tab[4,1] <- rbinom(1,contSums[4],trendSums[1]/sum(trendSums))
  rnd.tab[4,2] <- rbinom(1,contSums[4],trendSums[2]/sum(trendSums))
  rnd.tab[4,3] <- contSums[4] - sum(rnd.tab[4,1], rnd.tab[4,2])
  mam.less <- ifelse(rnd.tab[4,3] < 0, rbinom(1, abs(rnd.tab[4,3]), trendSums[2]/sum(trendSums[1:2])), 0)
  ave.less <- ifelse(mam.less < abs(rnd.tab[5,3]), abs(rnd.tab[4,3]) - mam.less, 0)
  rnd.tab[4,2] <- rnd.tab[4,2] - mam.less
  rnd.tab[4,1] <- ifelse(rnd.tab[4,3] > 0, rnd.tab[4,1], rnd.tab[4,1] - ave.less)
  rnd.tab[4,3] <- ifelse(rnd.tab[4,3] < 0, 0, rnd.tab[4,3])
  
  # OCE
  rnd.tab[5,1] <- rbinom(1,contSums[5],trendSums[1]/sum(trendSums))
  rnd.tab[5,2] <- rbinom(1,contSums[5],trendSums[2]/sum(trendSums))
  rnd.tab[5,3] <- contSums[5] - sum(rnd.tab[5,1], rnd.tab[5,2])
  mam.less <- ifelse(rnd.tab[5,3] < 0, rbinom(1, abs(rnd.tab[5,3]), trendSums[2]/sum(trendSums[1:2])), 0)
  ave.less <- ifelse(mam.less < abs(rnd.tab[5,3]), abs(rnd.tab[5,3]) - mam.less, 0)
  rnd.tab[5,2] <- rnd.tab[5,2] - mam.less
  rnd.tab[5,1] <- ifelse(rnd.tab[5,3] > 0, rnd.tab[5,1], rnd.tab[5,1] - ave.less)
  rnd.tab[5,3] <- ifelse(rnd.tab[5,3] < 0, 0, rnd.tab[5,3])
  
  chi.tab.rnd <- (rnd.tab - exp.tab)^2/exp.tab
  inf.sub <- which(is.infinite(chi.tab.rnd)==T, arr.ind=T)
  if (length(inf.sub) > 0) {
    chi.tab.rnd[inf.sub] <- 0
    chi.tab.adj[inf.sub] <- 0
  } # end if  
  
  chi.rnd[c] <- sum(chi.tab.rnd, na.rm=T)
  comp.vec <- ifelse(chi.rnd[c] >= sum(chi.tab.adj, na.rm=T), 1, 0)
  
  if (c %% itdiv==0) print(c)
  
} # end c
pr <- sum(comp.vec)/iter
pr # low type I error for hypothesis that distribution of species by major trend category differs by continent



## spatial congruence
## generalised linear mixed-effects models
# model set
m1 <- "N_Species_In_Country ~ ClassComb + (1|Continent)"
m2 <- "N_Species_In_Country ~ 1 + (1|Continent)"

## make model vector
mod.vec <- as.character(mget(paste("m",rep(1:2,1),sep="")))

## define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glmer(as.formula(mod.vec[i]),data=spp.dat, family=binomial(link = "logit"), na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  LL.vec[i] <- as.numeric(logLik(fit))
  k.vec[i] <- attr(logLik(fit),"df")
  AICc.vec[i] <- r.squared(fit)$AIC
  BIC.vec[i] <- BIC(fit)
  
  Rm[i] <- 100*r.squared(fit)$Marginal # marginal R-squared
  Rc[i] <- 100*r.squared(fit)$Conditional # conditional R-squared
  
  print(i)
}
dAICc <- delta.IC(AICc.vec)
wAICc <- weight.IC(dAICc)
dBIC <- delta.IC(BIC.vec)
wBIC <- weight.IC(dBIC)

sumtable <- data.frame(mod.num,k.vec,LL.vec,AICc.vec,round(dAICc,3),round(wAICc,4),round(Rm,4),round(Rc,4))
colnames(sumtable) <- c("model","k","LL","AICc","dAICc","wAICc","Rm","Rc")
row.names(sumtable) <- as.character(mod.vec)
summary.table <- sumtable[order(sumtable[,6],decreasing=T),]
summary.table

## evidence for variation in proportion of species in country by class and continent

plot_model(mod.list[[1]], show.values=T, type="re", vline.color = "purple")
plot_model(mod.list[[1]], show.values=T, type="est", vline.color = "purple")

probs <- plogis(predict(mod.list[[1]]))
length(probs)
dim(spp.dat)
spp.dat$Classification_Class
spp.dat$ClassComb
spp.dat$N_Species_In_Country
missdatsub <- which(is.na(spp.dat$N_Species_In_Country))


# predicted probabilities by category
probs.dat <- data.frame("cont"=spp.dat$Continent[-missdatsub], "Class"=spp.dat$Classification_Class[-missdatsub],
                        "ClassComb"=spp.dat$ClassComb[-missdatsub], "probs"=probs)
head(probs.dat)

# by continent
EUR.mn <- mean(probs.dat$probs[probs.dat$cont=="Europe"], na.rm=T)
EUR.lo <- as.numeric(quantile(probs.dat$probs[probs.dat$cont=="Europe"], probs=0.025, na.rm = T))
EUR.up <- as.numeric(quantile(probs.dat$probs[probs.dat$cont=="Europe"], probs=0.975, na.rm = T))
EUR.sd <- sd(probs.dat$probs[probs.dat$cont=="Europe"], na.rm=T)
print(c(EUR.lo, EUR.mn, EUR.up))

AMR.mn <- mean(probs.dat$probs[probs.dat$cont=="Americas"], na.rm=T)
AMR.lo <- as.numeric(quantile(probs.dat$probs[probs.dat$cont=="Americas"], probs=0.025, na.rm = T))
AMR.up <- as.numeric(quantile(probs.dat$probs[probs.dat$cont=="Americas"], probs=0.975, na.rm = T))
AMR.sd <- sd(probs.dat$probs[probs.dat$cont=="Americas"], na.rm=T)
print(c(AMR.lo, AMR.mn, AMR.up))

OCN.mn <- mean(probs.dat$probs[probs.dat$cont=="Oceania"], na.rm=T)
OCN.lo <- as.numeric(quantile(probs.dat$probs[probs.dat$cont=="Oceania"], probs=0.025, na.rm = T))
OCN.up <- as.numeric(quantile(probs.dat$probs[probs.dat$cont=="Oceania"], probs=0.975, na.rm = T))
OCN.sd <- sd(probs.dat$probs[probs.dat$cont=="Oceania"], na.rm=T)
print(c(OCN.lo, OCN.mn, OCN.up))

ASIA.mn <- mean(probs.dat$probs[probs.dat$cont=="Asia"], na.rm=T)
ASIA.lo <- as.numeric(quantile(probs.dat$probs[probs.dat$cont=="Asia"], probs=0.025, na.rm = T))
ASIA.up <- as.numeric(quantile(probs.dat$probs[probs.dat$cont=="Asia"], probs=0.975, na.rm = T))
ASIA.sd <- sd(probs.dat$probs[probs.dat$cont=="Asia"], na.rm=T)
print(c(ASIA.lo, ASIA.mn, ASIA.up))

AFR.mn <- mean(probs.dat$probs[probs.dat$cont=="Africa"], na.rm=T)
AFR.lo <- as.numeric(quantile(probs.dat$probs[probs.dat$cont=="Africa"], probs=0.025, na.rm = T))
AFR.up <- as.numeric(quantile(probs.dat$probs[probs.dat$cont=="Africa"], probs=0.975, na.rm = T))
AFR.sd <- sd(probs.dat$probs[probs.dat$cont=="Africa"], na.rm=T) 
print(c(AFR.lo, AFR.mn, AFR.up))

prob.mns <- c(EUR.mn, AMR.mn, OCN.mn, ASIA.mn, AFR.mn)
prob.los <- c(EUR.lo, AMR.lo, OCN.lo, ASIA.lo, AFR.lo)
prob.ups <- c(EUR.up, AMR.up, OCN.up, ASIA.up, AFR.up)
prob.sds <- c(EUR.sd, AMR.sd, OCN.sd, ASIA.sd, AFR.sd)

prob.cont.summ <- data.frame("mean"=prob.mns, "sd"=prob.sds, "lo"=prob.los, "up"=prob.ups)
rownames(prob.cont.summ) <- c("EUR", "AMR", "OCN", "ASIA", "AFR")
prob.cont.summ

prob.cont.plot <- ggplot(prob.cont.summ) +
  geom_bar(aes(x=reorder(row.names(prob.cont.summ), mean), y=mean), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar( aes(x=row.names(prob.cont.summ), ymin=mean-sd, ymax=mean+sd),
                 linewidth=0.4, colour="black", alpha=0.9, size=0.3)
prob.cont.plot + coord_flip() +
  xlab("continent") + ylab("probability")

## probability of species chosen being in country increasing EUR < AMR < OCN < ASIA < AFR

# by Class
MAM.mn <- mean(probs.dat$probs[probs.dat$ClassComb=="Mammalia"], na.rm=T)
MAM.lo <- as.numeric(quantile(probs.dat$probs[probs.dat$ClassComb=="Mammalia"], probs=0.025, na.rm = T))
MAM.up <- as.numeric(quantile(probs.dat$probs[probs.dat$ClassComb=="Mammalia"], probs=0.975, na.rm = T))
MAM.sd <- sd(probs.dat$probs[probs.dat$ClassComb=="Mammalia"], na.rm=T) 
print(c(MAM.lo, MAM.mn, MAM.up))

AVE.mn <- mean(probs.dat$probs[probs.dat$ClassComb=="Aves"], na.rm=T)
AVE.lo <- as.numeric(quantile(probs.dat$probs[probs.dat$ClassComb=="Aves"], probs=0.025, na.rm = T))
AVE.up <- as.numeric(quantile(probs.dat$probs[probs.dat$ClassComb=="Aves"], probs=0.975, na.rm = T))
AVE.sd <- sd(probs.dat$probs[probs.dat$ClassComb=="Aves"], na.rm=T) 
print(c(AVE.lo, AVE.mn, AVE.up))

OTH.mn <- mean(probs.dat$probs[probs.dat$ClassComb=="other"], na.rm=T)
OTH.lo <- as.numeric(quantile(probs.dat$probs[probs.dat$ClassComb=="other"], probs=0.025, na.rm = T))
OTH.up <- as.numeric(quantile(probs.dat$probs[probs.dat$ClassComb=="other"], probs=0.975, na.rm = T))
OTH.sd <- sd(probs.dat$probs[probs.dat$ClassComb=="other"], na.rm=T) 
print(c(OTH.lo, OTH.mn, OTH.up))

probcl.mns <- c(MAM.mn, AVE.mn, OTH.mn)
probcl.los <- c(MAM.lo, AVE.lo, OTH.lo)
probcl.ups <- c(MAM.up, AVE.up, OTH.up)
probcl.sds <- c(MAM.sd, AVE.sd, OTH.sd)

prob.class.summ <- data.frame("mean"=probcl.mns, "sd"=probcl.sds, "lo"=probcl.los, "up"=probcl.ups)
rownames(prob.class.summ) <- c("MAM", "AVE", "OTH")
prob.class.summ

prob.cont.plot <- ggplot(prob.class.summ) +
  geom_bar(aes(x=reorder(row.names(prob.class.summ), mean), y=mean), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar( aes(x=row.names(prob.class.summ), ymin=mean-sd, ymax=mean+sd),
                 linewidth=0.4, colour="black", alpha=0.9, size=0.3)
prob.cont.plot + coord_flip() +
  xlab("Class") + ylab("probability")

## lower representation of mammals within country than Aves or 'other'




