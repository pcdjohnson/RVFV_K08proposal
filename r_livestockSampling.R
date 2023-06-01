rm(list = ls())
# script to estimate power to compare RVFV seroprevalence between predicted high and low risk areas
start.time <- Sys.time()

# load packages
library(GLMMmisc) # available via devtools::install_github("pcdjohnson/GLMMmisc")
library(lme4)
library(parallel)

# no of data sets to simulate per species (1000 takes ~ 4 min)
nsim <- 1000

# sampling / design choices
n.village <- c(32, 48)[1]
n.hh <- c(15, 10)[1]         # per village
n <- 3                       # no of animals/humans per household

n.village * n.hh * n

# effect size: OR representing difference in seroprevalence 
# between high and low risk villages
OR <- 3

# Livestock serology data:
# rvf_livestock_data_2719.csv
# Human serology data:
# rvf_human_data_2719.csv and humanhouseholdlink_nocoords.csv

# note that in the human data, because of the low number of positives
# and the low number of individuals sampled per HH (91/242 HH had only one
# indivuadal sampled, and 166 HH had 1-2 individuals), there wasn't a lot 
# of information on the HH random effect, so I didn't fit it (i.e. assumed 
# it was zero).

# no of positives:
# cattle 156/3582 (4.4%)
# goat    45/3303 (1.4%)
# sheep   67/2584 (2.6%)
# human   48/574  (8.4%)


# load parameter estimates
par.tab <- read.csv("parameter.estimates.csv", row.names = 1)
species <- colnames(par.tab)


# simulate RVFV serology data in each species

# function to simulate data and estimate p-value for null hypothesis 
# that high and low risk areas have the same seroprevalence
res.tab.fn <- function(rand.seed = NULL, ...) {
  sapply(species, function(sp) {
    
    # create template data set
    dat <- expand.grid(hh = 1:n.hh, village = 1:n.village, n = n)
    print(sum(dat$n))
    # allocate villages to high and low prevalence in 1:1 ratio 
    dat$risk.level <- dat$village %% 2 - 0.5
    # simulate seropositives
    set.seed(rand.seed)
    simdat <-
      sim.glmm(
        design.data = dat, 
        fixed.eff = 
          list(
            intercept = par.tab["(Intercept)", sp],
            risk.level = log(OR)),
        distribution = "binomial",
        rand.V = c(hh = par.tab["barcode_hh", sp], 
                   village = par.tab["village", sp]))
    form <- 
      if(par.tab["barcode_hh", sp] == 0) {
        cbind(response, n - response) ~ risk.level + (1 | village)
      } else {
        cbind(response, n - response) ~ risk.level + (1 | hh) + (1 | village)
      }
    
    
    fit <- glmer(form, family = binomial, data = simdat, control = glmerControl(optimizer = "bobyqa"))
    fit0 <- update(fit, ~ . - risk.level)
    
    ICC.est <- sum(unlist(VarCorr(fit))) / (sum(unlist(VarCorr(fit))) + pi^2 / 3)
    
    #coef(summary(fit))["risk.level", "Pr(>|z|)"] # Wald P not reliable - gives inflated type 1 error
    c(p = anova(fit, fit0)[2, "Pr(>Chisq)"], 
      ci = plogis(confint(fit, method = "Wald")["(Intercept)", ]),
      ICC.est = ICC.est)
  })
}


# repeat simulations many times and calculate p-value
init.rand.seed <- round((as.numeric(Sys.time()) - floor(as.numeric(Sys.time())))*10000000)  # 
#init.rand.seed <- 6569149
set.seed(init.rand.seed)
rand.seed.list <- sample(1e9, nsim)

sim.res <- mclapply(rand.seed.list, function(rand.seed) {
  res.tab.fn(rand.seed = rand.seed)  
}, mc.cores = detectCores())
print(Sys.time() - start.time)

stopifnot(all(sapply(sim.res, is.numeric)))

# estimate power
out.tab <-
  sapply(species, function(sp) {
    re.var <- sum(par.tab[c("barcode_hh", "village"), sp])
    ICC <- re.var/(re.var + pi^2 / 3)
    tab <- t(sapply(sim.res, function(x) x[, sp]))
    power <-mean(tab[, "p"] < 0.05)
    mean.moe <- 
      mean(apply(tab[, c("ci.2.5 %", "ci.97.5 %")], 1, diff))/2
    mean.icc <-mean(tab[, "ICC.est"])
    c(n.village = n.village, n.hh = n.hh, n.per.hh = n, OR = OR,
      nsim = nsim, power = power, mean.moe = mean.moe, 
      mean.icc = mean.icc, true.icc = ICC)
  })
round(out.tab, 3)
