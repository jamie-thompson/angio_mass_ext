library(CRABS)
library(ggplot2)
library(phytools)
library(TESS)

# TESS analyses
# First Smith and Brown tree

tree = read.tree("ultra_Angio_SB_prune.tre")

samplingFraction <- length(tree$tip.label)/290000
numExpectedRateChanges <- 25
numExpectedMassExtinctions <- 1
expectedSurvivalProbability <- 0.25
pMassExtinctionPriorShape2 <- 100
pMassExtinctionPriorShape1 <- - pMassExtinctionPriorShape2 *
  expectedSurvivalProbability /
  (expectedSurvivalProbability - 1)

set.seed(73612) # random

tess.analysis(tree, 
              empiricalHyperPriors = TRUE, 
              samplingProbability = samplingFraction,
              numExpectedRateChanges = numExpectedRateChanges, 
              numExpectedMassExtinctions = numExpectedMassExtinctions, 
              pMassExtinctionPriorShape1 = pMassExtinctionPriorShape1, 
              pMassExtinctionPriorShape2 = pMassExtinctionPriorShape2, 
              MAX_ITERATIONS = 25000000,
              MAX_TIME = Inf, MIN_ESS = 300,
              dir = "SB_TESS_1_OP")


output <- tess.process.output("SB_TESS_1_OP",
                              numExpectedRateChanges = numExpectedRateChanges,
                              numExpectedMassExtinctions = numExpectedMassExtinctions)

saveRDS(output,"SB_TESS_1.RDS")

set.seed(6879) # random

tess.analysis(tree, 
              empiricalHyperPriors = TRUE, 
              samplingProbability = samplingFraction,
              numExpectedRateChanges = numExpectedRateChanges, 
              numExpectedMassExtinctions = numExpectedMassExtinctions, 
              pMassExtinctionPriorShape1 = pMassExtinctionPriorShape1, 
              pMassExtinctionPriorShape2 = pMassExtinctionPriorShape2, 
              MAX_ITERATIONS = 25000000,
              MAX_TIME = Inf, MIN_ESS = 300,
              dir = "SB_TESS_2_OP")


output2 <- tess.process.output("SB_TESS_2_OP",
                               numExpectedRateChanges = numExpectedRateChanges,
                               numExpectedMassExtinctions = numExpectedMassExtinctions)

saveRDS(output2,"SB_TESS_2.RDS")

# Now Qian and Jin tree

tree = read.tree("Qian_Jin_angios.tre")
tree = force.ultrametric(tree,method="extend")

samplingFraction <- length(tree$tip.label)/290000
numExpectedRateChanges <- 25
numExpectedMassExtinctions <- 1
expectedSurvivalProbability <- 0.25
pMassExtinctionPriorShape2 <- 100
pMassExtinctionPriorShape1 <- - pMassExtinctionPriorShape2 *
  expectedSurvivalProbability /
  (expectedSurvivalProbability - 1)

set.seed(73612) # random

tess.analysis(tree, 
              empiricalHyperPriors = TRUE, 
              samplingProbability = samplingFraction,
              numExpectedRateChanges = numExpectedRateChanges, 
              numExpectedMassExtinctions = numExpectedMassExtinctions, 
              pMassExtinctionPriorShape1 = pMassExtinctionPriorShape1, 
              pMassExtinctionPriorShape2 = pMassExtinctionPriorShape2, 
              MAX_ITERATIONS = 25000000,
              MAX_TIME = Inf, MIN_ESS = 300,
              dir = "QJ_TESS_1_OP")


output <- tess.process.output("QJ_TESS_1_OP",
                              numExpectedRateChanges = numExpectedRateChanges,
                              numExpectedMassExtinctions = numExpectedMassExtinctions)

saveRDS(output,"QJ_TESS_1.RDS")

set.seed(6879) # random

tess.analysis(tree, 
              empiricalHyperPriors = TRUE, 
              samplingProbability = samplingFraction,
              numExpectedRateChanges = numExpectedRateChanges, 
              numExpectedMassExtinctions = numExpectedMassExtinctions, 
              pMassExtinctionPriorShape1 = pMassExtinctionPriorShape1, 
              pMassExtinctionPriorShape2 = pMassExtinctionPriorShape2, 
              MAX_ITERATIONS = 25000000,
              MAX_TIME = Inf, MIN_ESS = 300,
              dir = "QJ_TESS_2_OP")


output2 <- tess.process.output("QJ_TESS_2_OP",
                               numExpectedRateChanges = numExpectedRateChanges,
                               numExpectedMassExtinctions = numExpectedMassExtinctions)

saveRDS(output2,"QJ_TESS_2.RDS")

# CRABS, testing models of constant extinction - following tutorial available here: https://afmagee.github.io/CRABS/overview.html
# QJ - read tree and TESS data from 25 shift analysis

tree = read.tree("Qian_Jin_angios.tre")
tree = force.ultrametric(tree,method="extend")

qjrun2 = read_rds("QJ_TESS_1.RDS")
qjrun3 = read_rds("QJ_TESS_2.RDS")

spec = rbind(qjrun1[["speciation rates"]],qjrun2[["speciation rates"]])
ext = rbind(qjrun1[["extinction rates"]],qjrun2[["extinction rates"]])

meanspec <- colMeans(spec)
meanext <- colMeans(ext)

intervals = qjrun1$intervals

mean_qj = data.frame(lambda = c(meanspec[1], meanspec), mu = c(meanext[1], meanext), time = intervals)
mean_qj = mean_qj[order(mean_qj$time),]

# CRABS code

btimes <- sort( as.numeric( branching.times( tree ) ) )
max_t <- max(btimes)

# 4.3 in tutorial

max_t <- max(mean_qj[["time"]])
times_fine <- seq(0, max_t, length.out = 100)

times <- mean_qj[["time"]]

lambda <- approxfun(times, mean_qj[["lambda"]])
mu <- approxfun(times, mean_qj[["mu"]])


my_model <- create.model(lambda, 
                         mu, 
                         times_fine)
plot(my_model)


# 5.1.1 Constant model

mu_vals <- seq(0,3,0.05) # Mu values chosen as 0 - 3
mu1 <- list()
for (i in 1:length(mu_vals)) {
  mu1[[i]] = local({
    mu = mu_vals[i]
    function(t) { rep(mu,length(t)) }
  })
}

alt_models <- congruent.models(my_model, mus = mu1)

plot( alt_models )

a <- summarize.trends(alt_models, threshold = 0.02)
plot(a)


# SB- read tree and TESS data from 25 shift analysis

tree = read.tree("ultra_Angio_SB_prune.tre")

sbrun1 = read_rds("SB_TESS_1.RDS")
sbrun2 = read_rds("SB_TESS_2.RDS")

spec = rbind(sbrun1[["speciation rates"]],sbrun2[["speciation rates"]])
ext = rbind(sbrun1[["extinction rates"]],sbrun2[["extinction rates"]])

meanspec <- colMeans(spec)
meanext <- colMeans(ext)

intervals = sbrun1$intervals

mean_sb = data.frame(lambda = c(meanspec[1], meanspec), mu = c(meanext[1], meanext), time = intervals)
mean_sb = mean_sb[order(mean_sb$time),]

# CRABS code

btimes <- sort( as.numeric( branching.times( tree ) ) )
max_t <- max(btimes)

# 4.3 in tutorial

max_t <- max(mean_sb[["time"]])
times_fine <- seq(0, max_t, length.out = 100)

times <- mean_sb[["time"]]

lambda <- approxfun(times, mean_sb[["lambda"]])
mu <- approxfun(times, mean_sb[["mu"]])


my_model <- create.model(lambda, 
                         mu, 
                         times_fine)
plot(my_model)

# 5.1.1 Constant model

mu_vals <- seq(0,3,0.05) # Mu values chosen as 0 - 3
mu1 <- list()
for (i in 1:length(mu_vals)) {
  mu1[[i]] = local({
    mu = mu_vals[i]
    function(t) { rep(mu,length(t)) }
  })
}

alt_models <- congruent.models(my_model, mus = mu1)

plot( alt_models )

b <- summarize.trends(alt_models, threshold = 0.02)
plot(b)

# Sensitivity testing with Qian and Jin tree

tree = read.tree("Qian_Jin_angios.tre")
tree = force.ultrametric(tree,method="extend")

# 1 shift

samplingFraction <- length(tree$tip.label)/290000
numExpectedRateChanges <- 1
numExpectedMassExtinctions <- 1
expectedSurvivalProbability <- 0.25
pMassExtinctionPriorShape2 <- 100
pMassExtinctionPriorShape1 <- - pMassExtinctionPriorShape2 *
  expectedSurvivalProbability /
  (expectedSurvivalProbability - 1)

tess.analysis(tree, 
              empiricalHyperPriors = TRUE, 
              samplingProbability = samplingFraction,
              numExpectedRateChanges = numExpectedRateChanges, 
              numExpectedMassExtinctions = numExpectedMassExtinctions, 
              pMassExtinctionPriorShape1 = pMassExtinctionPriorShape1, 
              pMassExtinctionPriorShape2 = pMassExtinctionPriorShape2, 
              MAX_ITERATIONS = 25000000,
              MAX_TIME = Inf, MIN_ESS = 300,
              dir = "1_shift_Qian_Jin_300_ESS")


output <- tess.process.output("1_shift_Qian_Jin_300_ESS",
                              numExpectedRateChanges = numExpectedRateChanges,
                              numExpectedMassExtinctions = numExpectedMassExtinctions)

saveRDS(output,"1_shift_Qian_Jin_300_ESS.RDS")


layout.mat <- matrix(1:6,nrow=3,ncol=2,byrow=TRUE)
layout(layout.mat)

pdf("1_shift_Qian_Jin_300_ESS.pdf")

tess.plot.output(output,
                 fig.types = c("speciation rates",
                               "speciation shift times",
                               "extinction rates",
                               "extinction shift times",
                               "mass extinction Bayes factors",
                               "mass extinction times"),
                 las=2)

dev.off()


# 100 shifts

samplingFraction <- length(tree$tip.label)/290000
numExpectedRateChanges <- 100
numExpectedMassExtinctions <- 1
expectedSurvivalProbability <- 0.25
pMassExtinctionPriorShape2 <- 100
pMassExtinctionPriorShape1 <- - pMassExtinctionPriorShape2 *
  expectedSurvivalProbability /
  (expectedSurvivalProbability - 1)

tess.analysis(tree, 
              empiricalHyperPriors = TRUE, 
              samplingProbability = samplingFraction,
              numExpectedRateChanges = numExpectedRateChanges, 
              numExpectedMassExtinctions = numExpectedMassExtinctions, 
              pMassExtinctionPriorShape1 = pMassExtinctionPriorShape1, 
              pMassExtinctionPriorShape2 = pMassExtinctionPriorShape2, 
              MAX_ITERATIONS = 25000000,
              MAX_TIME = Inf, MIN_ESS = 300,
              dir = "100_shifts_Qian_Jin_300_ESS")


output <- tess.process.output("100_shifts_Qian_Jin_300_ESS",
                              numExpectedRateChanges = numExpectedRateChanges,
                              numExpectedMassExtinctions = numExpectedMassExtinctions)

saveRDS(output,"100_shifts_Qian_Jin_300_ESS.RDS")

# 250 shifts

samplingFraction <- length(tree$tip.label)/290000
numExpectedRateChanges <- 250
numExpectedMassExtinctions <- 1
expectedSurvivalProbability <- 0.25
pMassExtinctionPriorShape2 <- 100
pMassExtinctionPriorShape1 <- - pMassExtinctionPriorShape2 *
  expectedSurvivalProbability /
  (expectedSurvivalProbability - 1)

tess.analysis(tree, 
              empiricalHyperPriors = TRUE, 
              samplingProbability = samplingFraction,
              numExpectedRateChanges = numExpectedRateChanges, 
              numExpectedMassExtinctions = numExpectedMassExtinctions, 
              pMassExtinctionPriorShape1 = pMassExtinctionPriorShape1, 
              pMassExtinctionPriorShape2 = pMassExtinctionPriorShape2, 
              MAX_ITERATIONS = 25000000,
              MAX_TIME = Inf, MIN_ESS = 300,
              dir = "250_shifts_Qian_Jin_300_ESS")


output <- tess.process.output("250_shifts_Qian_Jin_300_ESS",
                              numExpectedRateChanges = numExpectedRateChanges,
                              numExpectedMassExtinctions = numExpectedMassExtinctions)

saveRDS(output,"250_shifts_Qian_Jin_300_ESS.RDS")

# 500 shifts

samplingFraction <- length(tree$tip.label)/290000
numExpectedRateChanges <- 500
numExpectedMassExtinctions <- 1
expectedSurvivalProbability <- 0.25
pMassExtinctionPriorShape2 <- 100
pMassExtinctionPriorShape1 <- - pMassExtinctionPriorShape2 *
  expectedSurvivalProbability /
  (expectedSurvivalProbability - 1)

tess.analysis(tree, 
              empiricalHyperPriors = TRUE, 
              samplingProbability = samplingFraction,
              numExpectedRateChanges = numExpectedRateChanges, 
              numExpectedMassExtinctions = numExpectedMassExtinctions, 
              pMassExtinctionPriorShape1 = pMassExtinctionPriorShape1, 
              pMassExtinctionPriorShape2 = pMassExtinctionPriorShape2, 
              MAX_ITERATIONS = 25000000,
              MAX_TIME = Inf, MIN_ESS = 300,
              dir = "500_shifts_Qian_Jin_300_ESS")


output <- tess.process.output("500_shifts_Qian_Jin_300_ESS",
                              numExpectedRateChanges = numExpectedRateChanges,
                              numExpectedMassExtinctions = numExpectedMassExtinctions)

saveRDS(output,"500_shifts_Qian_Jin_300_ESS.RDS")

# 1000 shifts

samplingFraction <- length(tree$tip.label)/290000
numExpectedRateChanges <- 1000
numExpectedMassExtinctions <- 1
expectedSurvivalProbability <- 0.25
pMassExtinctionPriorShape2 <- 100
pMassExtinctionPriorShape1 <- - pMassExtinctionPriorShape2 *
  expectedSurvivalProbability /
  (expectedSurvivalProbability - 1)

tess.analysis(tree, 
              empiricalHyperPriors = TRUE, 
              samplingProbability = samplingFraction,
              numExpectedRateChanges = numExpectedRateChanges, 
              numExpectedMassExtinctions = numExpectedMassExtinctions, 
              pMassExtinctionPriorShape1 = pMassExtinctionPriorShape1, 
              pMassExtinctionPriorShape2 = pMassExtinctionPriorShape2, 
              MAX_ITERATIONS = 25000000,
              MAX_TIME = Inf, MIN_ESS = 300,
              dir = "1000_shifts_Qian_Jin_300_ESS")


output <- tess.process.output("1000_shifts_Qian_Jin_300_ESS",
                              numExpectedRateChanges = numExpectedRateChanges,
                              numExpectedMassExtinctions = numExpectedMassExtinctions)

saveRDS(output,"1000_shifts_Qian_Jin_300_ESS.RDS")


layout.mat <- matrix(1:6,nrow=3,ncol=2,byrow=TRUE)
layout(layout.mat)

pdf("1000_shifts_Qian_Jin_300_ESS.pdf")

tess.plot.output(output,
                 fig.types = c("speciation rates",
                               "speciation shift times",
                               "extinction rates",
                               "extinction shift times",
                               "mass extinction Bayes factors",
                               "mass extinction times"),
                 las=2)

dev.off()