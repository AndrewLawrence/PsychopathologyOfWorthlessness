library(tidyverse)
library(bootnet)
library(NetworkComparisonTest)

# Script Options ----------------------------------------------------------

# Number of cores to use:
core_opt <- 3
# Number of iterations (10 to test, 10000 to match paper):
niter_opt <- 10 #10000
# Tuning parameter:
gamma_opt <- 0.25

# Read in data ------------------------------------------------------------

study1 <- read.csv("input/study1.csv")
study2 <- read.csv("input/study2.csv")

# Threshold ---------------------------------------------------------------

study1 <- study1 %>% mutate_all(function(x) {as.numeric(x > 2)})
study2 <- study2 %>% mutate_all(function(x) {as.numeric(x > 2)})

# Bootstrap individual networks -------------------------------------------

model_file <- "model_objects/s1_bootstrap_D14_Ising_k10000_0.25.RData"
if ( ! file.exists(model_file) ) {

  bootD14_s1 <- bootnet(data = study1,
                     default = "IsingFit",
                     tuning = gamma_opt,
                     nBoots = niter_opt,
                     nCores = core_opt,
                     type = "nonparametric",
                     statistics = c("edge", "strength"),
                     weighted = TRUE,
                     directed = FALSE,
                     memorysaver = TRUE)

  save(bootD14_s1,
       file = model_file)

  # clear memory.
  rm(bootD14_s1)
  gc()
}

model_file <- "model_objects/s2_bootstrap_D14_Ising_k10000_0.25.RData"
if ( ! file.exists(model_file) ) {

  bootD14_s2 <- bootnet(data = study2,
                     default = "IsingFit",
                     tuning = gamma_opt,
                     nBoots = niter_opt,
                     nCores = core_opt,
                     type = "nonparametric",
                     statistics = c("edge", "strength"),
                     weighted = TRUE,
                     directed = FALSE,
                     memorysaver = TRUE)

  save(bootD14_s2,
       file = model_file)

  # clear memory.
  rm(bootD14_s2)
  gc()
}


# Permutation test --------------------------------------------------------

model_file <- "model_objects/nctD14_s1s2_bootstrap_SCL90_Ising_k10000_0.25.RData"
if ( ! file.exists(model_file) ) {

  nctD14_s1s2 <-
    NetworkComparisonTest::NCT(
      data1 = study1,
      data2 = study2,
      gamma = gamma_opt,
      it = niter_opt,
      binary.data = TRUE,
      paired = FALSE,
      weighted = TRUE,
      AND = TRUE,
      abs = TRUE,
      test.edges = TRUE,
      verbose = TRUE,
      progressbar = TRUE)

  save(nctD14_s1s2,
       file = model_file)
}


# Also write out data ----------------------------------------------------

save(study1, study2, file = "model_objects/s1s2_input.RData")


# Done.

