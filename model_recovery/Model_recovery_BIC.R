## Model comparison

library("pracma")
library("rstan")
library("truncnorm")
library("MASS")
library("TruncatedNormal")
library("loo")
rm(list=ls()) # clear all
setwd('./model_recovery')

M1_counts <- vector()
M3_counts <- vector()
M4_counts <- vector()
M5_counts <- vector()

## Generated Model 1
load("fit_M1_GM1.RData")
lik_M1_GM1    <- loo::extract_log_lik(stanfit = fit_simulated_M1_GM1, parameter_name = 'log_lik')
lik_M1_GM1 <- apply(lik_M1_GM1,2,max)
rm(fit_simulated_M1_GM1)

load("fit_M3_GM1.RData")
lik_M3_GM1    <- loo::extract_log_lik(stanfit = fit_simulated_M3_GM1, parameter_name = 'log_lik')
lik_M3_GM1 <- apply(lik_M3_GM1,2,max)
rm(fit_simulated_M3_GM1)

load("fit_M4_GM1.RData")
lik_M4_GM1    <- loo::extract_log_lik(stanfit = fit_simulated_M4_GM1, parameter_name = 'log_lik')
lik_M4_GM1 <- apply(lik_M4_GM1,2,max)
rm(fit_simulated_M4_GM1)

load("fit_M5_GM1.RData")
lik_M5_GM1    <- loo::extract_log_lik(stanfit = fit_simulated_M5_GM1, parameter_name = 'log_lik')
lik_M5_GM1 <- apply(lik_M5_GM1,2,max)
rm(fit_simulated_M5_GM1)


M1_GM1_BIC <- 2*log(250)-2*lik_M1_GM1
M3_GM1_BIC <- 3*log(250)-2*lik_M3_GM1
M4_GM1_BIC <- 3*log(250)-2*lik_M4_GM1
M5_GM1_BIC <- 3*log(250)-2*lik_M5_GM1

M1_count <- 0
M3_count <- 0
M4_count <- 0
M5_count <- 0
for (i in 1:100) {
  if (M1_GM1_BIC[i]<M3_GM1_BIC[i] & M1_GM1_BIC[i]<M4_GM1_BIC[i] & M1_GM1_BIC[i]<M5_GM1_BIC[i]) {
    M1_count <- M1_count+1
  }
  
  if (M3_GM1_BIC[i]<M1_GM1_BIC[i] & M3_GM1_BIC[i]<M4_GM1_BIC[i] & M3_GM1_BIC[i]<M5_GM1_BIC[i]) {
    M3_count <- M3_count+1
  }
  
  if (M4_GM1_BIC[i]<M3_GM1_BIC[i] & M4_GM1_BIC[i]<M1_GM1_BIC[i] & M4_GM1_BIC[i]<M5_GM1_BIC[i]) {
    M4_count <- M4_count+1
  }
  
  if (M5_GM1_BIC[i]<M3_GM1_BIC[i] & M5_GM1_BIC[i]<M4_GM1_BIC[i] & M5_GM1_BIC[i]<M1_GM1_BIC[i]) {
    M5_count <- M5_count+1
  }
  
}

M1_counts <- append(M1_counts,M1_count)
M3_counts <- append(M3_counts,M3_count)
M4_counts <- append(M4_counts,M4_count)
M5_counts <- append(M5_counts,M5_count)


## Generated Model 3
load("fit_M1_GM3.RData")
lik_M1_GM3    <- loo::extract_log_lik(stanfit = fit_simulated_M1_GM3, parameter_name = 'log_lik')
lik_M1_GM3 <- apply(lik_M1_GM3,2,max)
rm(fit_simulated_M1_GM3)

load("fit_M3_GM3.RData")
lik_M3_GM3    <- loo::extract_log_lik(stanfit = fit_simulated_M3_GM3, parameter_name = 'log_lik')
lik_M3_GM3 <- apply(lik_M3_GM3,2,max)
rm(fit_simulated_M3_GM3)

load("fit_M4_GM3.RData")
lik_M4_GM3    <- loo::extract_log_lik(stanfit = fit_simulated_M4_GM3, parameter_name = 'log_lik')
lik_M4_GM3 <- apply(lik_M4_GM3,2,max)
rm(fit_simulated_M4_GM3)

load("fit_M5_GM3.RData")
lik_M5_GM3    <- loo::extract_log_lik(stanfit = fit_simulated_M5_GM3, parameter_name = 'log_lik')
lik_M5_GM3 <- apply(lik_M5_GM3,2,max)
rm(fit_simulated_M5_GM3)




M1_GM3_BIC <- 2*log(250)-2*lik_M1_GM3
M3_GM3_BIC <- 3*log(250)-2*lik_M3_GM3
M4_GM3_BIC <- 3*log(250)-2*lik_M4_GM3
M5_GM3_BIC <- 3*log(250)-2*lik_M5_GM3

M1_count <- 0
M3_count <- 0
M4_count <- 0
M5_count <- 0
for (i in 1:100) {
  if (M1_GM3_BIC[i]<M3_GM3_BIC[i] & M1_GM3_BIC[i]<M4_GM3_BIC[i] & M1_GM3_BIC[i]<M5_GM3_BIC[i]) {
    M1_count <- M1_count+1
  }
  
  if (M3_GM3_BIC[i]<M1_GM3_BIC[i] & M3_GM3_BIC[i]<M4_GM3_BIC[i] & M3_GM3_BIC[i]<M5_GM3_BIC[i]) {
    M3_count <- M3_count+1
  }
  
  if (M4_GM3_BIC[i]<M3_GM3_BIC[i] & M4_GM3_BIC[i]<M1_GM3_BIC[i] & M4_GM3_BIC[i]<M5_GM3_BIC[i]) {
    M4_count <- M4_count+1
  }
  
  if (M5_GM3_BIC[i]<M3_GM3_BIC[i] & M5_GM3_BIC[i]<M4_GM3_BIC[i] & M5_GM3_BIC[i]<M1_GM3_BIC[i]) {
    M5_count <- M5_count+1
  }
  
}

M1_counts <- append(M1_counts,M1_count)
M3_counts <- append(M3_counts,M3_count)
M4_counts <- append(M4_counts,M4_count)
M5_counts <- append(M5_counts,M5_count)

# Data generated from model 3 AIC and BIC

2*2-2*mean(lik_M1_GM3)
2*3-2*mean(lik_M3_GM3)
2*3-2*mean(lik_M4_GM3)
2*3-2*mean(lik_M5_GM3)


# BIC
2*log(250)-2*mean(lik_M1_GM3)
3*log(250)-2*mean(lik_M3_GM3)
3*log(250)-2*mean(lik_M4_GM3)
3*log(250)-2*mean(lik_M5_GM3)

## Generated Model 4
load("fit_M1_GM4.RData")
lik_M1_GM4    <- loo::extract_log_lik(stanfit = fit_simulated_M1_GM4, parameter_name = 'log_lik')
lik_M1_GM4 <- apply(lik_M1_GM4,2,max)
rm(fit_simulated_M1_GM4)

load("fit_M3_GM4.RData")
lik_M3_GM4    <- loo::extract_log_lik(stanfit = fit_simulated_M3_GM4, parameter_name = 'log_lik')
lik_M3_GM4 <- apply(lik_M3_GM4,2,max)
rm(fit_simulated_M3_GM4)

load("fit_M4_GM4.RData")
lik_M4_GM4    <- loo::extract_log_lik(stanfit = fit_simulated_M4_GM4, parameter_name = 'log_lik')
lik_M4_GM4 <- apply(lik_M4_GM4,2,max)
rm(fit_simulated_M4_GM4)

load("fit_M5_GM4.RData")
lik_M5_GM4    <- loo::extract_log_lik(stanfit = fit_simulated_M5_GM4, parameter_name = 'log_lik')
lik_M5_GM4 <- apply(lik_M5_GM4,2,max)
rm(fit_simulated_M5_GM4)


M1_GM4_BIC <- 2*log(250)-2*lik_M1_GM4
M3_GM4_BIC <- 3*log(250)-2*lik_M3_GM4
M4_GM4_BIC <- 3*log(250)-2*lik_M4_GM4
M5_GM4_BIC <- 3*log(250)-2*lik_M5_GM4

M1_count <- 0
M3_count <- 0
M4_count <- 0
M5_count <- 0
for (i in 1:100) {
  if (M1_GM4_BIC[i]<M3_GM4_BIC[i] & M1_GM4_BIC[i]<M4_GM4_BIC[i] & M1_GM4_BIC[i]<M5_GM4_BIC[i]) {
    M1_count <- M1_count+1
  }
  
  if (M3_GM4_BIC[i]<M1_GM4_BIC[i] & M3_GM4_BIC[i]<M4_GM4_BIC[i] & M3_GM4_BIC[i]<M5_GM4_BIC[i]) {
    M3_count <- M3_count+1
  }
  
  if (M4_GM4_BIC[i]<M3_GM4_BIC[i] & M4_GM4_BIC[i]<M1_GM4_BIC[i] & M4_GM4_BIC[i]<M5_GM4_BIC[i]) {
    M4_count <- M4_count+1
  }
  
  if (M5_GM4_BIC[i]<M3_GM4_BIC[i] & M5_GM4_BIC[i]<M4_GM4_BIC[i] & M5_GM4_BIC[i]<M1_GM4_BIC[i]) {
    M5_count <- M5_count+1
  }
  
}

M1_counts <- append(M1_counts,M1_count)
M3_counts <- append(M3_counts,M3_count)
M4_counts <- append(M4_counts,M4_count)
M5_counts <- append(M5_counts,M5_count)

# Data generated from model 4 AIC and BIC

2*2-2*mean(lik_M1_GM4)
2*3-2*mean(lik_M3_GM4)
2*3-2*mean(lik_M4_GM4)
2*3-2*mean(lik_M5_GM4)


# BIC
2*log(250)-2*mean(lik_M1_GM4)
3*log(250)-2*mean(lik_M3_GM4)
3*log(250)-2*mean(lik_M4_GM4)
3*log(250)-2*mean(lik_M5_GM4)

## Generated Model 5
load("fit_M1_GM5.RData")
lik_M1_GM5    <- loo::extract_log_lik(stanfit = fit_simulated_M1_GM5, parameter_name = 'log_lik')
lik_M1_GM5 <- apply(lik_M1_GM5,2,max)
rm(fit_simulated_M1_GM5)

load("fit_M3_GM5.RData")
lik_M3_GM5    <- loo::extract_log_lik(stanfit = fit_simulated_M3_GM5, parameter_name = 'log_lik')
lik_M3_GM5 <- apply(lik_M3_GM5,2,max)
rm(fit_simulated_M3_GM5)

load("fit_M4_GM5.RData")
lik_M4_GM5    <- loo::extract_log_lik(stanfit = fit_simulated_M4_GM5, parameter_name = 'log_lik')
lik_M4_GM5 <- apply(lik_M4_GM5,2,max)
rm(fit_simulated_M4_GM5)

load("fit_M5_GM5.RData")
lik_M5_GM5    <- loo::extract_log_lik(stanfit = fit_simulated_M5_GM5, parameter_name = 'log_lik')
lik_M5_GM5 <- apply(lik_M5_GM5,2,max)
rm(fit_simulated_M5_GM5)

M1_GM5_BIC <- 2*log(250)-2*lik_M1_GM5
M3_GM5_BIC <- 3*log(250)-2*lik_M3_GM5
M4_GM5_BIC <- 3*log(250)-2*lik_M4_GM5
M5_GM5_BIC <- 3*log(250)-2*lik_M5_GM5

M1_count <- 0
M3_count <- 0
M4_count <- 0
M5_count <- 0
for (i in 1:100) {
  if (M1_GM5_BIC[i]<M3_GM5_BIC[i] & M1_GM5_BIC[i]<M4_GM5_BIC[i] & M1_GM5_BIC[i]<M5_GM5_BIC[i]) {
    M1_count <- M1_count+1
  }
  
  if (M3_GM5_BIC[i]<M1_GM5_BIC[i] & M3_GM5_BIC[i]<M4_GM5_BIC[i] & M3_GM5_BIC[i]<M5_GM5_BIC[i]) {
    M3_count <- M3_count+1
  }
  
  if (M4_GM5_BIC[i]<M3_GM5_BIC[i] & M4_GM5_BIC[i]<M1_GM5_BIC[i] & M4_GM5_BIC[i]<M5_GM5_BIC[i]) {
    M4_count <- M4_count+1
  }
  
  if (M5_GM5_BIC[i]<M3_GM5_BIC[i] & M5_GM5_BIC[i]<M4_GM5_BIC[i] & M5_GM5_BIC[i]<M1_GM5_BIC[i]) {
    M5_count <- M5_count+1
  }
  
}

M1_counts <- append(M1_counts,M1_count)
M3_counts <- append(M3_counts,M3_count)
M4_counts <- append(M4_counts,M4_count)
M5_counts <- append(M5_counts,M5_count)

# Data generated from model 5 AIC and BIC

2*2-2*mean(lik_M1_GM5)
2*3-2*mean(lik_M3_GM5)
2*3-2*mean(lik_M4_GM5)
2*3-2*mean(lik_M5_GM5)


# BIC
2*log(250)-2*mean(lik_M1_GM5)
3*log(250)-2*mean(lik_M3_GM5)
3*log(250)-2*mean(lik_M4_GM5)
3*log(250)-2*mean(lik_M5_GM5)

########################## Create matrix
# Rows generated (true) model, column percentage best fit model
confusion_matrix <- cbind(M1_counts,M3_counts,M4_counts,M5_counts)
write.csv(confusion_matrix, "BIC_confusion_matrix.csv", row.names = FALSE)

