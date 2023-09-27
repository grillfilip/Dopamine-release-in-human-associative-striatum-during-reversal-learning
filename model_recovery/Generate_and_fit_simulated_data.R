

library("pracma")
library("rstan")
library("ggplot2")
library("ggpubr")
library("truncnorm")
library("MASS")
library("TruncatedNormal")
setwd('./model_recovery')


################################################################################## Generate Data M1
rm(list=ls()) # clear all
set.seed(999)

N <- 100
fit_params <- read.csv('fit_params_M1.csv') # Read fit real data

fit_params_cov <- cov(fit_params) # Get covariance from parameters

Alpha <- fit_params$Alpha
Beta <- fit_params$Beta
params <- mvrandn(c(min(Alpha),min(Beta)),c(max(Alpha),max(Beta)), fit_params_cov, N, c(mean(Alpha), mean(Beta))) # Generate new parameters bound by min, max, mean, and covariance of parameters

Alpha_sim <- params[1,]
Beta_sim <- params[2,]

outcome <- vector()
choice <- vector()
PEs <- vector()
evA <- 0 # Initial expected value
evB <- 0


all_evA <- array(0, dim = c(251,N)) # For saving data in loop
all_evB <- array(0, dim = c(251,N))
all_choice <- array(0, dim = c(250,N))
all_outcome <- array(0, dim = c(250,N))
all_PEs <- array(0, dim = c(250,N))


# Task structure
for (i in 1:N){
  for (t in 1:250){
    
    probA <- exp(Beta_sim[i]*evA[t])/(exp(Beta_sim[i]*evA[t])+exp(Beta_sim[i]*evB[t])) # Prob. choose A based on Beta_sim and expected value
    
    # if choice = 1 choice was A, else choice was B
    choice[t] <- rbinom(1, size = 1, prob = probA)
    
    
    if (t <= 150){
      if (choice[t] == 1){ # if choice is A
        outcome[t] <- rbinom(1, size = 1, prob = 0.8) # if outcome is 1 reward was given,
        PE <- outcome[t] - evA[t]
      } else { # if choice is B
        outcome[t] <- rbinom(1, size = 1, prob = 0.2) # if outcome is 1 reward was given
        PE <- outcome[t] - evB[t]
      }
    }
    
    if (t > 150 && t <= 175){
      if (choice[t] == 1){ # if choice is A
        outcome[t] <- rbinom(1, size = 1, prob = 0.2) # if outcome is 1 reward was given,
        PE <- outcome[t] - evA[t]
      } else { # if choice is B
        outcome[t] <- rbinom(1, size = 1, prob = 0.8) # if outcome is 1 reward was given
        PE <- outcome[t] - evB[t]
      }
    }
    
    if (t > 175 && t <= 200){
      if (choice[t] == 1){ # if choice is A
        outcome[t] <- rbinom(1, size = 1, prob = 0.8) # if outcome is 1 reward was given,
        PE <- outcome[t] - evA[t]
      } else { # if choice is B
        outcome[t] <- rbinom(1, size = 1, prob = 0.2) # if outcome is 1 reward was given
        PE <- outcome[t] - evB[t]
      }
    }
    
    if (t > 200 && t <= 225){
      if (choice[t] == 1){ # if choice is A
        outcome[t] <- rbinom(1, size = 1, prob = 0.2) # if outcome is 1 reward was given,
        PE <- outcome[t] - evA[t]
      } else { # if choice is B
        outcome[t] <- rbinom(1, size = 1, prob = 0.8) # if outcome is 1 reward was given
        PE <- outcome[t] - evB[t]
      }
    }
    
    if (t > 225){
      if (choice[t] == 1){ # if choice is A
        outcome[t] <- rbinom(1, size = 1, prob = 0.8) # if outcome is 1 reward was given,
        PE <- outcome[t] - evA[t]
      } else { # if choice is B
        outcome[t] <- rbinom(1, size = 1, prob = 0.2) # if outcome is 1 reward was given
        PE <- outcome[t] - evB[t]
      }
    }
    
    if (choice[t] == 1){ # if choice is A
      evA[t+1] = evA[t] + Alpha_sim[i] * PE # Update expected value
      evB[t+1] = evB[t]
    } else { # if choice is B
      evB[t+1] = evB[t] + Alpha_sim[i] * PE
      evA[t+1] = evA[t]
    }
    PEs[t] <- PE
  }
  all_evA[,i] <- evA
  all_evB[,i] <- evB
  all_choice[,i] <- choice
  all_outcome[,i] <- outcome
  all_PEs[,i] <- PEs
  
}


################################ prepare data for stan models
IDs = N
choice <- t(all_choice)
outcome <- t(all_outcome)

choice[choice == 0] <- 2
outcome[outcome == 0] <- -1

stan_data = list(N = as.numeric(IDs),
                 T = as.numeric(250),
                 choice = choice,
                 outcome = outcome)

################################ FM1_GM1
fit_simulated_M1_GM1 <- stan(file = "bandit_2arm_M1.stan",
                             data = stan_data,
                             cores = 4,
                             chain = 4,
                             thin = 1,
                             warmup = 5000,
                             iter   = 6000,
                             save_warmup = F)
#
save(fit_simulated_M1_GM1, file = "fit_M1_GM1.RData")
rm(fit_simulated_M1_GM1)

################################ FM3_GM1
fit_simulated_M3_GM1 <- stan(file = "bandit_2arm_M3.stan",
                             data = stan_data,
                             cores = 4,
                             chain = 4,
                             thin = 1,
                             warmup = 5000,
                             iter   = 6000,
                             save_warmup = F)
#
save(fit_simulated_M3_GM1, file = "fit_M3_GM1.RData")
rm(fit_simulated_M3_GM1)

################################ FM4_GM1
fit_simulated_M4_GM1 <- stan(file = "bandit_2arm_M4.stan",
                             data = stan_data,
                             cores = 4,
                             chain = 4,
                             thin = 1,
                             warmup = 5000,
                             iter   = 6000,
                             save_warmup = F)
#
save(fit_simulated_M4_GM1, file = "fit_M4_GM1.RData")
rm(fit_simulated_M4_GM1)

################################ FM5_GM1
fit_simulated_M5_GM1 <- stan(file = "bandit_2arm_M5.stan",
                             data = stan_data,
                             cores = 4,
                             chain = 4,
                             thin = 1,
                             warmup = 5000,
                             iter   = 6000,
                             save_warmup = F)
#
save(fit_simulated_M5_GM1, file = "fit_M5_GM1.RData")
rm(fit_simulated_M5_GM1)


################################################################################## Generate Data M3
rm(list=ls()) # clear all
set.seed(999)

N <- 100
fit_params <- read.csv('fit_params_M3.csv') # Read fit real data

fit_params_cov <- cov(fit_params)
Alphapos <- fit_params$Alphapos
Alphaneg <- fit_params$Alphaneg
Beta <- fit_params$Beta

params <- mvrandn(c(min(Alphapos),min(Alphaneg),min(Beta)),c(max(Alphapos),max(Alphaneg),max(Beta)), fit_params_cov, N, c(mean(Alphapos), mean(Alphaneg), mean(Beta)))

Alphapos_sim <- params[1,]
Alphaneg_sim <- params[2,]
Beta_sim <- params[3,]

outcome <- vector()
choice <- vector()
PEs <- vector()
evA <- 0
evB <- 0


all_evA <- array(0, dim = c(251,N))
all_evB <- array(0, dim = c(251,N))
all_choice <- array(0, dim = c(250,N))
all_outcome <- array(0, dim = c(250,N))
all_PEs <- array(0, dim = c(250,N))


for (i in 1:N){
  for (t in 1:250){
    
    probA <- exp(Beta_sim[i]*evA[t])/(exp(Beta_sim[i]*evA[t])+exp(Beta_sim[i]*evB[t]))
    
    # if choice = 1 choice was A, else choice was B
    choice[t] <- rbinom(1, size = 1, prob = probA)
    
    
    if (t <= 150){
      if (choice[t] == 1){ # if choice is A
        outcome[t] <- rbinom(1, size = 1, prob = 0.8) # if outcome is 1 reward was given,
        PE <- outcome[t] - evA[t]
      } else { # if choice is B
        outcome[t] <- rbinom(1, size = 1, prob = 0.2) # if outcome is 1 reward was given
        PE <- outcome[t] - evB[t]
      }
    }
    
    if (t > 150 && t <= 175){
      if (choice[t] == 1){ # if choice is A
        outcome[t] <- rbinom(1, size = 1, prob = 0.2) # if outcome is 1 reward was given,
        PE <- outcome[t] - evA[t]
      } else { # if choice is B
        outcome[t] <- rbinom(1, size = 1, prob = 0.8) # if outcome is 1 reward was given
        PE <- outcome[t] - evB[t]
      }
    }
    
    if (t > 175 && t <= 200){
      if (choice[t] == 1){ # if choice is A
        outcome[t] <- rbinom(1, size = 1, prob = 0.8) # if outcome is 1 reward was given,
        PE <- outcome[t] - evA[t]
      } else { # if choice is B
        outcome[t] <- rbinom(1, size = 1, prob = 0.2) # if outcome is 1 reward was given
        PE <- outcome[t] - evB[t]
      }
    }
    
    if (t > 200 && t <= 225){
      if (choice[t] == 1){ # if choice is A
        outcome[t] <- rbinom(1, size = 1, prob = 0.2) # if outcome is 1 reward was given,
        PE <- outcome[t] - evA[t]
      } else { # if choice is B
        outcome[t] <- rbinom(1, size = 1, prob = 0.8) # if outcome is 1 reward was given
        PE <- outcome[t] - evB[t]
      }
    }
    
    if (t > 225){
      if (choice[t] == 1){ # if choice is A
        outcome[t] <- rbinom(1, size = 1, prob = 0.8) # if outcome is 1 reward was given,
        PE <- outcome[t] - evA[t]
      } else { # if choice is B
        outcome[t] <- rbinom(1, size = 1, prob = 0.2) # if outcome is 1 reward was given
        PE <- outcome[t] - evB[t]
      }
    }
    
    # Positive outcome
    if (choice[t] == 1 & outcome[t] == 1){ # if choice is A and positive outcome
      evA[t+1] = evA[t] + Alphapos_sim[i] * PE
      evB[t+1] = evB[t]
    } else if (choice[t] == 0 & outcome[t] == 1) { # if choice is B and positive outcome
      evB[t+1] = evB[t] + Alphapos_sim[i] * PE
      evA[t+1] = evA[t]
    }
    # Negative outcome
    if (choice[t] == 1 & outcome[t] == 0){ # if choice is A and negative outcome
      evA[t+1] = evA[t] + Alphaneg_sim[i] * PE
      evB[t+1] = evB[t]
    } else if (choice[t] == 0 & outcome[t] == 0) { # if choice is B and negative outcome
      evB[t+1] = evB[t] + Alphaneg_sim[i] * PE
      evA[t+1] = evA[t]
    }
    PEs[t] <- PE
  }
  all_evA[,i] <- evA
  all_evB[,i] <- evB
  all_choice[,i] <- choice
  all_outcome[,i] <- outcome
  all_PEs[,i] <- PEs
  
}

all_outcome_mean <- colMeans(all_outcome)
df <- as.data.frame(all_outcome_mean)
df$alpha <- Alphaneg_sim
ggplot(data=df, aes(alpha,all_outcome_mean), fill=ID) +
  geom_point() +
  geom_smooth(method="lm",formula=y~I(x)+I(x^2), color="black") +
  theme_classic()
model <- lm(df$all_outcome_mean ~ df$alpha + I(df$alpha^2))
summary(model)

################################ prepare data
IDs = N
choice <- t(all_choice)
outcome <- t(all_outcome)

choice[choice == 0] <- 2
outcome[outcome == 0] <- -1

stan_data = list(N = as.numeric(IDs),
                 T = as.numeric(250),
                 choice = choice,
                 outcome = outcome)

################################ FM1_GM3
fit_simulated_M1_GM3 <- stan(file = "bandit_2arm_M1.stan",
                             data = stan_data,
                             cores = 4,
                             chain = 4,
                             thin = 1,
                             warmup = 5000,
                             iter   = 6000,
                             save_warmup = F)
#
save(fit_simulated_M1_GM3, file = "fit_M1_GM3.RData")
rm(fit_simulated_M1_GM3)

################################ FM3_GM3
fit_simulated_M3_GM3 <- stan(file = "bandit_2arm_M3.stan",
                             data = stan_data,
                             cores = 4,
                             chain = 4,
                             thin = 1,
                             warmup = 5000,
                             iter   = 6000,
                             save_warmup = F)
#
save(fit_simulated_M3_GM3, file = "fit_M3_GM3.RData")
rm(fit_simulated_M3_GM3)

################################ FM4_GM3
fit_simulated_M4_GM3 <- stan(file = "bandit_2arm_M4.stan",
                             data = stan_data,
                             cores = 4,
                             chain = 4,
                             thin = 1,
                             warmup = 5000,
                             iter   = 6000,
                             save_warmup = F)
#
save(fit_simulated_M4_GM3, file = "fit_M4_GM3.RData")
rm(fit_simulated_M4_GM3)

################################ FM5_GM3
fit_simulated_M5_GM3 <- stan(file = "bandit_2arm_M5.stan",
                             data = stan_data,
                             cores = 4,
                             chain = 4,
                             thin = 1,
                             warmup = 5000,
                             iter   = 6000,
                             save_warmup = F)
#
save(fit_simulated_M5_GM3, file = "fit_M5_GM3.RData")
rm(fit_simulated_M5_GM3)


################################################################################## Generate Data M4
rm(list=ls()) # clear all
set.seed(999)

N <- 100
fit_params <- read.csv('fit_params_M4.csv') # Read fit real data
fit_params_cov <- cov(fit_params)
beta <- fit_params$beta
Alpha0 <- fit_params$Alpha0
kappa <- fit_params$kappa

params <- mvrandn(c(min(beta),min(Alpha0),min(kappa)),c(max(beta),max(Alpha0),max(kappa)), fit_params_cov, N, c(mean(beta),mean(Alpha0),mean(kappa)))
beta <- params[1,]
Alpha0 <- params[2,]
kappa <- params[3,]

initV <- 0
C0 <- 0.5
gamma <- 0.5
C <- C0
outcome <- vector()
choice <- vector()
PEs <- vector()
alpha <- vector()
evA <- 0
evB <- 0


all_evA <- array(0, dim = c(251,N))
all_evB <- array(0, dim = c(251,N))
all_choice <- array(0, dim = c(250,N))
all_outcome <- array(0, dim = c(250,N))
all_PEs <- array(0, dim = c(250,N))
all_alpha <- array(0, dim = c(250,N))

for (i in 1:N){
  for (t in 1:250){
    
    probA <- exp(beta[i]*evA[t])/(exp(beta[i]*evA[t])+exp(beta[i]*evB[t]))
    
    # if choice = 1 choice was A, else choice was B
    choice[t] <- rbinom(1, size = 1, prob = probA)
    
    
    if (t <= 150){
      if (choice[t] == 1){ # if choice is A
        outcome[t] <- rbinom(1, size = 1, prob = 0.8) # if outcome is 1 reward was given,
        PE <- outcome[t] - evA[t]
      } else { # if choice is B
        outcome[t] <- rbinom(1, size = 1, prob = 0.2) # if outcome is 1 reward was given
        PE <- outcome[t] - evB[t]
      }
    }
    
    if (t > 150 && t <= 175){
      if (choice[t] == 1){ # if choice is A
        outcome[t] <- rbinom(1, size = 1, prob = 0.2) # if outcome is 1 reward was given,
        PE <- outcome[t] - evA[t]
      } else { # if choice is B
        outcome[t] <- rbinom(1, size = 1, prob = 0.8) # if outcome is 1 reward was given
        PE <- outcome[t] - evB[t]
      }
    }
    
    if (t > 175 && t <= 200){
      if (choice[t] == 1){ # if choice is A
        outcome[t] <- rbinom(1, size = 1, prob = 0.8) # if outcome is 1 reward was given,
        PE <- outcome[t] - evA[t]
      } else { # if choice is B
        outcome[t] <- rbinom(1, size = 1, prob = 0.2) # if outcome is 1 reward was given
        PE <- outcome[t] - evB[t]
      }
    }
    
    if (t > 200 && t <= 225){
      if (choice[t] == 1){ # if choice is A
        outcome[t] <- rbinom(1, size = 1, prob = 0.2) # if outcome is 1 reward was given,
        PE <- outcome[t] - evA[t]
      } else { # if choice is B
        outcome[t] <- rbinom(1, size = 1, prob = 0.8) # if outcome is 1 reward was given
        PE <- outcome[t] - evB[t]
      }
    }
    
    if (t > 225){
      if (choice[t] == 1){ # if choice is A
        outcome[t] <- rbinom(1, size = 1, prob = 0.8) # if outcome is 1 reward was given,
        PE <- outcome[t] - evA[t]
      } else { # if choice is B
        outcome[t] <- rbinom(1, size = 1, prob = 0.2) # if outcome is 1 reward was given
        PE <- outcome[t] - evB[t]
      }
    }
    
    
    # Update confidence
    C = C + gamma * ((2 - abs(PE))/2 - C)
    
    # Adapt learning rate
    if (PE >= 0){
      alpha[t] = (Alpha0[i] + kappa[i] * C)/(1 + kappa[i] * C)
    } else {
      alpha[t] = Alpha0[i]/(1 + kappa[i] * C)
    }
    
    # Update value
    if (choice[t] == 1){ # if choice is A
      evA[t+1] = evA[t] + alpha[i] * PE
      evB[t+1] = evB[t]
    } else { # if choice is B
      evB[t+1] = evB[t] + alpha[i] * PE
      evA[t+1] = evA[t]
    }
    PEs[t] <- PE
  }
  
  all_evA[,i] <- evA
  all_evB[,i] <- evB
  all_choice[,i] <- choice
  all_outcome[,i] <- outcome
  all_PEs[,i] <- PEs
  all_alpha[,i] <- alpha 
  
}

################################ prepare data
IDs = N
choice <- t(all_choice)
outcome <- t(all_outcome)

choice[choice == 0] <- 2
outcome[outcome == 0] <- -1

stan_data = list(N = as.numeric(IDs),
                 T = as.numeric(250),
                 choice = choice,
                 outcome = outcome)


################################ FM1_GM4
fit_simulated_M1_GM4 <- stan(file = "bandit_2arm_M1.stan",
                             data = stan_data,
                             cores = 4,
                             chain = 4,
                             thin = 1,
                             warmup = 5000,
                             iter   = 6000,
                             save_warmup = F)
#
save(fit_simulated_M1_GM4, file = "fit_M1_GM4.RData")
rm(fit_simulated_M1_GM4)

################################ FM3_GM4
fit_simulated_M3_GM4 <- stan(file = "bandit_2arm_M3.stan",
                             data = stan_data,
                             cores = 4,
                             chain = 4,
                             thin = 1,
                             warmup = 5000,
                             iter   = 6000,
                             save_warmup = F)
#
save(fit_simulated_M3_GM4, file = "fit_M3_GM4.RData")
rm(fit_simulated_M3_GM4)

################################ FM4_GM4
fit_simulated_M4_GM4 <- stan(file = "bandit_2arm_M4.stan",
                             data = stan_data,
                             cores = 4,
                             chain = 4,
                             thin = 1,
                             warmup = 5000,
                             iter   = 6000,
                             save_warmup = F)
#
save(fit_simulated_M4_GM4, file = "fit_M4_GM4.RData")
rm(fit_simulated_M4_GM4)

################################ FM5_GM4
fit_simulated_M5_GM4 <- stan(file = "bandit_2arm_M5.stan",
                             data = stan_data,
                             cores = 4,
                             chain = 4,
                             thin = 1,
                             warmup = 5000,
                             iter   = 6000,
                             save_warmup = F)
#
save(fit_simulated_M5_GM4, file = "fit_M5_GM4.RData")
rm(fit_simulated_M5_GM4)

################################################################################## Generate Data M5
rm(list=ls()) # clear all
set.seed(999)

N <- 100
fit_params <- read.csv('fit_params_M5.csv') # Read fit real data
fit_params_cov <- cov(fit_params)
beta <- fit_params$beta
Alpha0 <- fit_params$Alpha0
kappa <- fit_params$kappa

params <- mvrandn(c(min(beta),min(Alpha0),min(kappa)),c(max(beta),max(Alpha0),max(kappa)), fit_params_cov, N, c(mean(beta),mean(Alpha0),mean(kappa)))
beta <- params[1,]
Alpha0 <- params[2,]
kappa <- params[3,]

initV <- 0
C0 <- 0.5
gamma <- 0.1
C <- C0
outcome <- vector()
choice <- vector()
PEs <- vector()
alpha <- vector()
evA <- 0
evB <- 0


all_evA <- array(0, dim = c(251,N))
all_evB <- array(0, dim = c(251,N))
all_choice <- array(0, dim = c(250,N))
all_outcome <- array(0, dim = c(250,N))
all_PEs <- array(0, dim = c(250,N))
all_alpha <- array(0, dim = c(250,N))


for (i in 1:N){
  for (t in 1:250){
    
    probA <- exp(beta[i]*evA[t])/(exp(beta[i]*evA[t])+exp(beta[i]*evB[t]))
    
    # if choice = 1 choice was A, else choice was B
    choice[t] <- rbinom(1, size = 1, prob = probA)
    
    
    if (t <= 150){
      if (choice[t] == 1){ # if choice is A
        outcome[t] <- rbinom(1, size = 1, prob = 0.8) # if outcome is 1 reward was given,
        PE <- outcome[t] - evA[t]
      } else { # if choice is B
        outcome[t] <- rbinom(1, size = 1, prob = 0.2) # if outcome is 1 reward was given
        PE <- outcome[t] - evB[t]
      }
    }
    
    if (t > 150 && t <= 175){
      if (choice[t] == 1){ # if choice is A
        outcome[t] <- rbinom(1, size = 1, prob = 0.2) # if outcome is 1 reward was given,
        PE <- outcome[t] - evA[t]
      } else { # if choice is B
        outcome[t] <- rbinom(1, size = 1, prob = 0.8) # if outcome is 1 reward was given
        PE <- outcome[t] - evB[t]
      }
    }
    
    if (t > 175 && t <= 200){
      if (choice[t] == 1){ # if choice is A
        outcome[t] <- rbinom(1, size = 1, prob = 0.8) # if outcome is 1 reward was given,
        PE <- outcome[t] - evA[t]
      } else { # if choice is B
        outcome[t] <- rbinom(1, size = 1, prob = 0.2) # if outcome is 1 reward was given
        PE <- outcome[t] - evB[t]
      }
    }
    
    if (t > 200 && t <= 225){
      if (choice[t] == 1){ # if choice is A
        outcome[t] <- rbinom(1, size = 1, prob = 0.2) # if outcome is 1 reward was given,
        PE <- outcome[t] - evA[t]
      } else { # if choice is B
        outcome[t] <- rbinom(1, size = 1, prob = 0.8) # if outcome is 1 reward was given
        PE <- outcome[t] - evB[t]
      }
    }
    
    if (t > 225){
      if (choice[t] == 1){ # if choice is A
        outcome[t] <- rbinom(1, size = 1, prob = 0.8) # if outcome is 1 reward was given,
        PE <- outcome[t] - evA[t]
      } else { # if choice is B
        outcome[t] <- rbinom(1, size = 1, prob = 0.2) # if outcome is 1 reward was given
        PE <- outcome[t] - evB[t]
      }
    }
    
    
    # Update confidence
    C = C + gamma * ((2 - abs(PE))/2 - C)
    
    # Adapt learning rate
    if (PE >= 0){
      alpha[t] = (Alpha0[i] + kappa[i] * C)/(1 + kappa[i] * C)
    } else {
      alpha[t] = Alpha0[i]/(1 + kappa[i] * C)
    }
    
    # Update value
    if (choice[t] == 1){ # if choice is A
      evA[t+1] = evA[t] + alpha[i] * PE
      evB[t+1] = evB[t]
    } else { # if choice is B
      evB[t+1] = evB[t] + alpha[i] * PE
      evA[t+1] = evA[t]
    }
    PEs[t] <- PE
  }
  
  all_evA[,i] <- evA
  all_evB[,i] <- evB
  all_choice[,i] <- choice
  all_outcome[,i] <- outcome
  all_PEs[,i] <- PEs
  all_alpha[,i] <- alpha 
  
}

################################ prepare data
IDs = N
choice <- t(all_choice)
outcome <- t(all_outcome)

choice[choice == 0] <- 2
outcome[outcome == 0] <- -1

stan_data = list(N = as.numeric(IDs),
                 T = as.numeric(250),
                 choice = choice,
                 outcome = outcome)


################################ FM1_GM5
fit_simulated_M1_GM5 <- stan(file = "bandit_2arm_M1.stan",
                             data = stan_data,
                             cores = 4,
                             chain = 4,
                             thin = 1,
                             warmup = 5000,
                             iter   = 6000,
                             save_warmup = F)
#
save(fit_simulated_M1_GM5, file = "fit_M1_GM5.RData")
rm(fit_simulated_M1_GM5)

################################ FM3_GM5
fit_simulated_M3_GM5 <- stan(file = "bandit_2arm_M3.stan",
                             data = stan_data,
                             cores = 4,
                             chain = 4,
                             thin = 1,
                             warmup = 5000,
                             iter   = 6000,
                             save_warmup = F)
#
save(fit_simulated_M3_GM5, file = "fit_M3_GM5.RData")
rm(fit_simulated_M3_GM5)

################################ FM4_GM5
fit_simulated_M4_GM5 <- stan(file = "bandit_2arm_M4.stan",
                             data = stan_data,
                             cores = 4,
                             chain = 4,
                             thin = 1,
                             warmup = 5000,
                             iter   = 6000,
                             save_warmup = F)
#
save(fit_simulated_M4_GM5, file = "fit_M4_GM5.RData")
rm(fit_simulated_M4_GM5)

################################ FM5_GM5
fit_simulated_M5_GM5 <- stan(file = "bandit_2arm_M5.stan",
                             data = stan_data,
                             cores = 4,
                             chain = 4,
                             thin = 1,
                             warmup = 5000,
                             iter   = 6000,
                             save_warmup = F)
#
save(fit_simulated_M5_GM5, file = "fit_M5_GM5.RData")
rm(fit_simulated_M5_GM5)



