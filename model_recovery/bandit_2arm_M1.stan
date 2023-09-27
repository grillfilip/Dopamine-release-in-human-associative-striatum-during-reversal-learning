// Simple delta rule model
data {
  int<lower=1> N;
  int<lower=1> T;
  // int<lower=1, upper=T> Tsubj[N];
  int<lower=-1, upper=2> choice[N, T];
  real outcome[N, T];  // no lower and upper bounds
}
transformed data {
  vector[2] initV;  // initial values for EV
  initV = rep_vector(0, 2);
}
parameters {
  // Declare all parameters as vectors for vectorizing
  // Hyper(group)-parameters
  vector[2] mu_pr;
  vector<lower=0>[2] sigma;
  
  // Subject-level raw parameters (for Matt trick)
  vector[N] alpha_pr;    // learning rate
  vector[N] beta_pr;  // inverse temperature
}
transformed parameters {
  // subject-level parameters
  vector<lower=0, upper=1>[N] alpha;
  vector<lower=0>[N] beta;
  
  for (i in 1:N) {
    alpha[i]   = Phi_approx(mu_pr[1]  + sigma[1]  * alpha_pr[i]);
    beta[i] = exp(mu_pr[2] + sigma[2] * beta_pr[i]);
  }
}
model {
  // Hyperparameters
  mu_pr  ~ normal(0, 1);
  sigma ~ normal(0, 0.2);
  
  // individual parameters
  alpha_pr   ~ normal(0, 1);
  beta_pr ~ normal(0, 1);
  
  // subject loop and trial loop
  for (i in 1:N) {
    vector[2] ev; // expected value
    real PE;      // prediction error
    
    ev = initV;
    
    for (t in 1:250) {
      if (outcome[i, t] != 0){
              // compute action probabilities
      choice[i, t] ~ categorical_logit(beta[i] * ev);
      
      // prediction error
      PE = outcome[i, t] - ev[choice[i, t]];
      
      // value updating (learning)
      ev[choice[i, t]] += alpha[i] * PE;
      }
    }
  }
}
generated quantities {
  // For group level parameters
  real<lower=0, upper=1> mu_alpha;
  real<lower=0> mu_beta;
  
  // For log likelihood calculation
  real log_lik[N];
  
  // For posterior predictive check
  real y_pred[N, T];
  real PEs[N, T];
  real EVAs[N, T];
  real EVBs[N, T];
  
  // Set all posterior predictions to 0 (avoids NULL values)
  for (i in 1:N) {
    for (t in 1:T) {
      y_pred[i, t] = -1;
      PEs[i, t] = -1;
      //Evs[i, t] = -1;
      EVAs[i, t] = 0;
      EVBs[i, t] = 0;
    }
  }
  
  mu_alpha  = Phi_approx(mu_pr[1]);
  mu_beta = exp(mu_pr[2]);
  
  { // local section, this saves time and space
    for (i in 1:N) {
      vector[2] ev; // expected value
      real PE;      // prediction error
      
      // Initialize values
      ev = initV;
      
      log_lik[i] = 0;
      
      for (t in 1:250) {
        if (outcome[i, t] != 0){
          
          // compute log likelihood of current trial
          log_lik[i] += categorical_logit_lpmf(choice[i, t] | beta[i] * ev);
          
          // generate posterior prediction for current trial
          y_pred[i, t] = categorical_rng(softmax(beta[i] * ev));
          
          // prediction error
          PE = outcome[i, t] - ev[choice[i, t]];
          PEs[i, t] = PE;
          
          // value updating (learning)
          ev[choice[i, t]] += alpha[i] * PE;
        
          EVAs[i,t] = ev[1];
          EVBs[i,t] = ev[2];

        }
      }
    }
  }
}
