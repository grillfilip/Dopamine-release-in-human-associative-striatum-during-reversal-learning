// Vinckier model with confidence updating learning rate dynamically C0=0.5 gamma=0.1
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
  vector[3] mu_pr;
  vector<lower=0>[3] sigma;
  
  // Subject-level raw parameters (for Matt trick)
  vector[N] alpha0_pr;    // learning rate
  vector[N] beta_pr;  // inverse temperature
  //vector[N] C0_pr;
  //vector[N] gamma_pr;
  vector[N] kappa_pr;
}
transformed parameters {
  // subject-level parameters
  vector<lower=0, upper=1>[N] alpha0;
  vector<lower=0>[N] beta;
  //vector<lower=0, upper=1>[N] C0;
  //vector<lower=0, upper=1>[N] gamma;
  vector<lower=0, upper=1>[N] kappa;
  
  for (i in 1:N) {
    alpha0[i]   = Phi_approx(mu_pr[1]  + sigma[1]  * alpha0_pr[i]);
    beta[i] = exp(mu_pr[2] + sigma[2] * beta_pr[i]);
    //C0[i]       = Phi_approx(mu_pr[3] + sigma[3] * C0_pr[i]);
    //gamma[i]    = Phi_approx(mu_pr[4] + sigma[4] * gamma_pr[i]);
    kappa[i]    = Phi_approx(mu_pr[3] + sigma[3] * kappa_pr[i]);
  }
}
model {
  // Hyperparameters
  mu_pr  ~ normal(0, 1);
  sigma ~ normal(0, 0.2);
  
  // individual parameters
  alpha0_pr   ~ normal(0, 1);
  beta_pr ~ normal(0, 1);
  //C0_pr        ~ normal(0, 1.0);
  //gamma_pr     ~ normal(0, 1.0);
  kappa_pr     ~ normal(0, 1.0);
  
  // subject loop and trial loop
  for (i in 1:N) {
    vector[2] ev; // expected value
    real PE;      // prediction error
    real C;
    real alpha;  // effective alpha
    real C0;
    real gamma;
    
    ev = initV;
    C0 = 0.5;
    gamma = 0.1;
    C = C0;
    
    for (t in 1:250) {
      if (outcome[i, t] != 0){
        // compute action probabilities
        choice[i, t] ~ categorical_logit(beta[i] * ev);
        
        // prediction error
        PE = outcome[i, t] - ev[choice[i, t]];
        
        // updates confidence based on PE
        C = C + gamma * ((2 - fabs(PE))/2 - C);
        
        // adapts learning rate according to confidence and feedback type
        if (PE >= 0){
          alpha = (alpha0[i] + kappa[i]*C)/(1 + kappa[i]*C);
        }
        if (PE < 0){
          alpha = alpha0[i]/(1 + kappa[i] * C);
        }
        
        // value updating (learning)
        ev[choice[i, t]] += alpha * PE;
      }
    }
  }
}
generated quantities {
  // For group level parameters
  real<lower=0, upper=1> mu_alpha;
  real<lower=0> mu_beta;
 // real<lower=0, upper=1> mu_C0;
 // real<lower=0, upper=1> mu_gamma;
  real<lower=0, upper=1> mu_kappa;
  
  // For log likelihood calculation
  real log_lik[N];
  
  // For posterior predictive check
  real y_pred[N, T];
  
  // Set all posterior predictions to 0 (avoids NULL values)
  for (i in 1:N) {
    for (t in 1:T) {
      y_pred[i, t] = -1;
    }
  }
  
  mu_alpha    = Phi_approx(mu_pr[1]);
  mu_beta     = exp(mu_pr[2]);
  //mu_C0       = Phi_approx(mu_pr[3]);
  //mu_gamma    = Phi_approx(mu_pr[4]);
  mu_kappa    = Phi_approx(mu_pr[3]);
  
  //real alphas[N,T];
  //real PEs[N,T];
  //real Cs[N,T];
   // local section, this saves time and space
   
    for (i in 1:N) {
      vector[2] ev; // expected value
      real PE;      // prediction error
      real C;
      real alpha;
      real C0;
      real gamma;

      
      // Initialize values
      ev = initV;
      C0 = 0.5;
      gamma = 0.5;
      C = C0;
      
      log_lik[i] = 0;
      
      for (t in 1:250) {
        // compute log likelihood of current trial
        log_lik[i] += categorical_logit_lpmf(choice[i, t] | beta[i] * ev);
        
        // generate posterior prediction for current trial
        y_pred[i, t] = categorical_rng(softmax(beta[i] * ev));
        
        // prediction error
        PE = outcome[i, t] - ev[choice[i, t]];
        
        // update confidence
        C = C + gamma * ((2 - fabs(PE))/2 - C);
      
        // adapts learning rate according to confidence
        if (PE >= 0){
          alpha = (alpha0[i] + kappa[i]*C)/(1 + kappa[i]*C);
        }
        if (PE < 0){
          alpha = alpha0[i]/(1 + kappa[i] * C);
        }
        
        // value updating (learning)
        ev[choice[i, t]] += alpha * PE;
        
        // saving variables
        //alphas[i, t] = alpha;
        //PEs[i, t] = PE;
        //Cs[i, t] = C;
      }
    }
  
}
