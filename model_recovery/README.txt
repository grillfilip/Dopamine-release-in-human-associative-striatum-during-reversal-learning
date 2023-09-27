# Filip Grill 2023-09-19
Script and data to generate simulated parameters for parameter and model
recovery.

R version 4.2.2 Patched (2022-11-10 r83330)
STAN version 2.21.0

Files:
Generate_and_fit_simulated_data.R - main script
fit_params_MX.csv - parameters fit under real data
bandit_2arm_MX.stan - STAN models to fit simulated data
Model_recovery_BIC.R - Perform model recovery


Generate_and_fit_simulated_data.R reads in the fit_params_MX.csv and
generates new N paramaters that are bound by the min, max, mean, and
covariance of the real parameters. The simulated parameters are then
used to generate simulated data for each model. This data is then fit
across each model. The script outputs fit_MX_GMX.RData for each
combination of generated and fit models e.g. fit_M1_GM1.RData fits
model M1 on data generated from M1, fit_M2_GM1.RData fits model M2 on
data generated from M1.

Model_recovery_BIC.R takes the fit_MX_GMX.RData as input and counts how
many times a model is the best fit given the generating (true) model and
outputs a confusion matrix (BIC_confusion_matrix.csv) where rows are the
generating model and columns are the fit model.
