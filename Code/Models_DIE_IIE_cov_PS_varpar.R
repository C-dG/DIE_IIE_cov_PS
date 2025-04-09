#===============================================================================

# Rstan models MS variance partitioning models PS DIE-IIE covariance 
# Producing & Scrounging sparrow project

# Author: Corné de Groot
# Readme file for variable details

#===============================================================================

# Save local package information:
#install.packages("renv")
#renv::init()

# Load local package information:
renv::restore()

# Load in neccessary packages
library(rstan)
library(shinystan)
library(dplyr)
library(tidyr)
library(parallelly)
library(beepr)

# Package versions for reproducibility
packageVersion("rstan") #2.32.6
packageVersion("shinystan") #2.6.0
packageVersion("dplyr") #1.1.3
packageVersion("tidyr") #1.3.0
packageVersion("parallelly") #1.36.0
packageVersion("beepr") #1.3

R.Version() #"R version 4.3.1 

# Load in data 
#===============================================================================

# Set wd if loading the df does not work
#setwd(rstudioapi::getActiveProject())

df <- read.csv("Data/Data_PS_DIE_IIE_cov.csv", header = TRUE)

# Make new directories to save the models and their output
ifelse(file.exists("Model_output"), "Dir already exists",
       dir.create("Model_output"))

ifelse(file.exists("Models"), "Dir already exists",
       dir.create("Models"))

#===============================================================================

# Legend of models

# Bivariate model with full fixed effect structure
# Bivariate model with only random effect structure
# Univariate models for cross-year repeatability

##### Bivariate model #####

# Create the dataset used in the model fit
stan_bivar_data <- list(No = nrow(df),
                   Ni = length(unique(df$RingNR)),
                  individual = as.integer(as.factor(df$RingNR)),
                  opponent1 = as.integer(as.factor(df$Opponent1_ID)),
                  opponent2 = as.integer(as.factor(df$Opponent2_ID)),
                  trialID = as.integer(as.factor(df$TrialID)),
                  n_trials = length(unique(df$TrialID)),
                  groupID = as.integer(as.factor(df$GroupID)),
                  n_groups = length(unique(df$GroupID)),
                  trial_order = as.vector(scale(df$TrialNR_Indiv, center = T, scale = T)),
                  trialday = ifelse(df$Trial_day_bin == 0, -0.5, 0.5),
                  group_trial_b4 = df$Grp_tr_b4,
                  year = ifelse(df$Year == 2022, -0.5, 0.5),
                  zs = as.matrix(cbind(as.vector(scale(df$Producing_events)),
                                       as.vector(scale(df$Scrounging_events)))))

# The model code
write(temp <- "data {
  // Number of clusters
  int<lower=1> No;  // Total number of observations
  int<lower=1> Ni; // Number of individuals
  int<lower=1> n_trials; // Number of trials
  int<lower=1> n_groups; // Number of groups
  // Cluster identifier
  int<lower=1> individual[No]; // Individual identity for each observation
  int<lower=1> opponent1[No]; // Opponent identity for each observation
  int<lower=1> opponent2[No]; // Opponent identity for each observation
  int<lower=1> trialID[No];  //  TrialID of each observation
  int<lower=1> groupID[No];  //  GroupID of each observation
  // Predictors
  vector[No] year; // Year 0/1 (2022/2023)for each observation
  vector[No] trial_order; // Trial order/sequence per individual for each observation
  vector[No] trialday; // Trialday 0/1 for each observation
  vector[No] group_trial_b4; // Group trial before 0/1 for each observation
  // Response variable
  array[No] vector[2] zs; // Phenotypic observations of producing & scrounging 
}
transformed data {
// transform interactions into singular vector for faster interaction estimation
  vector[No] trial_order_trialday = trial_order .* trialday;  
  vector[No] trial_order_group_trial_b4 = trial_order .* group_trial_b4; 
  vector[No] group_trial_b4_trialday = group_trial_b4 .* trialday;  
  vector[No] trial_order_trialday_group_trial_b4 = trial_order .* trialday .* group_trial_b4;  
}
parameters {
  // Fixed effects
  real mu; // overall intercept
  real B_year; // slope of year
  real B_trial_order;// slope of trial order
  real B_trialday; // Slope of trialday
  real B_group_trial_b4; // Slope of individuals that participated in group trials
  real mu_2; // Overall intercept
  real B_year_2; // Slope of year
  real B_trial_order_2;// Slope of trial order
  real B_trialday_2; // Slope of trialday
  real B_group_trial_b4_2; // Slope of individuals that participated in group trials
  // Interactions
  real int_trialday_group_trial_b4_trialorder;
  real int_trialday_trialorder; 
  real int_trialday_group_trial_b4; 
  real int_group_trial_b4_trialorder;
  real int_trialday_group_trial_b4_trialorder_2;
  real int_trialday_trialorder_2; 
  real int_trialday_group_trial_b4_2; 
  real int_group_trial_b4_trialorder_2;
  // Random effects individual level matrix I
  vector<lower=0>[4] sigma_I; // Individual standard deviations DIE1, IIE1, DIE2, IIE2
  matrix[Ni,4] Iz; // Standardised  individual effects
  cholesky_factor_corr[4] LI; // Correlation/covariance factor
  // Other random effects
  real <lower=0> sigma_trial;
  real <lower=0> sigma_group;
  vector[n_trials] zTrial_I;// standardised to speed up computation
  vector[n_groups] zGroup_I; 
  // for equation2
  real <lower=0> sigma_trial_2;
  real <lower=0> sigma_group_2;
  vector[n_trials] zTrial_I_2; 
  vector[n_groups] zGroup_I_2; 
  // Residual variance & correlation
  cholesky_factor_corr[2] LR; // Cholesky corr matrix for residuals
  vector<lower=0>[2] sigma_R; // Standard deviation (i.e. residual variance)
}
transformed parameters {
  // Make I-matrix of individual level effects (DIEs & IIEs) 
  matrix[Ni,4] I = Iz * diag_pre_multiply(sigma_I, LI)';// Individual effects
  vector[n_trials] trial_I  =  zTrial_I * sigma_trial; // Get the unscaled values for random effect
  vector[n_groups] group_I  =  zGroup_I * sigma_group; 
  vector[n_trials] trial_I_2  =  zTrial_I_2 * sigma_trial_2; 
  vector[n_groups] group_I_2  =  zGroup_I_2 * sigma_group_2;
  // Expected trait values for each observation
  vector[No] z_exp; // expected phenotypic values producing (eq1)
  vector[No] z2_exp; // expected phenotypic values scrounging (eq2)
  //Model equations:
  // Partition the phenotypes into DIEs and IIEs for producing and scrounging
  // Equation 1: Producing
  z_exp = mu + B_year*year + B_trialday*trialday + 
               B_trial_order*trial_order + B_group_trial_b4*group_trial_b4 + 
               int_trialday_trialorder*trial_order_trialday + 
               int_trialday_group_trial_b4*group_trial_b4_trialday + 
               int_group_trial_b4_trialorder*trial_order_group_trial_b4 + 
               int_trialday_group_trial_b4_trialorder*trial_order_trialday_group_trial_b4 + 
               I[individual,1] + I[opponent1,2] + I[opponent2,2] + 
               trial_I[trialID] + group_I[groupID];
  // Equation 2: Scrounging
  z2_exp = mu_2 + B_year_2*year + B_trialday_2*trialday + 
                  B_trial_order_2*trial_order + B_group_trial_b4_2*group_trial_b4 + 
                  int_trialday_trialorder_2*trial_order_trialday +
                  int_trialday_group_trial_b4_2*trial_order_group_trial_b4 +
                  int_group_trial_b4_trialorder_2*group_trial_b4_trialday +
                  int_trialday_group_trial_b4_trialorder_2*trial_order_trialday_group_trial_b4 + 
                  I[individual,3] + I[opponent1,4] + I[opponent2,4] + 
                  trial_I_2[trialID] + group_I_2[groupID];
}
model {
  // Priors
  // Fixed effects
  mu ~ normal(0, 1); // Weakly informative prior on the intercept
  B_year ~ normal(0,1); 
  B_trialday ~ normal(0,1); 
  B_trial_order ~ normal(0,1); 
  B_group_trial_b4 ~ normal(0,1); 
  int_trialday_trialorder ~ normal(0,1); 
  int_trialday_group_trial_b4 ~ normal(0,1); 
  int_group_trial_b4_trialorder ~ normal(0,1); 
  int_trialday_group_trial_b4_trialorder ~ normal(0,1);
  mu_2 ~ normal(0, 1); 
  B_year_2 ~ normal(0,1); 
  B_trialday_2 ~ normal(0,1); 
  B_trial_order_2 ~ normal(0,1); 
  B_group_trial_b4_2 ~ normal(0,1); 
  int_trialday_trialorder_2 ~ normal(0,1); 
  int_trialday_group_trial_b4_2 ~ normal(0,1); 
  int_group_trial_b4_trialorder_2 ~ normal(0,1); 
  int_trialday_group_trial_b4_trialorder_2 ~ normal(0,1); 
  // Random effects
  to_vector(Iz) ~ normal(0,1); // Implies i ~ normal(0,sigma_I)
  to_vector(zTrial_I) ~ normal(0, 1);
	to_vector(zGroup_I) ~ normal(0, 1);
  to_vector(zTrial_I_2) ~ normal(0, 1);
	to_vector(zGroup_I_2) ~ normal(0, 1);
  to_vector(sigma_I) ~ exponential(3); // Weakly informative prior for the variances
  to_vector(sigma_R) ~ exponential(3); 
  sigma_trial ~ exponential(3);
  sigma_group ~ exponential(3);
  sigma_trial_2 ~ exponential(3);
  sigma_group_2 ~ exponential(3);
  LI ~ lkj_corr_cholesky(1); //Prior individual level correlation
  LR ~ lkj_corr_cholesky(1); //Prior residual correlation
  // Residual correlation Choleksky factor
  matrix[2,2] L_sigma = diag_pre_multiply(sigma_R, LR);
  // Transform expected values to array
  array[No] vector[2] mus;
  for (o in 1:No) 
  mus[o] = [z_exp[o], z2_exp[o]]';
  // Cholesky multinormal likelihood function to estimate residual correlations
  zs ~ multi_normal_cholesky(mus, L_sigma);
}
generated quantities {
  // Variances individual effects
  real var_DIE = sigma_I[1]^2;
  real var_IIE = sigma_I[2]^2;
  real var_DIE_2 = sigma_I[3]^2;
  real var_IIE_2 = sigma_I[4]^2;
  // Variances other random effects
  real var_trialID = sigma_trial^2;
  real var_groupID = sigma_group^2;
  real var_trialID_2 = sigma_trial_2^2;
  real var_groupID_2 = sigma_group_2^2;
  // Residual variation
  real var_res = sigma_R[1]^2; 
  real var_res_2 = sigma_R[2]^2; 
  // Total variation 
  real<lower=0> var_P = var_DIE + (var_IIE*2) + var_trialID + var_groupID + var_res;
  real<lower=0> var_P_2 = var_DIE_2 + (var_IIE_2*2) + var_trialID_2 + var_groupID_2 +var_res_2;
  // Repeatabilities
  real<lower=0> rep_DIE = var_DIE/var_P; // rep DIE
  real<lower=0> rep_IIE = var_IIE/var_P; // rep IIE
  real<lower=0> rep_trialID = var_trialID/var_P; // rep trialID
  real<lower=0> rep_groupID = var_groupID/var_P; // rep groupID
  real<lower=0> rep_res = var_res/var_P; // rep residual
  real<lower=0> rep_DIE_2 = var_DIE_2/var_P_2; // rep DIE
  real<lower=0> rep_IIE_2 = var_IIE_2/var_P_2; // rep IIE
  real<lower=0> rep_trialID_2 = var_trialID_2/var_P_2; // rep trialID
  real<lower=0> rep_groupID_2 = var_groupID_2/var_P_2; // rep groupID
  real<lower=0> rep_res_2 = var_res_2/var_P_2; // rep residual
  // derive P matrix
  matrix[4,4] Omega_P = LI * LI'; // Correlation matrix
  matrix[4,4] D_I = diag_matrix(sigma_I); // Diagonal SD matrix
  matrix[4,4] Cov_I = D_I*Omega_P*D_I; // Covariance matrix  
  // Covariance
  real cov_P1 = Cov_I[1,2];
  real cov_P2 = Cov_I[1,3];
  real cov_P3 = Cov_I[1,4];
  real cov_P4 = Cov_I[2,3];
  real cov_P5 = Cov_I[2,4];
  real cov_P6 = Cov_I[3,4];
  // Correlations
  real cor_P1 = Omega_P[1,2];
  real cor_P2 = Omega_P[1,3];
  real cor_P3 = Omega_P[1,4];
  real cor_P4 = Omega_P[2,3];
  real cor_P5 = Omega_P[2,4];
  real cor_P6 = Omega_P[3,4];
  // Residual correlation
  matrix[2,2] res_cor_mat = LR * LR';
  real<lower=-1,upper=1> res_cor = res_cor_mat[1,2];
  real res_cov = res_cor_mat[1,2]*sqrt(var_res*var_res_2);
  // Total phenotypic effect
  real total_P_eff_1 = var_DIE + (4*cov_P1) + (4*var_IIE);
  real total_P_eff_2 = var_DIE_2 + (4*cov_P6) + (4*var_IIE_2);
  real rep_P_eff_1 = total_P_eff_1/var_P;
  real rep_P_eff_2 = total_P_eff_2/var_P_2;
}", file = "Models/Cholesky_bivar_PS.stan")

# Fit the model
md_PS_bivar_full <- stan("Models/Cholesky_bivar_PS.stan", data = stan_bivar_data, 
                         chains = 5, iter = 5000,  
                         warmup = 2000, thin = 1,
                         save_warmup = F,
                         cores = parallelly::availableCores()-1,
                         seed = 534094228
                         )
beep(0)# Plays a sound when the model has been fit

# Derive seed if necessary to reproduce the same model output
rstan::get_seed(md_PS_bivar_full)#534094228

# Set the model parameters that need to be extracted from the model fit:
params_bivar <- c("mu", "B_year", "B_trial_order", "B_trialday", "B_group_trial_b4",
                  "int_trialday_group_trial_b4", "int_trialday_trialorder","int_group_trial_b4_trialorder", "int_trialday_group_trial_b4_trialorder",
                  "mu_2", "B_year_2", "B_trial_order_2", "B_trialday_2", "B_group_trial_b4_2",
                  "int_trialday_group_trial_b4_2", "int_trialday_trialorder_2","int_group_trial_b4_trialorder_2", "int_trialday_group_trial_b4_trialorder_2",
                  "var_DIE", "var_IIE", "var_trialID", "var_groupID", "var_res","var_P",
                  "var_DIE_2", "var_IIE_2", "var_trialID_2", "var_groupID_2","var_res_2", "var_P_2",
                  "rep_DIE", "rep_IIE", "rep_trialID", "rep_groupID", "rep_res",
                  "rep_DIE_2", "rep_IIE_2", "rep_trialID_2", "rep_groupID_2","rep_res_2", 
                  "cov_P1", "cov_P2", "cov_P3", "cov_P4" , "cov_P5", "cov_P6",
                  "cor_P1", "cor_P2", "cor_P3", "cor_P4" , "cor_P5", "cor_P6", 
                  "res_cov" ,"res_cor",
                  "total_P_eff_1", "total_P_eff_2", "rep_P_eff_1", "rep_P_eff_2","lp__")

# Derive the model output from the posterior (Mean + CI + median + ESS + Rhat)
round(summary(md_PS_bivar_full, pars = params_bivar)$summary[,c(1,4,6,8,9,10)],3)
Sum_md_PS_bivar_full <- round(summary(md_PS_bivar_full, pars = params_bivar)$summary[,c(1,4,6,8,9,10)],3)

# Look at posterior distribtions:
#launch_shinystan(md_PS_bivar_full)

# Extracting BLUPS
BLUPS_I_bivar_PS <- round(summary(md_PS_bivar_full, pars = "I")$summary[,c(1,4,6, 8, 9,10)],3)

# Save the model, summary & BLUPS
saveRDS(Sum_md_PS_bivar_full, file = paste("Model_output/Sum_md_PS_bivar_full.rds", sep = ""))
saveRDS(md_PS_bivar_full, file = paste("Model_output/md_PS_bivar_full.rds", sep = ""))
saveRDS(BLUPS_I_bivar_PS, file = paste("Model_output/BLUPS_I_bivar_PS.rds", sep = ""))

#===============================================================================

##### Unadjusted repeatability model #####

# Create the dataset used in the model fit
stan_bivar_data_rep_unadjusted <- list(No = nrow(df),
                                       Ni = length(unique(df$RingNR)),
                                       individual = as.integer(as.factor(df$RingNR)),
                                       opponent1 = as.integer(as.factor(df$Opponent1_ID)),
                                       opponent2 = as.integer(as.factor(df$Opponent2_ID)),
                                       trialID = as.integer(as.factor(df$TrialID)),
                                       n_trials = length(unique(df$TrialID)),
                                       groupID = as.integer(as.factor(df$GroupID)),
                                       n_groups = length(unique(df$GroupID)),
                                       zs = as.matrix(cbind(as.vector(scale(df$Producing_events)),
                                                            as.vector(scale(df$Scrounging_events)))))
# The model code
write(temp <- "data {
  // Number of clusters
  int<lower=1> No;  // Total number of observations
  int<lower=1> Ni; // Number of individuals
  int<lower=1> n_trials; // Number of trials
  int<lower=1> n_groups; // Number of groups
  // Cluster identifier
  int individual[No]; // Individual identity for each observation
  int<lower=1> opponent1[No]; // Opponent identity for each observation
  int<lower=1> opponent2[No]; // Opponent identity for each observation
  int<lower=1> trialID[No];  //  TrialID of each observation
  int<lower=1> groupID[No];  //  GroupID of each observation
  // Response variable
  array[No] vector[2] zs; // Phenotypic observations of producing & scrounging 
}
parameters {
  // Fixed effects
  real mu; // overall intercept
  real mu_2; // Overall intercept
  // Random effects individual level matrix I
  vector<lower=0>[4] sigma_I; // Individual standard deviations DIE1, IIE1, DIE2, IIE2
  matrix[Ni,4] Iz; // Standardised  individual effects
  cholesky_factor_corr[4] LI; // Correlation/covariance factor
  // Other random effects
  real <lower=0> sigma_trial;
  real <lower=0> sigma_group;
  vector[n_trials] zTrial_I;// standardised to speed up computation
  vector[n_groups] zGroup_I; 
  // for equation2
  real <lower=0> sigma_trial_2;
  real <lower=0> sigma_group_2;
  vector[n_trials] zTrial_I_2; 
  vector[n_groups] zGroup_I_2; 
  // Residual variance & correlation
  cholesky_factor_corr[2] LR; // Cholesky corr matrix for residuals
  vector<lower=0>[2] sigma_R; // Standard deviation (i.e. residual variance)
}
transformed parameters {
  // Make I-matrix of individual level effects (DIEs & IIEs) 
  matrix[Ni,4] I = Iz * diag_pre_multiply(sigma_I, LI)';// Individual effects
  vector[n_trials] trial_I  =  zTrial_I * sigma_trial; // Get the unscaled values for random effect
  vector[n_groups] group_I  =  zGroup_I * sigma_group; 
  vector[n_trials] trial_I_2  =  zTrial_I_2 * sigma_trial_2; 
  vector[n_groups] group_I_2  =  zGroup_I_2 * sigma_group_2;
  // Expected trait values for each observation
  vector[No] z_exp; // expected phenotypic values producing (eq1)
  vector[No] z2_exp; // expected phenotypic values scrounging (eq2)
  //Model equations:
  // Partition the phenotypes into DIEs and IIEs for producing and scrounging
  // Equation 1: Producing
  z_exp = mu + I[individual,1] + I[opponent1,2] + I[opponent2,2] + 
          trial_I[trialID] + group_I[groupID];
  // Equation 2: Scrounging
  z2_exp = mu_2 + I[individual,3] + I[opponent1,4] + I[opponent2,4] + 
           trial_I_2[trialID] + group_I_2[groupID];
}
model {
  // Priors
  // Fixed effects
  mu ~ normal(0, 1); // Weakly informative prior on the intercept
  // Random effects
  to_vector(Iz) ~ normal(0,1); // Implies i ~ normal(0,sigma_I)
  to_vector(zTrial_I) ~ normal(0, 1);
	to_vector(zGroup_I) ~ normal(0, 1);
  to_vector(zTrial_I_2) ~ normal(0, 1);
	to_vector(zGroup_I_2) ~ normal(0, 1);
  to_vector(sigma_I) ~ exponential(3); // Weakly informative prior for the variances
  to_vector(sigma_R) ~ exponential(3); 
  sigma_trial ~ exponential(3);
  sigma_group ~ exponential(3);
  sigma_trial_2 ~ exponential(3);
  sigma_group_2 ~ exponential(3);
  LI ~ lkj_corr_cholesky(1); //Prior individual level correlation
  LR ~ lkj_corr_cholesky(1); //Prior residual correlation
  // Residual correlation Choleksky factor
  matrix[2,2] L_sigma = diag_pre_multiply(sigma_R, LR);
  // Transform expected values to array
  array[No] vector[2] mus;
  for (o in 1:No) 
  mus[o] = [z_exp[o], z2_exp[o]]';
  // Cholesky multinormal likelihood function to estimate residual correlations
  zs ~ multi_normal_cholesky(mus, L_sigma);
}
generated quantities {
  // Variances individual effects
  real var_DIE = sigma_I[1]^2;
  real var_IIE = sigma_I[2]^2;
  real var_DIE_2 = sigma_I[3]^2;
  real var_IIE_2 = sigma_I[4]^2;
  // Variances other random effects
  real var_trialID = sigma_trial^2;
  real var_groupID = sigma_group^2;
  real var_trialID_2 = sigma_trial_2^2;
  real var_groupID_2 = sigma_group_2^2;
  // Residual variation
  real var_res = sigma_R[1]^2; 
  real var_res_2 = sigma_R[2]^2; 
  // Total variation 
  real<lower=0> var_P = var_DIE + (var_IIE*2) + var_trialID + var_groupID + var_res;
  real<lower=0> var_P_2 = var_DIE_2 + (var_IIE_2*2) + var_trialID_2 + var_groupID_2 +var_res_2;
  // Repeatabilities
  real<lower=0> rep_DIE = var_DIE/var_P; // rep DIE
  real<lower=0> rep_IIE = var_IIE/var_P; // rep IIE
  real<lower=0> rep_trialID = var_trialID/var_P; // rep trialID
  real<lower=0> rep_groupID = var_groupID/var_P; // rep groupID
  real<lower=0> rep_res = var_res/var_P; // rep residual
  real<lower=0> rep_DIE_2 = var_DIE_2/var_P_2; // rep DIE
  real<lower=0> rep_IIE_2 = var_IIE_2/var_P_2; // rep IIE
  real<lower=0> rep_trialID_2 = var_trialID_2/var_P_2; // rep trialID
  real<lower=0> rep_groupID_2 = var_groupID_2/var_P_2; // rep groupID
  real<lower=0> rep_res_2 = var_res_2/var_P_2; // rep residual
  // derive P matrix
  matrix[4,4] Omega_P = LI * LI'; // Correlation matrix
  matrix[4,4] D_I = diag_matrix(sigma_I); // Diagonal SD matrix
  matrix[4,4] Cov_I = D_I*Omega_P*D_I; // Covariance matrix  
  // Covariance
  real cov_P1 = Cov_I[1,2];
  real cov_P2 = Cov_I[1,3];
  real cov_P3 = Cov_I[1,4];
  real cov_P4 = Cov_I[2,3];
  real cov_P5 = Cov_I[2,4];
  real cov_P6 = Cov_I[3,4];
  // Correlations
  real cor_P1 = Omega_P[1,2];
  real cor_P2 = Omega_P[1,3];
  real cor_P3 = Omega_P[1,4];
  real cor_P4 = Omega_P[2,3];
  real cor_P5 = Omega_P[2,4];
  real cor_P6 = Omega_P[3,4];
  // Residual correlation
  matrix[2,2] res_cor_mat = LR * LR';
  real<lower=-1,upper=1> res_cor = res_cor_mat[1,2];
  real res_cov = res_cor_mat[1,2]*sqrt(var_res*var_res_2);
  // Total phenotypic effect
  real total_P_eff_1 = var_DIE + (4*cov_P1) + (4*var_IIE);
  real total_P_eff_2 = var_DIE_2 + (4*cov_P6) + (4*var_IIE_2);
  real rep_P_eff_1 = total_P_eff_1/var_P;
  real rep_P_eff_2 = total_P_eff_2/var_P_2;
}", file = "Models/Cholesky_bivar_PS_unadjusted_rep.stan")

# Fit the model
md_PS_bivar_full_rep_unadjusted <- stan("Models/Cholesky_bivar_PS_unadjusted_rep.stan", data = stan_bivar_data_rep_unadjusted, 
                                        chains = 5, iter = 5000,  
                                        warmup = 2000, thin = 1,
                                        save_warmup = F,
                                        cores = parallelly::availableCores()-1,
                                        seed = 2052728742
                                        )

beep(0)# Plays a sound when the model has been fit

# Derive seed if necessary to reproduce the same model output
rstan::get_seed(md_PS_bivar_full_rep_unadjusted)#2052728742

# Set the model parameters that need to be extracted from the model fit:
Pars_rep_unadjusted <- c("mu", "mu_2", 
                  "var_DIE", "var_IIE", "var_trialID", "var_groupID", "var_res","var_P",
                  "var_DIE_2", "var_IIE_2", "var_trialID_2", "var_groupID_2","var_res_2", "var_P_2",
                  "rep_DIE", "rep_IIE", "rep_trialID", "rep_groupID", "rep_res",
                  "rep_DIE_2", "rep_IIE_2", "rep_trialID_2", "rep_groupID_2","rep_res_2", 
                  "cov_P1", "cov_P2", "cov_P3", "cov_P4" , "cov_P5", "cov_P6",
                  "cor_P1", "cor_P2", "cor_P3", "cor_P4" , "cor_P5", "cor_P6", 
                  "res_cov", "res_cor",
                  "total_P_eff_1", "total_P_eff_2", "rep_P_eff_1", "rep_P_eff_2","lp__")

# Derive the model output from the posterior (Mean + CI + median + ESS + Rhat)
round(summary(md_PS_bivar_full_rep_unadjusted, pars = Pars_rep_unadjusted)$summary[,c(1,4,6,8,9,10)],3)
Sum_md_PS_bivar_full_rep_unadjusted <- round(summary(md_PS_bivar_full_rep_unadjusted, pars = Pars_rep_unadjusted)$summary[,c(1,4,6,8,9,10)],3)

# Look at posterior distribtions:
#launch_shinystan(md_PS_bivar_full)

# Extracting BLUPS & excpected values
#BLUPS_I_bivar_PS <- as.data.frame(rstan::extract(md_PS_bivar_full)[["I"]])
BLUPS_I_bivar_PS_rep_unadjusted <- round(summary(md_PS_bivar_full_rep_unadjusted, pars = "I")$summary[,c(1,4,6, 8, 9,10)],3)

# Save the model, summary & BLUPS
saveRDS(Sum_md_PS_bivar_full_rep_unadjusted, file = paste("Model_output/Sum_md_PS_bivar_full_unadjusted_rep.rds", sep = ""))
saveRDS(md_PS_bivar_full_rep_unadjusted, file = paste("Model_output/md_PS_bivar_full_unadjusted_rep.rds", sep = ""))
saveRDS(BLUPS_I_bivar_PS_rep_unadjusted, file = paste("Model_output/BLUPS_I_bivar_PS_unadjusted_rep.rds", sep = ""))

#===============================================================================

##### Cross-year repeatability model #####

# Create the dataset used in the model fit

# For producing
stan_data_short_rep_prd <- list(No = nrow(df),
                        Ni = length(unique(df$RingNR)),
                        individual = as.integer(as.factor(df$RingNR)),
                        opponent1 = as.integer(as.factor(df$Opponent1_ID)),
                        opponent2 = as.integer(as.factor(df$Opponent2_ID)),
                        trialID = as.integer(as.factor(df$TrialID)),
                        n_trials = length(unique(df$TrialID)),
                        groupID = as.integer(as.factor(df$GroupID)),
                        n_groups = length(unique(df$GroupID)),
                        trial_order = as.vector(scale(df$TrialNR_Indiv, center = T, scale = T)),
                        trialday = ifelse(df$Trial_day_bin == 0, -0.5, 0.5),
                        group_trial_b4 = df$Grp_tr_b4,
                        year = ifelse(df$Year == 2022, -0.5, 0.5),
                        z = as.vector(scale(df$Producing_events)),
                        series = as.integer(as.factor(df$Series_ID_year)),
                        n_series = length(unique(df$Series_ID_year)),
                        series_opp1 = as.integer(as.factor(df$Series_opp1_year)),
                        series_opp2 = as.integer(as.factor(df$Series_opp2_year)))

# For scrounging
stan_data_short_rep_scr <- list(No = nrow(df),
                                    Ni = length(unique(df$RingNR)),
                                    individual = as.integer(as.factor(df$RingNR)),
                                    opponent1 = as.integer(as.factor(df$Opponent1_ID)),
                                    opponent2 = as.integer(as.factor(df$Opponent2_ID)),
                                    trialID = as.integer(as.factor(df$TrialID)),
                                    n_trials = length(unique(df$TrialID)),
                                    groupID = as.integer(as.factor(df$GroupID)),
                                    n_groups = length(unique(df$GroupID)),
                                    trial_order = as.vector(scale(df$TrialNR_Indiv, center = T, scale = T)),
                                    trialday = ifelse(df$Trial_day_bin == 0, -0.5, 0.5),
                                    group_trial_b4 = df$Grp_tr_b4,
                                    year = ifelse(df$Year == 2022, -0.5, 0.5),
                                    z = as.vector(scale(df$Scrounging_events)),
                                    series = as.integer(as.factor(df$Series_ID_year)),
                                    n_series = length(unique(df$Series_ID_year)),
                                    series_opp1 = as.integer(as.factor(df$Series_opp1_year)),
                                    series_opp2 = as.integer(as.factor(df$Series_opp2_year)))

# The model code
write(
  temp <-   "data {
  // Number of clusters
   int<lower=0> No; // number of observations for phenotypes or number of rows
   int<lower=0> Ni; // number of individuals
   int<lower=0> n_trials; // number of trials
   int<lower=0> n_groups; // number of groups
   int<lower=0> n_series; // number of series (unique individual + year)
  // Clusters identifiers
   array[No] int<lower=1> individual;  //  Individual ID repeated obs
   array[No] int<lower=1> opponent1;  //  Individual ID opponent1 repeated obs
   array[No] int<lower=1> opponent2;  //  Individual ID opponent2 repeated obs
   array[No] int<lower=1> trialID;  //  Individual ID opponent repeated obs
   array[No] int<lower=1> groupID;  //  Individual ID opponent repeated obs
   array[No] int<lower=1> series;  //  Individual ID opponent repeated obs
   array[No] int<lower=1> series_opp1;  //  Individual ID opponent repeated obs
   array[No] int<lower=1> series_opp2;  //  Individual ID opponent repeated obs
  // Predictors
  vector[No] year; // Year 0/1 (2022/2023)for each observation
  vector[No] trial_order; // Trial order/sequence per individual for each observation
  vector[No] trialday; // Trialday 0/1 for each observation
  vector[No] group_trial_b4; // Group trial before 0/1 for each observation
  // Continuous variables
 	vector[No] z;  // phenotypic observations
  }
 transformed data {
// transform interactions into singular vector for faster interaction estimation
  vector[No] trial_order_trialday = trial_order .* trialday;  
  vector[No] trial_order_group_trial_b4 = trial_order .* group_trial_b4; 
  vector[No] group_trial_b4_trialday = group_trial_b4 .* trialday;  
  vector[No] trial_order_trialday_group_trial_b4 = trial_order .* trialday .* group_trial_b4;  
}
 parameters {
  // Fixed effects
 	real B_0; //intercept
  real B_year; // slope of year
  real B_trial_order;// slope of trial order
  real B_trialday; // slope of trialday
  real B_group_trial_b4; // slope of trialday
  real int_trialday_group_trial_b4_trialorder;
  real int_trialday_trialorder; 
  real int_trialday_group_trial_b4; 
  real int_group_trial_b4_trialorder;
  // Random effects
  matrix[Ni,2] zI; //(intercepts and opponent effect for each individual)
  vector<lower=0>[2] sigma_I; // sd  intercepts and slopes
  cholesky_factor_corr[2] LI;  // factor to estimate covariance int-slopes
  matrix[n_series,2] zIT; //(BLUPS time series)
  vector<lower=0>[2] sigma_IT; // sd  BLUPS
  cholesky_factor_corr[2] L_IT;  // factor to estimate covariance int-slopes  vector[n_trials] zTrial_I; //
  vector[n_groups] zGroup_I; //
  vector[n_trials] zTrial_I; //
  real <lower=0> sigma_e;
  real <lower=0> sigma_trial;
  real <lower=0> sigma_group;
 }
 transformed parameters{
	matrix[Ni,2] I  = zI * diag_pre_multiply(sigma_I, LI)'; //  Unscaled blups intercept and res_impact for each individual
  matrix[n_series,2] IT = zIT * diag_pre_multiply(sigma_IT, L_IT)'; // get the unscaled value  vector[n_trials] trial_I  =  zTrial_I * sigma_trial; // get the unscaled value
  vector[n_groups] group_I  =  zGroup_I * sigma_group; // get the unscaled value
  vector[n_trials] trial_I  =  zTrial_I * sigma_trial; // Get the unscaled values for random effect
 vector[No] e_z; // predicted values for phenotype
 // Model equation	
 e_z = B_0 + B_year*year + B_trialday * trialday + 
       B_trial_order * trial_order + B_group_trial_b4 * group_trial_b4 +
      int_trialday_trialorder *trial_order_trialday + 
      int_trialday_group_trial_b4 * group_trial_b4_trialday   +
      int_group_trial_b4_trialorder * trial_order_group_trial_b4 + 
      int_trialday_group_trial_b4_trialorder * trial_order_trialday_group_trial_b4 +
      I[individual,1] + I[opponent1,2] + I[opponent2,2] +
      IT[series,1] + IT[series_opp1,2] + IT[series_opp2,2] +
      trial_I[trialID] + group_I[groupID];
}
model {
 // Fixed effects prior distributions
 B_0 ~ normal(0,1); 
 B_year ~ normal(0,1); 
 B_trialday ~ normal(0,1); 
 B_trial_order ~ normal(0,1); 
 B_group_trial_b4 ~ normal(0,1); 
 int_trialday_trialorder ~ normal(0,1); 
 int_trialday_group_trial_b4 ~ normal(0,1); 
 int_group_trial_b4_trialorder ~ normal(0,1); 
 int_trialday_group_trial_b4_trialorder ~ normal(0,1); 
 // Random effects prior distribution
	to_vector(zI) ~ normal(0, 1);
  to_vector(zIT) ~ normal(0, 1);	
  to_vector(zTrial_I) ~ normal(0, 1);
	to_vector(zGroup_I) ~ normal(0, 1);
  to_vector(sigma_I) ~ exponential(3);
  to_vector(sigma_IT) ~ exponential(3);
  sigma_trial ~ exponential(3);
  sigma_group ~ exponential(3);
  sigma_e ~ exponential(3);
  LI ~ lkj_corr_cholesky(1);
  L_IT ~ lkj_corr_cholesky(1);
 // Likelihood function
 	z ~ normal(e_z, sigma_e);
}
generated quantities{
real<lower=0> Sigma2_intercept = sigma_I[1]^2;
real<lower=0> Sigma2_res_impact = sigma_I[2]^2;
real<lower=0> Sigma2_trialID = sigma_trial^2;  
real<lower=0> Sigma2_groupID = sigma_group^2;
real<lower=0> Sigma2_series = sigma_IT[1]^2;
real<lower=0> Sigma2_series_opp = sigma_IT[2]^2;
real<lower=0> Sigma2_res = sigma_e^2;
real<lower=0> Sigma2_total = Sigma2_intercept + (2*Sigma2_res_impact) + 
                             Sigma2_series + (2*Sigma2_series_opp) +
                             Sigma2_groupID + Sigma2_trialID + 
                             Sigma2_res;
matrix[2, 2]  Omega_I = LI * LI';
real cov_1 = Omega_I[1,2]*sqrt(Sigma2_res_impact*Sigma2_intercept); // DIE*IIE
real cor_1 = Omega_I[1,2];//    inte_res_impact
matrix[2, 2]  Omega_IT = L_IT * L_IT';
real cov_1_IT = Omega_IT[1,2]*sqrt(Sigma2_series_opp*Sigma2_series); // DIE*IIE IT
real cor_1_IT = Omega_IT[1,2];//    inte_res_impact
real var_comp_focal = Sigma2_intercept/Sigma2_total;
real var_comp_opponent = Sigma2_res_impact/Sigma2_total;
real var_comp_group = Sigma2_groupID/Sigma2_total;
real var_comp_trial = Sigma2_trialID/Sigma2_total;
real var_comp_res = Sigma2_res/Sigma2_total;
real Short_term_rep = (Sigma2_intercept + Sigma2_series)/Sigma2_total;
real Intercept_rep = Sigma2_intercept/(Sigma2_intercept + Sigma2_series);
real var_comp_series = Sigma2_series/Sigma2_total;
real Short_term_rep_opp = (Sigma2_res_impact + Sigma2_series_opp)/Sigma2_total;
real Intercept_rep_opp = Sigma2_res_impact/(Sigma2_res_impact + Sigma2_series_opp);
real var_comp_series_opp = Sigma2_series_opp/Sigma2_total;
}"
, file="Models/Short_rep_PS_full.stan")

### Producing ###

# Fit the model
md_short_rep_prd <- stan("Models/Short_rep_PS_full.stan", data = stan_data_short_rep_prd, 
                             chains = 5, iter = 5000,  
                             warmup = 2000, thin = 1, 
                             save_warmup = F,
                             cores = parallelly::availableCores()-1,
                             seed = 1503340858
                             )
beep(0)# Plays a sound when the model has been fit

# Derive seed if necessary to reproduce the same model output
rstan::get_seed(md_short_rep_prd)#1503340858

params_short_rep <- c("B_0", "B_year", "B_trial_order", "B_trialday", "B_group_trial_b4", 
                      "int_trialday_group_trial_b4", "int_trialday_trialorder","int_group_trial_b4_trialorder", "int_trialday_group_trial_b4_trialorder",  
                      "Sigma2_intercept", "Sigma2_res_impact", "Sigma2_trialID", "Sigma2_groupID", "Sigma2_series", "Sigma2_series_opp", "Sigma2_res", "Sigma2_total", 
                      "cov_1", "cor_1", "cov_1_IT", "cor_1_IT",
                      "var_comp_focal", "var_comp_opponent", 
                      "var_comp_group", "var_comp_trial", "var_comp_series", "var_comp_series_opp", "var_comp_res",
                      "Intercept_rep", "Short_term_rep", "Intercept_rep_opp", "Short_term_rep_opp", "lp__")

# Derive the model output from the posterior (Mean + CI + median + ESS + Rhat)
round(summary(md_short_rep_prd, pars = params_short_rep)$summary[,c(1,4,6,8,9,10)],3)
Sum_md_short_rep_prd <- round(summary(md_short_rep_prd, pars = params_short_rep)$summary[,c(1,4,6,8,9,10)],3)

# Extracting BLUPS
BLUPS_I_short_rep_prd <- round(summary(md_short_rep_prd, pars = "I")$summary[,c(1,4,6, 8, 9,10)],3)
BLUPS_IT_short_rep_prd <- round(summary(md_short_rep_prd, pars = "IT")$summary[,c(1,4,6, 8, 9,10)],3)

# Save the model, summary & BLUPS
saveRDS(Sum_md_short_rep_prd, file = paste("Model_output/Sum_md_short_rep_prd.rds", sep = ""))
saveRDS(md_short_rep_prd, file = paste("Model_output/md_short_rep_prd.rds", sep = ""))
saveRDS(BLUPS_I_short_rep_prd, file = paste("Model_output/BLUPS_I_short_rep_prd.rds", sep = ""))
saveRDS(BLUPS_IT_short_rep_prd, file = paste("Model_output/BLUPS_IT_short_rep_prd.rds", sep = ""))

### Scrounging ###

# Fit the model
md_short_rep_scr <- stan("Models/Short_rep_PS_full.stan", data = stan_data_short_rep_scr, 
                             chains = 5, iter = 5000,  
                             warmup = 2000, thin = 1, 
                             save_warmup = F,
                             cores = parallelly::availableCores()-1,
                             seed = 167006717
                             )

#2123140085
beep(0)# Plays a sound when the model has been fit

# Derive seed if necessary to reproduce the same model output
rstan::get_seed(md_short_rep_scr)#1961525646

# Derive the model output from the posterior (Mean + CI + median + ESS + Rhat)
round(summary(md_short_rep_scr, pars = params_short_rep)$summary[,c(1,4,6,8,9,10)],3)
Sum_md_short_rep_scr <- round(summary(md_short_rep_scr, pars = params_short_rep)$summary[,c(1,4,6,8,9,10)],3)

#launch_shinystan(md_short_rep_scr)

# Extracting BLUPS
BLUPS_I_short_rep_scr <- round(summary(md_short_rep_scr, pars = "I")$summary[,c(1,4,6, 8, 9,10)],3)
BLUPS_IT_short_rep_scr <- round(summary(md_short_rep_scr, pars = "IT")$summary[,c(1,4,6, 8, 9,10)],3)

# Save the model, summary & BLUPS
saveRDS(Sum_md_short_rep_scr, file = paste("Model_output/Sum_md_short_rep_scr.rds", sep = ""))
saveRDS(md_short_rep_scr, file = paste("Model_output/md_short_rep_scr.rds", sep = ""))
saveRDS(BLUPS_I_short_rep_scr, file = paste("Model_output/BLUPS_I_short_rep_scr.rds", sep = ""))
saveRDS(BLUPS_IT_short_rep_scr, file = paste("Model_output/BLUPS_IT_short_rep_scr.rds", sep = ""))

##### End of script ######
