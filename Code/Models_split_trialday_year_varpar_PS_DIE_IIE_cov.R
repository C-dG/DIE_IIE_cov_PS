#===============================================================================

# DIE-IIE variance partitioning IGE models
# Split per year & trial day

# Author: Corné de Groot
# Readme file for variable details

#===============================================================================

#install.packages("renv")
#renv::init()

library(rstan)
library(shinystan)
library(dplyr)
library(parallelly)
library(beepr)

# Package versions for reproducibility
packageVersion("rstan") #2.32.6
packageVersion("shinystan") #2.6.0
packageVersion("dplyr") #1.1.3
packageVersion("parallelly") #1.36.0
packageVersion("beepr") #1.3

R.Version() #"R version 4.3.1 

# Load in data 
#===============================================================================

# Set wd if loading the df does not work
#setwd(rstudioapi::getActiveProject())

df <- read.csv("Data/Data_PS_DIE_IIE_cov.csv", header = TRUE)

# Create a directory to save the split models
ifelse(file.exists("Model_output/Split_models"), "Dir already exists",
       dir.create("Model_output/Split_models"))

#===============================================================================
# Split up datasets

#2022
df_22 <- df %>% filter(Year == "2022")
#2022 Trialday 1
df_22_td1 <- df %>% filter(Year == 2022 & Trial_day_bin == 0)
#2022 Trialday 2
df_22_td2 <- df %>% filter(Year == 2022 & Trial_day_bin == 1)
#2023 
df_23 <- df %>% filter(Year == 2023)
#2023 Trialday 1
df_23_td1 <- df %>% filter(Year == 2023 & Trial_day_bin == 0)
#2023 Trialday 2
df_23_td2 <- df %>% filter(Year == 2023 & Trial_day_bin == 1)

#===============================================================================

# 2022 full

# Scrounging

params_22 <- c("B_0", "B_bmr","B_trialday", "B_trial_order", "Int_trialday_trialorder","Sigma2_intercept", "Sigma2_res_impact", 
               "Sigma2_trialID", "Sigma2_groupID",  "Sigma2_res", "Sigma2_total", "cov_1", "cor_1",
               "var_comp_focal", "var_comp_opponent", "var_comp_group", "var_comp_trial","var_comp_res")

stan_data_scr_22 = list(n_obs = nrow(df_22),
                        n_ind = length(unique(df_22$RingNR)),
                        individual = as.integer(as.factor(df_22$RingNR)),
                        opponent1 = as.integer(as.factor(df_22$Opponent1_ID)),
                        opponent2 = as.integer(as.factor(df_22$Opponent2_ID)),
                        trialID = as.integer(as.factor(df_22$TrialID)),
                        n_trials = length(unique(df_22$TrialID)),
                        groupID = as.integer(as.factor(df_22$GroupID)),
                        n_groups = length(unique(df_22$GroupID)),
                        trialday = ifelse(df_22$Trial_day_bin == 0, -0.5, 0.5),
                        trial_order = as.vector(scale(df_22$TrialNR_Indiv, center = T, scale = T)),
                        bmr = df_22$bmr_YN, 
                        z = as.vector(scale(df_22$Scrounging_events)))

stan_data_prd_22 = list(n_obs = nrow(df_22),
                        n_ind = length(unique(df_22$RingNR)),
                        individual = as.integer(as.factor(df_22$RingNR)),
                        opponent1 = as.integer(as.factor(df_22$Opponent1_ID)),
                        opponent2 = as.integer(as.factor(df_22$Opponent2_ID)),
                        trialID = as.integer(as.factor(df_22$TrialID)),
                        n_trials = length(unique(df_22$TrialID)),
                        groupID = as.integer(as.factor(df_22$GroupID)),
                        n_groups = length(unique(df_22$GroupID)),
                        trialday = ifelse(df_22$Trial_day_bin == 0, -0.5, 0.5),
                        trial_order = as.vector(scale(df_22$TrialNR_Indiv, center = T, scale = T)),
                        bmr = df_22$bmr_YN, 
                        z = as.vector(scale(df_22$Producing_events)))

write(
  temp <-   "data {
//Phenotype dataset
 
  // Number of clusters (an integer)
   int<lower=0> n_obs; // number of observations for phenotypes or number of rows
   int<lower=0> n_ind; // number of individuals

   int<lower=0> n_trials; // number of trials
   int<lower=0> n_groups; // number of groups

  // Clusters identifiers (an integer)
   int<lower=0> individual[n_obs];  //  Individual ID repeated obs
   int<lower=0> opponent1[n_obs];  //  Individual ID opponent1 repeated obs
   int<lower=0> opponent2[n_obs];  //  Individual ID opponent2 repeated obs
   
   int<lower=0> trialID[n_obs];  //  Individual ID opponent repeated obs
   int<lower=0> groupID[n_obs];  //  Individual ID opponent repeated obs

  //Predictors
  int<lower=0,upper=1> bmr[n_obs]; // fixed effect of bmr

  real trial_order[n_obs]; // fixed effect of trialday 1/2
  real trialday[n_obs]; // fixed effect of trial order/sequence per individual

  // Continuous variables
 	real z[n_obs];  // phenotypic observations
 
 }
 
 parameters {
// Define parameters to estimate
  // Parameters for model of phenotype
 	real B_0; //intercept

  real B_bmr; // beta for individuals that have been bmr'ed before the trials
  real B_trialday; // slope of trialday
  real B_trial_order;// slope of trial order
  real Int_trialday_trialorder;// 

   // Random effects
   matrix[2,n_ind]     	zI; //(intercepts and opponent effect for each individual)
   vector<lower=0>[2]  	sigma_I; // sd  intercepts and slopes
   cholesky_factor_corr[2] L;  // factor to estimate covariance int-slopes
   
   vector[n_trials] zTrial_I; //
   vector[n_groups] zGroup_I; //

   real <lower=0> sigma_e;
   
   real <lower=0> sigma_trial;
   real <lower=0> sigma_group;

 }
 
 transformed parameters{
	matrix[2,n_ind] I; //  Unscaled blups intercept and res_impact for each individual
	real e_z[n_obs]; // predicted values for phenotype
	
  vector[n_trials] trial_I; //
  vector[n_groups] group_I; //

  I  = diag_pre_multiply(sigma_I, L) * zI; // get the unscaled value
  
  trial_I  =  zTrial_I * sigma_trial; // get the unscaled value
  group_I  =  zGroup_I * sigma_group; // get the unscaled value

   for (i in 1:n_obs) {
 	e_z[i]  = B_0  + B_bmr*bmr[i] +B_trialday*trialday[i] + B_trial_order*trial_order[i] + 
  Int_trialday_trialorder*trialday[i]*trial_order[i] +  
  I[1, individual[i]] + I[2, opponent1[i]] + I[2, opponent2[i]] +  
  trial_I[trialID[i]] + group_I[groupID[i]];
   }
   
}
 
model {
// Create vector of predicted values
  B_0 ~ normal(0, 1); 

  B_trialday ~ normal(0, 1); 
  B_trial_order ~ normal(0, 1); 
  B_bmr ~ normal(0, 1); 
  Int_trialday_trialorder ~ normal(0, 1);

 // Random effects distribution
	to_vector(zI) ~ normal(0,1);

	to_vector(zTrial_I) ~ normal(0,1);
	to_vector(zGroup_I) ~ normal(0,1);

 to_vector(sigma_I) ~ normal(0,1);
 
 sigma_trial ~ exponential(3);
 sigma_group ~ exponential(3);
 sigma_e ~ exponential(3);

 L ~ lkj_corr_cholesky(1);
    
 // Likelihood function
	for (i in 1:n_obs)
 	z[i]~normal(e_z[i], sigma_e);
}

generated quantities{
real<lower=0> Sigma2_intercept;
real<lower=0> Sigma2_res_impact;
real<lower=0> Sigma2_trialID; // variance trialID
real<lower=0> Sigma2_groupID;

real<lower=0> Sigma2_res;
real<lower=0> Sigma2_total;

real cov_1; //   int_res_impact
real cor_1;//    inte_res_impact
matrix[2, 2]  Omega_I;

Sigma2_intercept = sigma_I[1]^2;
Sigma2_res_impact= sigma_I[2]^2;

Sigma2_groupID = sigma_group^2;
Sigma2_trialID = sigma_trial^2;

Sigma2_res = sigma_e^2;
Sigma2_total = Sigma2_intercept + Sigma2_res_impact + Sigma2_groupID + Sigma2_trialID + Sigma2_res;

Omega_I = L * L';
cov_1 = Omega_I[1,2]*sqrt(Sigma2_res_impact*Sigma2_intercept);

cor_1 = cov_1 /sqrt(Sigma2_intercept*Sigma2_res_impact);

// Estimate variance components of random effects

real var_comp_focal;
real var_comp_opponent;
real var_comp_group;
real var_comp_trial;
real var_comp_res;

var_comp_focal = Sigma2_intercept/Sigma2_total;
var_comp_opponent = Sigma2_res_impact/Sigma2_total ;
var_comp_group = Sigma2_groupID/Sigma2_total;
var_comp_trial = Sigma2_trialID/Sigma2_total;
var_comp_res = Sigma2_res/Sigma2_total;
}"

, file="Models/univ_22.stan")

#===============================================================================
# Scrounging

md_scr_22 <- stan("Models/univ_22.stan", data = stan_data_scr_22, 
                  pars = params_22, chains = 5, 
                  iter = 5000,  warmup = 2000, thin = 1, 
                  cores = parallel::detectCores()-1,
                  seed = 636313818
                  )
beep(0)

rstan::get_seed(md_scr_22)#636313818

# Mean + CI + median + ESS + Rhat
round(summary(md_scr_22)$summary[,c(1,4,6,8,9,10)],3)
Sum_md_scr_22 <- round(summary(md_scr_22)$summary[,c(1,4,6,8,9,10)],3)

#launch_shinystan(md_scr_22)
saveRDS(Sum_md_scr_22, file = paste("Model_output/Split_models/Sum_md_scr_22.rds", sep = ""))
saveRDS(md_scr_22, file = paste("Model_output/Split_models/md_scr_22.rds", sep = ""))

#===============================================================================
# Producing

md_prd_22 <- stan("Models/univ_22.stan", data = stan_data_prd_22, 
                  pars = params_22, chains = 5, 
                  iter = 5000,  warmup = 2000, 
                  thin = 1, cores = parallel::detectCores()-1,
                  seed = 1127420074
                  )
beep(0)

rstan::get_seed(md_prd_22)#1127420074

# Mean + CI + median + ESS + Rhat
round(summary(md_prd_22)$summary[,c(1,4,6,8,9,10)],3)
Sum_md_prd_22 <- round(summary(md_prd_22)$summary[,c(1,4,6,8,9,10)],3)

#launch_shinystan(md_prd_22)
saveRDS(Sum_md_prd_22, file = paste("Model_output/Split_models/Sum_md_prd_22.rds", sep = ""))
saveRDS(md_prd_22, file = paste("Model_output/Split_models/md_prd_22.rds", sep = ""))


#===============================================================================

# 2022 Trial day 1 & 2

params_22_td1 <- c("B_0", "B_bmr","B_trial_order","Sigma2_intercept", "Sigma2_res_impact", 
                   "Sigma2_trialID", "Sigma2_groupID", "Sigma2_res", "Sigma2_total", "cov_1", "cor_1",
                   "var_comp_focal", "var_comp_opponent", "var_comp_group", "var_comp_trial","var_comp_res")

stan_data_scr_22_td1 = list(n_obs = nrow(df_22_td1),
                            n_ind = length(unique(df_22_td1$RingNR)),
                            individual = as.integer(as.factor(df_22_td1$RingNR)),
                            opponent1 = as.integer(as.factor(df_22_td1$Opponent1_ID)),
                            opponent2 = as.integer(as.factor(df_22_td1$Opponent2_ID)),
                            trialID = as.integer(as.factor(df_22_td1$TrialID)),
                            n_trials = length(unique(df_22_td1$TrialID)),
                            groupID = as.integer(as.factor(df_22_td1$GroupID)),
                            n_groups = length(unique(df_22_td1$GroupID)),
                            trial_order = as.vector(scale(df_22_td1$TrialNR_Indiv, center = T, scale = T)),
                            bmr = df_22_td1$bmr_YN,
                            z = as.vector(scale(df_22_td1$Scrounging_events)))


stan_data_scr_22_td2 = list(n_obs = nrow(df_22_td2),
                            n_ind = length(unique(df_22_td2$RingNR)),
                            individual = as.integer(as.factor(df_22_td2$RingNR)),
                            opponent1 = as.integer(as.factor(df_22_td2$Opponent1_ID)),
                            opponent2 = as.integer(as.factor(df_22_td2$Opponent2_ID)),
                            trialID = as.integer(as.factor(df_22_td2$TrialID)),
                            n_trials = length(unique(df_22_td2$TrialID)),
                            groupID = as.integer(as.factor(df_22_td2$GroupID)),
                            n_groups = length(unique(df_22_td2$GroupID)),
                            trial_order = as.vector(scale(df_22_td2$TrialNR_Indiv, center = T, scale = T)),
                            bmr = df_22_td2$bmr_YN,
                            z = as.vector(scale(df_22_td2$Scrounging_events)))                           


stan_data_prd_22_td1 = list(n_obs = nrow(df_22_td1),
                            n_ind = length(unique(df_22_td1$RingNR)),
                            individual = as.integer(as.factor(df_22_td1$RingNR)),
                            opponent1 = as.integer(as.factor(df_22_td1$Opponent1_ID)),
                            opponent2 = as.integer(as.factor(df_22_td1$Opponent2_ID)),
                            trialID = as.integer(as.factor(df_22_td1$TrialID)),
                            n_trials = length(unique(df_22_td1$TrialID)),
                            groupID = as.integer(as.factor(df_22_td1$GroupID)),
                            n_groups = length(unique(df_22_td1$GroupID)),
                            trial_order = as.vector(scale(df_22_td1$TrialNR_Indiv, center = T, scale = T)),
                            bmr = df_22_td1$bmr_YN,
                            z = as.vector(scale(df_22_td1$Producing_events)))


stan_data_prd_22_td2 = list(n_obs = nrow(df_22_td2),
                            n_ind = length(unique(df_22_td2$RingNR)),
                            individual = as.integer(as.factor(df_22_td2$RingNR)),
                            opponent1 = as.integer(as.factor(df_22_td2$Opponent1_ID)),
                            opponent2 = as.integer(as.factor(df_22_td2$Opponent2_ID)),
                            trialID = as.integer(as.factor(df_22_td2$TrialID)),
                            n_trials = length(unique(df_22_td2$TrialID)),
                            groupID = as.integer(as.factor(df_22_td2$GroupID)),
                            n_groups = length(unique(df_22_td2$GroupID)),
                            trial_order = as.vector(scale(df_22_td2$TrialNR_Indiv, center = T, scale = T)),
                            bmr = df_22_td2$bmr_YN,
                            z = as.vector(scale(df_22_td2$Producing_events)))                          


write(
  temp <-   "data {
//Phenotype dataset
 
  // Number of clusters (an integer)
   int<lower=0> n_obs; // number of observations for phenotypes or number of rows
   int<lower=0> n_ind; // number of individuals

   int<lower=0> n_trials; // number of trials
   int<lower=0> n_groups; // number of groups

  // Clusters identifiers (an integer)
   int<lower=0> individual[n_obs];  //  Individual ID repeated obs
   int<lower=0> opponent1[n_obs];  //  Individual ID opponent1 repeated obs
   int<lower=0> opponent2[n_obs];  //  Individual ID opponent2 repeated obs
   
   int<lower=0> trialID[n_obs];  //  Individual ID opponent repeated obs
   int<lower=0> groupID[n_obs];  //  Individual ID opponent repeated obs

  //Predictors
  real trial_order[n_obs]; // fixed effect of trial order/sequence per individual
  int<lower=0,upper=1> bmr[n_obs]; // fixed effect of bmr

  // Continuous variables
 	real z[n_obs];  // phenotypic observations
 
 }
 
 parameters {
// Define parameters to estimate
  // Parameters for model of phenotype
 	real B_0; //intercept

  real B_trial_order;// slope of trial order
  real  B_bmr;//

   // Random effects
   matrix[2,n_ind]     	zI; //(intercepts and opponent effect for each individual)
   vector<lower=0>[2]  	sigma_I; // sd  intercepts and slopes
   cholesky_factor_corr[2] L;  // factor to estimate covariance int-slopes
   
   vector[n_trials] zTrial_I; //
   vector[n_groups] zGroup_I; //

   real <lower=0> sigma_e;
   
   real <lower=0> sigma_trial;
   real <lower=0> sigma_group;

 }
 
 transformed parameters{
	matrix[2,n_ind] I; //  Unscaled blups intercept and res_impact for each individual
	real e_z[n_obs]; // predicted values for phenotype
	
  vector[n_trials] trial_I; //
  vector[n_groups] group_I; //

  I  = diag_pre_multiply(sigma_I, L) * zI; // get the unscaled value
  
  trial_I  =  zTrial_I * sigma_trial; // get the unscaled value
  group_I  =  zGroup_I * sigma_group; // get the unscaled value

   for (i in 1:n_obs) {
 	e_z[i]  = B_0  + B_bmr*bmr[i] + B_trial_order*trial_order[i] +  I[1, individual[i]] + I[2, opponent1[i]] + I[2, opponent2[i]] +  trial_I[trialID[i]] + group_I[groupID[i]];
   }
   
}
 
model {
// Create vector of predicted values
 B_0 ~ normal(0, 1); 

  B_trial_order ~ normal(0, 1);
  B_bmr ~ normal(0, 1);

 // Random effects distribution
	to_vector(zI) ~ normal(0,1);

	to_vector(zTrial_I) ~ normal(0,1);
	to_vector(zGroup_I) ~ normal(0,1);

 to_vector(sigma_I) ~ normal(0,1);
 
 sigma_trial ~ exponential(3);
 sigma_group ~ exponential(3);
 sigma_e ~exponential(3);

 L ~ lkj_corr_cholesky(1);
    
 // Likelihood function
	for (i in 1:n_obs)
 	z[i]~normal(e_z[i], sigma_e);
}

generated quantities{
real<lower=0> Sigma2_intercept;
real<lower=0> Sigma2_res_impact;
real<lower=0> Sigma2_trialID; // variance trialID
real<lower=0> Sigma2_groupID;

real<lower=0> Sigma2_res;
real<lower=0> Sigma2_total;

real cov_1; //   int_res_impact
real cor_1;//    inte_res_impact
matrix[2, 2]  Omega_I;

Sigma2_intercept = sigma_I[1]^2;
Sigma2_res_impact= sigma_I[2]^2;

Sigma2_groupID = sigma_group^2;
Sigma2_trialID = sigma_trial^2;

Sigma2_res = sigma_e^2;
Sigma2_total = Sigma2_intercept + Sigma2_res_impact + Sigma2_groupID + Sigma2_trialID + Sigma2_res ;

Omega_I = L * L';
cov_1 = Omega_I[1,2]*sqrt(Sigma2_res_impact*Sigma2_intercept);

cor_1 = cov_1 /sqrt(Sigma2_intercept*Sigma2_res_impact);

// Estimate variance components of random effects

real var_comp_focal;
real var_comp_opponent;
real var_comp_group;
real var_comp_trial;
real var_comp_res;

var_comp_focal = Sigma2_intercept/Sigma2_total;
var_comp_opponent = Sigma2_res_impact/Sigma2_total ;
var_comp_group = Sigma2_groupID/Sigma2_total;
var_comp_trial = Sigma2_trialID/Sigma2_total;
var_comp_res = Sigma2_res/Sigma2_total;


}"

, file="Models/univ_22_td1.stan")

#===============================================================================
#  Scrounging 

md_scr_22_td1 <- stan("Models/univ_22_td1.stan", data = stan_data_scr_22_td1, 
                      pars = params_22_td1, chains = 5, 
                      iter = 5000,  warmup = 2000, 
                      thin = 1, cores = parallel::detectCores()-1,
                      seed = 422181711
                        )
beep(0)

rstan::get_seed(md_scr_22_td1)#422181711

# Mean + CI + median + ESS + Rhat
round(summary(md_scr_22_td1)$summary[,c(1,4,6,8,9,10)],3)
Sum_md_scr_22_td1 <- round(summary(md_scr_22_td1)$summary[,c(1,4,6,8,9,10)],3)

#launch_shinystan(md_scr_22_td1)
saveRDS(Sum_md_scr_22_td1, file = paste("Model_output/Split_models/Sum_md_scr_22_td1.rds", sep = ""))
saveRDS(md_scr_22_td1, file = paste("Model_output/Split_models/md_scr_22_td1.rds", sep = ""))


md_scr_22_td2 <- stan("Models/univ_22_td1.stan", data = stan_data_scr_22_td2, 
                      pars = params_22_td1, chains = 5, 
                      iter = 5000,  warmup = 2000, 
                      thin = 1, cores = parallel::detectCores()-1,
                      seed = 844629292
                        )
beep(0)

rstan::get_seed(md_scr_22_td2)#844629292

# Mean + CI + median + ESS + Rhat
round(summary(md_scr_22_td2)$summary[,c(1,4,6,8,9,10)],3)
Sum_md_scr_22_td2 <- round(summary(md_scr_22_td2)$summary[,c(1,4,6,8,9,10)],3)

#launch_shinystan(md_scr_22_td2)
saveRDS(Sum_md_scr_22_td2, file = paste("Model_output/Split_models/Sum_md_scr_22_td2.rds", sep = ""))
saveRDS(md_scr_22_td2, file = paste("Model_output/Split_models/md_scr_22_td2.rds", sep = ""))


#===============================================================================
# Producing

md_prd_22_td1 <- stan("Models/univ_22_td1.stan", data = stan_data_prd_22_td1, 
                      pars = params_22_td1, 
                      chains = 5, iter = 5000,  
                      warmup = 2000, thin = 1, 
                      cores = parallel::detectCores()-1,
                      seed = 899696659
                        )
beep(0)

rstan::get_seed(md_prd_22_td1)#899696659

# Mean + CI + median + ESS + Rhat
round(summary(md_prd_22_td1)$summary[,c(1,4,6,8,9,10)],3)
Sum_md_prd_22_td1 <- round(summary(md_prd_22_td1)$summary[,c(1,4,6,8,9,10)],3)

#launch_shinystan(md_prd_22_td1)
saveRDS(Sum_md_prd_22_td1, file = paste("Model_output/Split_models/Sum_md_prd_22_td1.rds", sep = ""))
saveRDS(md_prd_22_td1, file = paste("Model_output/Split_models/md_prd_22_td1.rds", sep = ""))


md_prd_22_td2 <- stan("Models/univ_22_td1.stan", data = stan_data_prd_22_td2, 
                      pars = params_22_td1, chains = 5, 
                      iter = 5000,  warmup = 2000, 
                      thin = 1, cores = parallel::detectCores()-1,
                      seed = 606701846
                        )
beep(0)

rstan::get_seed(md_prd_22_td2)#606701846

# Mean + CI + median + ESS + Rhat
round(summary(md_prd_22_td2)$summary[,c(1,4,6,8,9,10)],3)
Sum_md_prd_22_td2 <- round(summary(md_prd_22_td2)$summary[,c(1,4,6,8,9,10)],3)

#launch_shinystan(md_prd_22_td2)
saveRDS(Sum_md_prd_22_td2, file = paste("Model_output/Split_models/Sum_md_prd_22_td2.rds", sep = ""))
saveRDS(md_prd_22_td2, file = paste("Model_output/Split_models/md_prd_22_td2.rds", sep = ""))


#===============================================================================

#2023

stan_data_scr_23 =  list(n_obs = nrow(df_23),
                         n_ind = length(unique(df_23$RingNR)),
                         individual = as.integer(as.factor(df_23$RingNR)),
                         opponent1 = as.integer(as.factor(df_23$Opponent1_ID)),
                         opponent2 = as.integer(as.factor(df_23$Opponent2_ID)),
                         trialID = as.integer(as.factor(df_23$TrialID)),
                         n_trials = length(unique(df_23$TrialID)),
                         groupID = as.integer(as.factor(df_23$GroupID)),
                         n_groups = length(unique(df_23$GroupID)),
                         trialday = ifelse(df_23$Trial_day_bin == 0, -0.5, 0.5),                         trial_order = as.vector(scale(df_23$TrialNR_Indiv, center = T, scale = T)),
                         group_trial_b4 = df_23$Grp_tr_b4,
                         z = as.vector(scale(df_23$Scrounging_events)))

stan_data_prd_23 =  list(n_obs = nrow(df_23),
                         n_ind = length(unique(df_23$RingNR)),
                         individual = as.integer(as.factor(df_23$RingNR)),
                         opponent1 = as.integer(as.factor(df_23$Opponent1_ID)),
                         opponent2 = as.integer(as.factor(df_23$Opponent2_ID)),
                         trialID = as.integer(as.factor(df_23$TrialID)),
                         n_trials = length(unique(df_23$TrialID)),
                         groupID = as.integer(as.factor(df_23$GroupID)),
                         n_groups = length(unique(df_23$GroupID)),
                         trial_order = as.vector(scale(df_23$TrialNR_Indiv, center = T, scale = T)),
                         trialday = ifelse(df_23$Trial_day_bin == 0, -0.5, 0.5),                         trial_order = as.vector(scale(df_23$TrialNR_Indiv, center = T, scale = T)),
                         group_trial_b4 = df_23$Grp_tr_b4,
                         z = as.vector(scale(df_23$Producing_events)))

params_23 <- c("B_0", "B_trial_order", "B_trialday", "B_group_trial_b4", "int_trialday_trialorder", "int_trialday_group_trial_b4", "int_group_trial_b4_trialorder", "int_trialday_group_trial_b4_trialorder", "Sigma2_intercept", "Sigma2_res_impact", 
                 "Sigma2_trialID", "Sigma2_groupID",  "Sigma2_res", "Sigma2_total", "cov_1", "cor_1",
                 "var_comp_focal", "var_comp_opponent", "var_comp_group", "var_comp_trial","var_comp_res")


write(
  temp <-   "data {
//Phenotype dataset
 
  // Number of clusters (an integer)
   int<lower=0> n_obs; // number of observations for phenotypes or number of rows
   int<lower=0> n_ind; // number of individuals

   int<lower=0> n_trials; // number of trials
   int<lower=0> n_groups; // number of groups

  // Clusters identifiers (an integer)
   int<lower=0> individual[n_obs];  //  Individual ID repeated obs
   int<lower=0> opponent1[n_obs];  //  Individual ID opponent1 repeated obs
   int<lower=0> opponent2[n_obs];  //  Individual ID opponent2 repeated obs
   
   int<lower=0> trialID[n_obs];  //  Individual ID opponent repeated obs
   int<lower=0> groupID[n_obs];  //  Individual ID opponent repeated obs

  //Predictors
  real trialday[n_obs]; // fixed effect of trial order/sequence per individual

  real trial_order[n_obs]; // fixed effect of trialday 1/2
  int<lower=0,upper=1> group_trial_b4[n_obs]; // fixed effect of trialday 1/2

  // Continuous variables
 	real z[n_obs];  // phenotypic observations
 
 }
 
 parameters {
// Define parameters to estimate
  // Parameters for model of phenotype
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
   matrix[2,n_ind]     	zI; //(intercepts and opponent effect for each individual)
   vector<lower=0>[2]  	sigma_I; // sd  intercepts and slopes
   cholesky_factor_corr[2] L;  // factor to estimate covariance int-slopes
   
   vector[n_trials] zTrial_I; //
   vector[n_groups] zGroup_I; //

   real <lower=0> sigma_e;
   
   real <lower=0> sigma_trial;
   real <lower=0> sigma_group;

 }
 
 transformed parameters{
	matrix[2,n_ind] I; //  Unscaled blups intercept and res_impact for each individual
	real e_z[n_obs]; // predicted values for phenotype
	
  vector[n_trials] trial_I; //
  vector[n_groups] group_I; //

  I  = diag_pre_multiply(sigma_I, L) * zI; // get the unscaled value
  
  trial_I  =  zTrial_I * sigma_trial; // get the unscaled value
  group_I  =  zGroup_I * sigma_group; // get the unscaled value

   for (i in 1:n_obs) {
 	e_z[i]  = B_0  + B_trialday*trialday[i] + B_trial_order*trial_order[i] + B_group_trial_b4*group_trial_b4[i] + int_trialday_trialorder*trialday[i]*trial_order[i] + int_trialday_group_trial_b4*trialday[i]*group_trial_b4[i] + int_group_trial_b4_trialorder*group_trial_b4[i]*trial_order[i] + int_trialday_group_trial_b4_trialorder*trialday[i]*group_trial_b4[i]*trial_order[i] + I[1, individual[i]] + I[2, opponent1[i]] + I[2, opponent2[i]] + trial_I[trialID[i]] + group_I[groupID[i]];
   }
   
}
 
model {
// Create vector of predicted values
 B_0 ~ normal(0, 1); 

  B_year ~ normal(0, 1);
  B_trial_order ~ normal(0, 1);

  B_trialday ~ normal(0, 1);
  B_group_trial_b4 ~ normal(0, 1);

  int_trialday_group_trial_b4_trialorder ~ normal(0, 1);

  int_trialday_trialorder ~ normal(0, 1);
  int_trialday_group_trial_b4 ~ normal(0, 1);
  int_group_trial_b4_trialorder ~ normal(0, 1);

 // Random effects distribution
	to_vector(zI) ~ normal(0,1);

	to_vector(zTrial_I) ~ normal(0,1);
	to_vector(zGroup_I) ~ normal(0,1);

 to_vector(sigma_I) ~ normal(0,1);
 
 sigma_trial ~ exponential(3);
 sigma_group ~ exponential(3);
 sigma_e ~exponential(3);

 L ~ lkj_corr_cholesky(1);
    
    
 // Likelihood function
	for (i in 1:n_obs)
 	z[i]~normal(e_z[i], sigma_e);
}

generated quantities{
real<lower=0> Sigma2_intercept;
real<lower=0> Sigma2_res_impact;

real<lower=0> Sigma2_trialID; // variance trialID
real<lower=0> Sigma2_groupID;

real<lower=0> Sigma2_res;
real<lower=0> Sigma2_total;

real cov_1; //   int_res_impact
real cor_1;//    int_res_impact
matrix[2, 2]  Omega_I;

Sigma2_intercept = sigma_I[1]^2;
Sigma2_res_impact= sigma_I[2]^2;

Sigma2_groupID = sigma_group^2;
Sigma2_trialID = sigma_trial^2;


Sigma2_res = sigma_e^2;

Sigma2_total = Sigma2_intercept + Sigma2_res_impact + Sigma2_groupID + Sigma2_trialID + Sigma2_res ;

Omega_I = L * L';
cov_1 = Omega_I[1,2]*sqrt(Sigma2_res_impact*Sigma2_intercept);

cor_1 = cov_1 /sqrt(Sigma2_intercept*Sigma2_res_impact);

// Estimate variance components of random effects

real var_comp_focal;
real var_comp_opponent;
real var_comp_group;
real var_comp_trial;
real var_comp_res;

var_comp_focal = Sigma2_intercept/Sigma2_total;
var_comp_opponent = Sigma2_res_impact/Sigma2_total ;
var_comp_group = Sigma2_groupID/Sigma2_total;
var_comp_trial = Sigma2_trialID/Sigma2_total;
var_comp_res = Sigma2_res/Sigma2_total;

}"

, file="Models/univ_23.stan")


#===============================================================================
#  Scrounging 

md_scr_23 <- stan("Models/univ_23.stan", data = stan_data_scr_23, 
                    pars = params_23, chains = 5, 
                    iter = 5000,  warmup = 2000, 
                    thin = 1, cores = parallel::detectCores()-1,
                    seed = 1740749604
                      )
beep(0)

rstan::get_seed(md_scr_23)#1740749604

# Mean + CI + median + ESS + Rhat
round(summary(md_scr_23)$summary[,c(1,4,6,8,9,10)],3)
Sum_md_scr_23 <- round(summary(md_scr_23)$summary[,c(1,4,6,8,9,10)],3)

#launch_shinystan(md_scr_23)
saveRDS(Sum_md_scr_23, file = paste("Model_output/Split_models/Sum_md_scr_23.rds", sep = ""))
saveRDS(md_scr_23, file = paste("Model_output/Split_models/md_scr_23.rds", sep = ""))

#===============================================================================
# Producing

md_prd_23 <- stan("Models/univ_23.stan", data = stan_data_prd_23, 
                    pars = params_23, chains = 5, 
                    iter = 5000,  warmup = 2000, thin = 1, 
                    cores = parallel::detectCores()-1,
                    seed = 445948041)
beep(0)

rstan::get_seed(md_prd_23)#445948041

# Mean + CI + median + ESS + Rhat
round(summary(md_prd_23)$summary[,c(1,4,6,8,9,10)],3)
Sum_md_prd_23 <- round(summary(md_prd_23)$summary[,c(1,4,6,8,9,10)],3)

#launch_shinystan(md_prd_23)
saveRDS(Sum_md_prd_23, file = "Model_output/Split_models/Sum_md_prd_23.rds")
saveRDS(md_prd_23, file = "Model_output/Split_models/md_prd_23.rds")


#===============================================================================

#2023 Trialday 1 & 2


params_23_td1 <- c("B_0", "B_trial_order", "B_group_trial_b4", "int_trialorder_group_trial_b4", "Sigma2_intercept", "Sigma2_res_impact", 
                   "Sigma2_trialID", "Sigma2_groupID",  "Sigma2_res", "Sigma2_total", "cov_1", "cor_1",
                   "var_comp_focal", "var_comp_opponent", "var_comp_group", "var_comp_trial","var_comp_res")

stan_data_scr_23_td1 =  list(n_obs = nrow(df_23_td1),
                             n_ind = length(unique(df_23_td1$RingNR)),
                             individual = as.integer(as.factor(df_23_td1$RingNR)),
                             opponent1 = as.integer(as.factor(df_23_td1$Opponent1_ID)),
                             opponent2 = as.integer(as.factor(df_23_td1$Opponent2_ID)),
                             trialID = as.integer(as.factor(df_23_td1$TrialID)),
                             n_trials = length(unique(df_23_td1$TrialID)),
                             groupID = as.integer(as.factor(df_23_td1$GroupID)),
                             n_groups = length(unique(df_23_td1$GroupID)),
                             trial_order = as.vector(scale(df_23_td1$TrialNR_Indiv, center = T, scale = T)),
                             group_trial_b4 = df_23_td1$Grp_tr_b4,
                             z = as.vector(scale(df_23_td1$Scrounging_events)))


stan_data_scr_23_td2 =  list(n_obs = nrow(df_23_td2),
                             n_ind = length(unique(df_23_td2$RingNR)),
                             individual = as.integer(as.factor(df_23_td2$RingNR)),
                             opponent1 = as.integer(as.factor(df_23_td2$Opponent1_ID)),
                             opponent2 = as.integer(as.factor(df_23_td2$Opponent2_ID)),
                             trialID = as.integer(as.factor(df_23_td2$TrialID)),
                             n_trials = length(unique(df_23_td2$TrialID)),
                             groupID = as.integer(as.factor(df_23_td2$GroupID)),
                             n_groups = length(unique(df_23_td2$GroupID)),
                             trial_order = as.vector(scale(df_23_td2$TrialNR_Indiv, center = T, scale = T)),
                             group_trial_b4 = df_23_td2$Grp_tr_b4,
                             z = as.vector(scale(df_23_td2$Scrounging_events)))

stan_data_prd_23_td1 =  list(n_obs = nrow(df_23_td1),
                             n_ind = length(unique(df_23_td1$RingNR)),
                             individual = as.integer(as.factor(df_23_td1$RingNR)),
                             opponent1 = as.integer(as.factor(df_23_td1$Opponent1_ID)),
                             opponent2 = as.integer(as.factor(df_23_td1$Opponent2_ID)),
                             trialID = as.integer(as.factor(df_23_td1$TrialID)),
                             n_trials = length(unique(df_23_td1$TrialID)),
                             groupID = as.integer(as.factor(df_23_td1$GroupID)),
                             n_groups = length(unique(df_23_td1$GroupID)),
                             trial_order = as.vector(scale(df_23_td1$TrialNR_Indiv, center = T, scale = T)),
                             group_trial_b4 = df_23_td1$Grp_tr_b4,
                             z = as.vector(scale(df_23_td1$Producing_events)))

stan_data_prd_23_td2 =  list(n_obs = nrow(df_23_td2),
                             n_ind = length(unique(df_23_td2$RingNR)),
                             individual = as.integer(as.factor(df_23_td2$RingNR)),
                             opponent1 = as.integer(as.factor(df_23_td2$Opponent1_ID)),
                             opponent2 = as.integer(as.factor(df_23_td2$Opponent2_ID)),
                             trialID = as.integer(as.factor(df_23_td2$TrialID)),
                             n_trials = length(unique(df_23_td2$TrialID)),
                             groupID = as.integer(as.factor(df_23_td2$GroupID)),
                             n_groups = length(unique(df_23_td2$GroupID)),
                             trial_order = as.vector(scale(df_23_td2$TrialNR_Indiv, center = T, scale = T)),
                             group_trial_b4 = df_23_td2$Grp_tr_b4,
                             z = as.vector(scale(df_23_td2$Producing_events)))

write(
  temp <-   "data {
//Phenotype dataset
 
  // Number of clusters (an integer)
   int<lower=0> n_obs; // number of observations for phenotypes or number of rows
   int<lower=0> n_ind; // number of individuals

   int<lower=0> n_trials; // number of trials
   int<lower=0> n_groups; // number of groups

  // Clusters identifiers (an integer)
   int<lower=0> individual[n_obs];  //  Individual ID repeated obs
   int<lower=0> opponent1[n_obs];  //  Individual ID opponent1 repeated obs
   int<lower=0> opponent2[n_obs];  //  Individual ID opponent2 repeated obs
   
   int<lower=0> trialID[n_obs];  //  Individual ID opponent repeated obs
   int<lower=0> groupID[n_obs];  //  Individual ID opponent repeated obs

  //Predictors
  real trial_order[n_obs]; // fixed effect of trial order/sequence per individual

  int<lower=0,upper=1> group_trial_b4[n_obs]; // fixed effect of trialday 1/2

  // Continuous variables
 	real z[n_obs];  // phenotypic observations
 
 }
 
 parameters {
// Define parameters to estimate
  // Parameters for model of phenotype
 	real B_0; //intercept

  real B_trial_order;// slope of trial order

  real B_group_trial_b4; // slope of trialday

  real int_trialorder_group_trial_b4;

   // Random effects
   matrix[2,n_ind]     	zI; //(intercepts and opponent effect for each individual)
   vector<lower=0>[2]  	sigma_I; // sd  intercepts and slopes
   cholesky_factor_corr[2] L;  // factor to estimate covariance int-slopes
   
   vector[n_trials] zTrial_I; //
   vector[n_groups] zGroup_I; //

   real <lower=0> sigma_e;
   
   real <lower=0> sigma_trial;
   real <lower=0> sigma_group;

 }
 
 transformed parameters{
	matrix[2,n_ind] I; //  Unscaled blups intercept and res_impact for each individual
	real e_z[n_obs]; // predicted values for phenotype
	
  vector[n_trials] trial_I; //
  vector[n_groups] group_I; //

  I  = diag_pre_multiply(sigma_I, L) * zI; // get the unscaled value
  
  trial_I  =  zTrial_I * sigma_trial; // get the unscaled value
  group_I  =  zGroup_I * sigma_group; // get the unscaled value

   for (i in 1:n_obs) {
 	e_z[i]  = B_0  + B_trial_order*trial_order[i] + B_group_trial_b4*group_trial_b4[i] + int_trialorder_group_trial_b4*trial_order[i]*group_trial_b4[i] + I[1, individual[i]] + I[2, opponent1[i]] + I[2, opponent2[i]] + trial_I[trialID[i]] + group_I[groupID[i]];
   }
   
}
 
model {
// Create vector of predicted values
 B_0 ~ normal(0, 1); 

  B_trial_order ~ normal(0, 1);

  B_group_trial_b4 ~ normal(0, 1);

  int_trialorder_group_trial_b4 ~ normal(0, 1);

 // Random effects distribution
	to_vector(zI) ~ normal(0,1);

	to_vector(zTrial_I) ~ normal(0,1);
	to_vector(zGroup_I) ~ normal(0,1);

 to_vector(sigma_I) ~ normal(0,1);
 
 sigma_trial ~ exponential(3);
 sigma_group ~ exponential(3);
 sigma_e ~exponential(3);

 L ~ lkj_corr_cholesky(1);
    
    
 // Likelihood function
	for (i in 1:n_obs)
 	z[i]~normal(e_z[i], sigma_e);
}

generated quantities{
real<lower=0> Sigma2_intercept;
real<lower=0> Sigma2_res_impact;

real<lower=0> Sigma2_trialID; // variance trialID
real<lower=0> Sigma2_groupID;

real<lower=0> Sigma2_res;
real<lower=0> Sigma2_total;

real cov_1; //   int_res_impact
real cor_1;//    int_res_impact
matrix[2, 2]  Omega_I;

Sigma2_intercept = sigma_I[1]^2;
Sigma2_res_impact= sigma_I[2]^2;

Sigma2_groupID = sigma_group^2;
Sigma2_trialID = sigma_trial^2;

Sigma2_res = sigma_e^2;

Sigma2_total = Sigma2_intercept + Sigma2_res_impact + Sigma2_groupID + Sigma2_trialID + Sigma2_res ;

Omega_I = L * L';
cov_1 = Omega_I[1,2]*sqrt(Sigma2_res_impact*Sigma2_intercept);

cor_1 = cov_1 /sqrt(Sigma2_intercept*Sigma2_res_impact);

// Estimate variance components of random effects

real var_comp_focal;
real var_comp_opponent;
real var_comp_group;
real var_comp_trial;
real var_comp_res;

var_comp_focal = Sigma2_intercept/Sigma2_total;
var_comp_opponent = Sigma2_res_impact/Sigma2_total ;
var_comp_group = Sigma2_groupID/Sigma2_total;
var_comp_trial = Sigma2_trialID/Sigma2_total;
var_comp_res = Sigma2_res/Sigma2_total;

}"

, file="Models/univ_23_td1.stan")


#===============================================================================
#  Scrounging 

md_scr_23_td1 <- stan("Models/univ_23_td1.stan", data = stan_data_scr_23_td1, 
                      pars = params_23_td1, chains = 5, 
                      iter = 5000,  warmup = 2000, 
                      thin = 1, cores = parallel::detectCores()-1,
                      seed = 1844982079
                        )
beep(0)

rstan::get_seed(md_scr_23_td1)#1844982079

# Mean + CI + median + ESS + Rhat
round(summary(md_scr_23_td1)$summary[,c(1,4,6,8,9,10)],3)
Sum_md_scr_23_td1 <- round(summary(md_scr_23_td1)$summary[,c(1,4,6,8,9,10)],3)

#launch_shinystan(md_scr_23_td1)
saveRDS(Sum_md_scr_23_td1, file = paste("Model_output/Split_models/Sum_md_scr_23_td1.rds", sep = ""))
saveRDS(md_scr_23_td1, file = paste("Model_output/Split_models/md_scr_23_td1.rds", sep = ""))


md_scr_23_td2 <- stan("Models/univ_23_td1.stan", data = stan_data_scr_23_td2, 
                      pars = params_23_td1, chains = 5, 
                      iter = 5000,  warmup = 2000, 
                      thin = 1, cores = parallel::detectCores()-1,
                      seed = 1356399243
                      )
beep(0)

rstan::get_seed(md_scr_23_td2)#1140292252

# Mean + CI + median + ESS + Rhat
round(summary(md_scr_23_td2)$summary[,c(1,4,6,8,9,10)],3)
Sum_md_scr_23_td2 <- round(summary(md_scr_23_td2)$summary[,c(1,4,6,8,9,10)],3)

#launch_shinystan(md_scr_23_td2)
saveRDS(Sum_md_scr_23_td2, file = paste("Model_output/Split_models/Sum_md_scr_23_td2.rds", sep = ""))
saveRDS(md_scr_23_td2, file = paste("Model_output/Split_models/md_scr_23_td2.rds", sep = ""))


#===============================================================================
# Producing

md_prd_23_td1 <- stan("Models/univ_23_td1.stan", data = stan_data_prd_23_td1, 
                      pars = params_23_td1, chains = 5, 
                      iter = 5000,  warmup = 2000, 
                      thin = 1, cores = parallel::detectCores()-1,
                      seed = 2026792821
                        )
beep(0)

rstan::get_seed(md_prd_23_td1)#2026792821

# Mean + CI + median + ESS + Rhat
round(summary(md_prd_23_td1)$summary[,c(1,4,6,8,9,10)],3)
Sum_md_prd_23_td1 <- round(summary(md_prd_23_td1)$summary[,c(1,4,6,8,9,10)],3)

#launch_shinystan(md_prd_23_td1)
saveRDS(Sum_md_prd_23_td1, file = "Model_output/Split_models/Sum_md_prd_23_td1.rds")
saveRDS(md_prd_23_td1, file = "Model_output/Split_models/md_prd_23_td1.rds")


md_prd_23_td2 <- stan("Models/univ_23_td1.stan", data = stan_data_prd_23_td2, 
                      pars = params_23_td1, chains = 5, 
                      iter = 5000,  warmup = 2000, 
                      thin = 1, cores = parallel::detectCores()-1,
                      seed = 385885096
                        )
beep(0)

rstan::get_seed(md_prd_23_td2)#385885096

# Mean + CI + median + ESS + Rhat
round(summary(md_prd_23_td2)$summary[,c(1,4,6,8,9,10)],3)
Sum_md_prd_23_td2 <- round(summary(md_prd_23_td2)$summary[,c(1,4,6,8,9,10)],3)

#launch_shinystan(md_prd_23_td2)
saveRDS(Sum_md_prd_23_td2, file = paste("Model_output/Split_models/Sum_md_prd_23_td2.rds", sep = ""))
saveRDS(md_prd_23_td2, file = paste("Model_output/Split_models/md_prd_23_td2.rds", sep = ""))


#===============================================================================


# Full 2022 + 2023 model

params_full <- c("B_0", "B_year", "B_trial_order", "B_trialday", "B_group_trial_b4", "int_trialday_trialorder", "int_trialday_group_trial_b4", "int_group_trial_b4_trialorder", "int_trialday_group_trial_b4_trialorder", "Sigma2_intercept", "Sigma2_res_impact", 
                  "Sigma2_trialID", "Sigma2_groupID",  "Sigma2_res", "Sigma2_total", "cov_1", "cor_1",
                  "var_comp_focal", "var_comp_opponent", "var_comp_group", "var_comp_trial","var_comp_res")


stan_data_full_scr =  list(n_obs = nrow(df),
                                 n_ind = length(unique(df$RingNR)),
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
                                 z = as.vector(scale(df$Scrounging_events)))

stan_data_full_prd =  list(n_obs = nrow(df),
                                 n_ind = length(unique(df$RingNR)),
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
                                 z = as.vector(scale(df$Producing_events)))


write(
  temp <-   "data {
//Phenotype dataset
 
  // Number of clusters (an integer)
   int<lower=0> n_obs; // number of observations for phenotypes or number of rows
   int<lower=0> n_ind; // number of individuals

   int<lower=0> n_trials; // number of trials
   int<lower=0> n_groups; // number of groups

  // Clusters identifiers (an integer)
   int<lower=0> individual[n_obs];  //  Individual ID repeated obs
   int<lower=0> opponent1[n_obs];  //  Individual ID opponent1 repeated obs
   int<lower=0> opponent2[n_obs];  //  Individual ID opponent2 repeated obs
   
   int<lower=0> trialID[n_obs];  //  Individual ID opponent repeated obs
   int<lower=0> groupID[n_obs];  //  Individual ID opponent repeated obs

  //Predictors
  real year[n_obs]; // fixed effect of year 1/2 (2022/2023)
  real trialday[n_obs]; // fixed effect of trial order/sequence per individual

  real trial_order[n_obs]; // fixed effect of trialday 1/2
  int<lower=0,upper=1> group_trial_b4[n_obs]; // fixed effect of trialday 1/2

  // Continuous variables
 	real z[n_obs];  // phenotypic observations

 }
 
 parameters {
// Define parameters to estimate
  // Parameters for model of phenotype
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
   matrix[2,n_ind]     	zI; //(intercepts and opponent effect for each individual)
   vector<lower=0>[2]  	sigma_I; // sd  intercepts and slopes
   cholesky_factor_corr[2] L;  // factor to estimate covariance int-slopes
   
   vector[n_trials] zTrial_I; //
   vector[n_groups] zGroup_I; //

   real <lower=0> sigma_e;
   
   real <lower=0> sigma_trial;
   real <lower=0> sigma_group;
 }
 
 transformed parameters{
	matrix[2,n_ind] I; //  Unscaled blups intercept and res_impact for each individual
	real e_z[n_obs]; // predicted values for phenotype
	
  vector[n_trials] trial_I; //
  vector[n_groups] group_I; //

  I  = diag_pre_multiply(sigma_I, L) * zI; // get the unscaled value
  
  trial_I  =  zTrial_I * sigma_trial; // get the unscaled value
  group_I  =  zGroup_I * sigma_group; // get the unscaled value

   for (i in 1:n_obs) {
 	e_z[i]  = B_0  + B_year*year[i] + B_trialday*trialday[i] + B_trial_order*trial_order[i] + B_group_trial_b4*group_trial_b4[i] + int_trialday_trialorder*trialday[i]*trial_order[i] + int_trialday_group_trial_b4*trialday[i]*group_trial_b4[i] + int_group_trial_b4_trialorder*group_trial_b4[i]*trial_order[i] + int_trialday_group_trial_b4_trialorder*trialday[i]*group_trial_b4[i]*trial_order[i] + I[1, individual[i]] + I[2, opponent1[i]] + I[2, opponent2[i]] + trial_I[trialID[i]] + group_I[groupID[i]];
   }
   
}
 
model {
// Create vector of predicted values
 B_0 ~ normal(0, 1);  

 B_year ~ normal(0, 1);
 B_trial_order ~ normal(0, 1);

 B_trialday ~ normal(0, 1);
 B_group_trial_b4 ~ normal(0, 1);

 int_trialday_group_trial_b4_trialorder ~ normal(0, 1);

 int_trialday_trialorder ~ normal(0, 1); 
 int_trialday_group_trial_b4 ~ normal(0, 1);
 int_group_trial_b4_trialorder ~ normal(0, 1);


 // Random effects distribution
	to_vector(zI) ~ normal(0,1);

	to_vector(zTrial_I) ~ normal(0,1);
	to_vector(zGroup_I) ~ normal(0,1);

 to_vector(sigma_I) ~ normal(0,1);
 
 sigma_trial ~ exponential(3);
 sigma_group ~ exponential(3);
 sigma_e ~exponential(3);

 L ~ lkj_corr_cholesky(1);
    
    
 // Likelihood function
	for (i in 1:n_obs)
 	z[i]~normal(e_z[i], sigma_e);
}

generated quantities{
real<lower=0> Sigma2_intercept;
real<lower=0> Sigma2_res_impact;

real<lower=0> Sigma2_trialID; // variance trialID
real<lower=0> Sigma2_groupID;

real<lower=0> Sigma2_res;
real<lower=0> Sigma2_total;

real cov_1; //   int_res_impact
real cor_1;//    int_res_impact
matrix[2, 2]  Omega_I;

Sigma2_intercept = sigma_I[1]^2;
Sigma2_res_impact= sigma_I[2]^2;

Sigma2_groupID = sigma_group^2;
Sigma2_trialID = sigma_trial^2;

Sigma2_res = sigma_e^2;

Sigma2_total = Sigma2_intercept + Sigma2_res_impact + Sigma2_groupID + Sigma2_trialID + Sigma2_res ;

Omega_I = L * L';
cov_1 = Omega_I[1,2]*sqrt(Sigma2_res_impact*Sigma2_intercept);

cor_1 = cov_1 /sqrt(Sigma2_intercept*Sigma2_res_impact);

// Estimate variance components of random effects

real var_comp_focal;
real var_comp_opponent;
real var_comp_group;
real var_comp_trial;
real var_comp_res;

var_comp_focal = Sigma2_intercept/Sigma2_total;
var_comp_opponent = Sigma2_res_impact/Sigma2_total ;
var_comp_group = Sigma2_groupID/Sigma2_total;
var_comp_trial = Sigma2_trialID/Sigma2_total;
var_comp_res = Sigma2_res/Sigma2_total;

}"

, file="Models/univ_full.stan")

#  Scrounging 

md_full_scr <- stan("Models/univ_full.stan", data = stan_data_full_scr, 
                     pars = params_full, chains = 5, iter = 5000,  
                     warmup = 2000, thin = 1, cores = parallel::detectCores()-1,
                     seed = 681484520
                    )
beep(0)

rstan::get_seed(md_full_scr)#681484520

# Mean + CI + median + ESS + Rhat
round(summary(md_full_scr)$summary[,c(1,4,6,8,9,10)],3)
Sum_md_scr_full <- round(summary(md_full_scr)$summary[,c(1,4,6,8,9,10)],3)

#launch_shinystan(md_full_scr)
saveRDS(md_full_scr, file = paste("Model_output/Split_models/md_scr_full.rds", sep = ""))
saveRDS(Sum_md_scr_full, file = paste("Model_output/Split_models/Sum_md_scr_full.rds", sep = ""))

#===============================================================================
# Producing

md_full_prd <- stan("Models/univ_full.stan", data = stan_data_full_prd, 
                     pars = params_full, chains = 5, iter = 5000,  
                     warmup = 2000, thin = 1, cores = parallel::detectCores()-1,
                     seed = 861086328
                      )
beep(0)

rstan::get_seed(md_full_prd2)#861086328


# Mean + CI + median + ESS + Rhat
round(summary(md_full_prd)$summary[,c(1,4,6,8,9,10)],3)
Sum_md_prd_full <- round(summary(md_full_prd)$summary[,c(1,4,6,8,9,10)],3)

#launch_shinystan(md_full_prd2)
saveRDS(md_full_prd, file = paste("Model_output/Split_models/md_prd_full.rds", sep = ""))
saveRDS(Sum_md_prd_full, file = paste("Model_output/Split_models/Sum_md_prd_full.rds", sep = ""))


#===============================================================================
### End of script ###
