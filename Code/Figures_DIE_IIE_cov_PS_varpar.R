#===============================================================================

# Figures MS variance partitioning models PS DIE-IIE covariance
# Producing & Scrounging sparrow project

# Author: Corné de Groot

#===============================================================================

# Save local package information:
#install.packages("renv")
#renv::init()

# Load local package information:
renv::restore()

#load in packages
library(ggplot2)
library(ggExtra)
library(stringr)
library(sjPlot)
library(dplyr)
library(tidyr)
library(rstan)
library(introdataviz)

# Package versions for reproducibility
packageVersion("ggplot2") #3.5.0
packageVersion("ggExtra") #0.10.1
packageVersion("stringr") #1.5.0
packageVersion("sjPlot") #2.8.15
packageVersion("dplyr") #1.1.3
packageVersion("tidyr") #1.3.0
packageVersion("rstan") #2.32.6
packageVersion("introdataviz") #0.0.0.9003

R.Version() #"R version 4.3.1 

#===============================================================================

# Set wd if loading the df does not work
#setwd(rstudioapi::getActiveProject())

# Load in the raw data for means and standard
df <- read.csv("Data/Data_PS_DIE_IIE_cov.csv", header = TRUE)

Res <- 1000 # Resolution used in the plots that are saved as png.

# Legend of figures: (All work independently)

# Main figures
# Fig 2AB: Cov/cor matrices
# Fig 2C: Derive psi from posterior distribution
# Fig 2C: Total phenotpyic efect table 
# Fig 2E: Stacked barplot

# Fig3: Cross-year repeatability plot
# Fig4: BLUP figures to show individual variation in traits and among individual correlations

# Supplemental figures & tables
# Suppl fig 5: Multinormal residual distribution
# Suppl tab2: Bivariate model tables for full model and adjusted repeatability model 
# Suppl tab3: Model tables for short & cross-year repeatability  
# Suppl tab4 & 5: Table for all split models

#===============================================================================

# Fig 2AB: Cov/cor matrices

# DGE-IGE correlations PS
Biv_PS_mod <- as.data.frame(readRDS("Model_output/Sum_md_PS_bivar_full.rds"))
Biv_PS_mod_unadjusted <- as.data.frame(readRDS("Model_output/Sum_md_PS_bivar_full_unadjusted_rep.rds"))
df <- read.csv("Data/Data_PS_DIE_IIE_cov.csv", header = TRUE)

Biv_PS_mod$Variable <- rownames(Biv_PS_mod)
Covs_1 <- Biv_PS_mod %>% filter(str_detect(Variable,"cov"))
Cors_1 <- Biv_PS_mod %>% filter(str_detect(Variable,"cor"))
Vars <- Biv_PS_mod %>% filter(str_detect(Variable,"var"))
Reps <- Biv_PS_mod %>% filter(str_detect(Variable,"rep"))

Vars_1 <- Vars %>% filter(Variable %in% c("var_DIE", "var_IIE", "var_DIE_2", "var_IIE_2"))
Reps_1 <- Reps %>% filter(Variable %in% c("rep_DIE", "rep_IIE", "rep_DIE_2", "rep_IIE_2"))

Res_var <- Biv_PS_mod %>% filter(Variable %in% c("var_res", "var_res_2", "res_cov"))
Res_rep <- Biv_PS_mod %>% filter(Variable %in% c("rep_res", "rep_res_2", "res_cor"))

# Back-transform where necessary
Covs_1 <- Covs_1 %>% mutate(across(c("2.5%", "50%", "97.5%"), 
                ~ case_when(Variable == "cov_P1" ~ .* sd(df$Producing_events), TRUE ~ .)))
Covs_1 <- Covs_1 %>% mutate(across(c("2.5%", "50%", "97.5%"), 
                     ~ case_when(Variable %in% c("cov_P2", "cov_P3", "cov_P4", "cov_P5") ~ 
                     .* (sd(df$Producing_events)*sd(df$Scrounging_events)), TRUE ~ .)))
Covs_1 <- Covs_1 %>% mutate(across(c("2.5%", "50%", "97.5%"), 
                      ~ case_when(Variable == "cov_P6" ~ .* sd(df$Scrounging_events), TRUE ~ .)))

Vars_1 <- Vars_1 %>% mutate(across(c("2.5%", "50%", "97.5%"), 
                ~ ifelse(grepl("_2", Variable), 
                .* var(df$Scrounging_events), 
                .* var(df$Producing_events))))

Res_var <- Res_var %>% mutate(across(c("2.5%", "50%", "97.5%"), 
                      ~ case_when(Variable == "var_res" ~ .* sd(df$Producing_events), TRUE ~ .)))
Res_var <- Res_var %>% mutate(across(c("2.5%", "50%", "97.5%"), 
                      ~ case_when(Variable == "var_res_2" ~ .* sd(df$Scrounging_events), TRUE ~ .)))
Res_var <- Res_var %>% mutate(across(c("2.5%", "50%", "97.5%"), 
                      ~ case_when(Variable == "res_cov" ~ .* (sd(df$Producing_events)*sd(df$Scrounging_events)), TRUE ~ .)))

# Derive the estimate + 95% credible interval 
Covs_1 <- Covs_1 %>% mutate(across(c("2.5%", "50%", "97.5%"), ~ format(format(round(., 3), nsmall = 3), nsmall = 3)))
Covs_1 <- Covs_1 %>% unite(`2.5%`, `97.5%`, col = "Credible interval", sep = " ; ", remove = F)
Covs_1 <- Covs_1 %>% mutate(`Credible interval`= paste0("(",`Credible interval`, ")"))
Covs_1 <- Covs_1 %>% unite(`50%`, `Credible interval`, col = "Estimate", sep = " ", remove = F)
Covs_1 <- Covs_1 %>% relocate("Estimate")

Cors_1 <- Cors_1 %>% mutate(across(c("2.5%", "50%", "97.5%"), ~ format(round(., 3), nsmall = 3)))
Cors_1 <- Cors_1 %>% unite(`2.5%`, `97.5%`, col = "Credible interval", sep = " ; ", remove = F)
Cors_1 <- Cors_1 %>% mutate(`Credible interval`= paste0("(",`Credible interval`, ")"))
Cors_1 <- Cors_1 %>% unite(`50%`, `Credible interval`, col = "Estimate", sep = " ", remove = F)
Cors_1 <- Cors_1 %>% relocate("Estimate")

Vars_1 <- Vars_1 %>% mutate(across(c("2.5%", "50%", "97.5%"), ~ format(round(., 3), nsmall = 3)))
Vars_1 <- Vars_1 %>% unite(`2.5%`, `97.5%`, col = "Credible interval", sep = " ; ", remove = F)
Vars_1 <- Vars_1 %>% mutate(`Credible interval`= paste0("(",`Credible interval`, ")"))
Vars_1 <- Vars_1 %>% unite(`50%`, `Credible interval`, col = "Estimate", sep = " ", remove = F)
Vars_1 <- Vars_1 %>% relocate("Estimate")

Reps_1 <- Reps_1 %>% mutate(across(c("2.5%", "50%", "97.5%"), ~ format(round(., 3), nsmall = 3)))
Reps_1 <- Reps_1 %>% unite(`2.5%`, `97.5%`, col = "Credible interval", sep = " ; ", remove = F)
Reps_1 <- Reps_1 %>% mutate(`Credible interval`= paste0("(",`Credible interval`, ")"))
Reps_1 <- Reps_1 %>% unite(`50%`, `Credible interval`, col = "Estimate", sep = " ", remove = F)
Reps_1 <- Reps_1 %>% relocate("Estimate")

Res_var <- Res_var %>% mutate(across(c("2.5%", "50%", "97.5%"), ~ format(round(., 3), nsmall = 3)))
Res_var <- Res_var %>% unite(`2.5%`, `97.5%`, col = "Credible interval", sep = " ; ", remove = F)
Res_var <- Res_var %>% mutate(`Credible interval`= paste0("(",`Credible interval`, ")"))
Res_var <- Res_var %>% unite(`50%`, `Credible interval`, col = "Estimate", sep = " ", remove = F)
Res_var <- Res_var %>% relocate("Estimate")

Res_rep <- Res_rep %>% mutate(across(c("2.5%", "50%", "97.5%"), ~ format(format(round(., 3), nsmall = 3), nsmall = 3)))
Res_rep <- Res_rep %>% unite(`2.5%`, `97.5%`, col = "Credible interval", sep = " ; ", remove = F)
Res_rep <- Res_rep %>% mutate(`Credible interval`= paste0("(",`Credible interval`, ")"))
Res_rep <- Res_rep %>% unite(`50%`, `Credible interval`, col = "Estimate", sep = " ", remove = F)
Res_rep <- Res_rep %>% relocate("Estimate")

# Covariance matrix
Tab1 <- as.data.frame(matrix(nrow = 4, ncol = 4))
row.names(Tab1) <- c("Producing DIE", "Producing IIE", "Scrounging DIE", "Scrounging IIE")
colnames(Tab1) <- c("Producing DIE", "Producing IIE", "Scrounging DIE", "Scrounging IIE")
Tab2 <- Tab1
Tab1_sym <- Tab1

# variances in table
Tab1[1,1] <- print(Vars_1["var_DIE","Estimate"])
Tab1[2,2] <- print(Vars_1["var_IIE","Estimate"])
Tab1[3,3] <- print(Vars_1["var_DIE_2","Estimate"])
Tab1[4,4] <- print(Vars_1["var_IIE_2","Estimate"])

# Covariances in table
Tab1[1,2] <- print(Covs_1["cov_P1","Estimate"])
Tab1[1,3] <- print(Covs_1["cov_P2","Estimate"])
Tab1[1,4] <- print(Covs_1["cov_P3","Estimate"])
Tab1[2,3] <- print(Covs_1["cov_P4","Estimate"])
Tab1[2,4] <- print(Covs_1["cov_P5","Estimate"])
Tab1[3,4] <- print(Covs_1["cov_P6","Estimate"])

Tab1[is.na(Tab1)] <- "-"

Tab1$. <- rownames(Tab1)
Tab1 <- Tab1 %>% relocate(".")

tab_df(Tab1, alternate.rows = T)
sjPlot::tab_df(Tab1, alternate.rows = T, file = paste("Plots/Table_cov_bivar_PS.doc"))

# Correlation matrix
# Repeatabilities in table
Tab2[1,1] <- print(Reps_1["rep_DIE","Estimate"])
Tab2[2,2] <- print(Reps_1["rep_IIE","Estimate"])
Tab2[3,3] <- print(Reps_1["rep_DIE_2","Estimate"])
Tab2[4,4] <- print(Reps_1["rep_IIE_2","Estimate"])

# Correlations in table
Tab2[2,1] <- print(Cors_1["cor_P1","Estimate"])
Tab2[3,1] <- print(Cors_1["cor_P2","Estimate"])
Tab2[4,1] <- print(Cors_1["cor_P3","Estimate"])
Tab2[3,2] <- print(Cors_1["cor_P4","Estimate"])
Tab2[4,2] <- print(Cors_1["cor_P5","Estimate"])
Tab2[4,3] <- print(Cors_1["cor_P6","Estimate"])

Tab2[is.na(Tab2)] <- "-"

Tab2$. <- rownames(Tab2)
Tab2 <- Tab2 %>% relocate(".")

tab_df(Tab2, alternate.rows = T)
sjPlot::tab_df(Tab2, alternate.rows = T, file = paste("Plots/Table_cor_bivar_PS.doc"))

# residual covarinace matrix
Tab3 <- as.data.frame(matrix(nrow = 2, ncol = 2))
row.names(Tab3) <- c("Residual producing", "Residual scrounging")
colnames(Tab3) <- c("Residual producing", "Residual scrounging")
Tab4 <- Tab3

# Variances in table
Tab3[1,1] <- print(Res_var["var_res","Estimate"])
Tab3[2,2] <- print(Res_var["var_res_2","Estimate"])

# Covariance in table
Tab3[1,2] <- print(Res_var["res_cov","Estimate"])

Tab3[is.na(Tab3)] <- "-"

Tab3$. <- rownames(Tab3)
Tab3 <- Tab3 %>% relocate(".")

tab_df(Tab3, alternate.rows = T)
sjPlot::tab_df(Tab3, alternate.rows = T, file = paste("Plots/Table_rescov_bivar_PS.doc"))

# Residual correlation matrix
# Repeatabilities in table
Tab4[1,1] <- print(Res_rep["rep_res","Estimate"])
Tab4[2,2] <- print(Res_rep["rep_res_2","Estimate"])

# Correlation in table
Tab4[2,1] <- print(Res_rep["res_cor","Estimate"])

Tab4[is.na(Tab4)] <- "-"

Tab4$. <- rownames(Tab4)
Tab4 <- Tab4 %>% relocate(".")

tab_df(Tab4, alternate.rows = T)
sjPlot::tab_df(Tab4, alternate.rows = T, file = paste("Plots/Table_rescor_bivar_PS.doc"))

#===============================================================================
# Fig 2C: Derive psi from posterior distribution

md_bivar <- readRDS("Model_output/md_PS_bivar_full.rds")
df <- read.csv("Data/Data_PS_DIE_IIE_cov.csv", header = TRUE)

posterior_bivar_covs <- as.data.frame(rstan::extract(md_bivar, pars = c("cov_P1", "cov_P2", "cov_P3", "cov_P4", "cov_P5", "cov_P6",
                                                                        "var_DIE", "var_DIE_2", "var_IIE", "var_IIE_2")))

posterior_bivar_covs$iteration <- seq_along(1:nrow(posterior_bivar_covs))

round(summary(md_bivar, pars = c("cov_P1", "cov_P2", "cov_P3", "cov_P4", "cov_P5", "cov_P6",
                                         "var_DIE", "var_DIE_2", "var_IIE", "var_IIE_2"))$summary[,c(1,4,6,8,9,10)],3)

# Backtransforming is only necessary if you want the psi's on the "original" unstandardised scale

# # Backtransform variances
# posterior_bivar_covs$var_DIE <- posterior_bivar_covs$var_DIE * var(df$Producing_events)
# posterior_bivar_covs$var_DIE_2 <- posterior_bivar_covs$var_DIE_2 * var(df$Scrounging_events)  
# posterior_bivar_covs$var_IIE <- posterior_bivar_covs$var_IIE * var(df$Producing_events)
# posterior_bivar_covs$var_IIE_2 <- posterior_bivar_covs$var_IIE_2 * var(df$Scrounging_events)  
# 
# # Backtransform covariances
# posterior_bivar_covs$cov_P1 <- posterior_bivar_covs$cov_P1 * sd(df$Producing_events)
# posterior_bivar_covs$cov_P2 <- posterior_bivar_covs$cov_P2 * (sd(df$Producing_events) * sd(df$Scrounging_events))   
# posterior_bivar_covs$cov_P3 <- posterior_bivar_covs$cov_P3 * (sd(df$Producing_events) * sd(df$Scrounging_events))
# posterior_bivar_covs$cov_P4 <- posterior_bivar_covs$cov_P4 * (sd(df$Producing_events) * sd(df$Scrounging_events))    
# posterior_bivar_covs$cov_P5 <- posterior_bivar_covs$cov_P5 * (sd(df$Producing_events) * sd(df$Scrounging_events))
# posterior_bivar_covs$cov_P6 <- posterior_bivar_covs$cov_P6 * sd(df$Scrounging_events)    

# Derive psi matrix per posterior sample (iteration)
DIE_mat <- list()

for(i in 1:nrow(posterior_bivar_covs)){
  DIE_mat[[i]] <- matrix(c(posterior_bivar_covs[i, "var_DIE"], 
                           posterior_bivar_covs[i, "cov_P2"], 
                           posterior_bivar_covs[i, "cov_P2"], 
                           posterior_bivar_covs[i, "var_DIE_2"]), 
                         nrow = 2, ncol =2, byrow = TRUE)
} 

IIE_DIE_cov_mat <- list()

for(i in 1:nrow(posterior_bivar_covs)){
  IIE_DIE_cov_mat[[i]] <- matrix(c(posterior_bivar_covs[i, "cov_P1"], 
                                   posterior_bivar_covs[i, "cov_P4"], 
                                   posterior_bivar_covs[i, "cov_P3"], 
                                   posterior_bivar_covs[i, "cov_P6"]), 
                                 nrow = 2, ncol =2, byrow = TRUE)
} 


# Adjusted for group size (equation 26 Brodie & McGlothlin 2009)
n <- 3
psi_mat_posterior_adj <- array(NA, dim = c(nrow(posterior_bivar_covs), 2, 2))

for(i in 1:nrow(posterior_bivar_covs)){
  psi_mat_posterior_adj[i,,] <- IIE_DIE_cov_mat[[i]] %*% solve(DIE_mat[[i]] + (n - 2) * IIE_DIE_cov_mat[[i]])
} 

# Summarize posterior results
psi_mat_adj <- as.data.frame(apply(psi_mat_posterior_adj, c(2,3), quantile,probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
psi_mat_adj <- as.data.frame(t(psi_mat_adj))
psi_mat_adj <- psi_mat_adj %>% mutate(across(c("2.5%", "50%", "97.5%"), ~ format(round(., 3), nsmall = 3)))
psi_mat_adj <- psi_mat_adj %>% unite(`2.5%`, `97.5%`, col = "Credible interval", sep = " ; ", remove = F)
psi_mat_adj <- psi_mat_adj %>% mutate(`Credible interval`= paste0("(",`Credible interval`, ")"))
psi_mat_adj <- psi_mat_adj %>% unite(`50%`, `Credible interval`, col = "Estimate", sep = " ", remove = F)
psi_mat_adj <- psi_mat_adj %>% relocate("Estimate")

# Make (capital letter) psi table
Tab_psi_adj <- as.data.frame(matrix(nrow = 2, ncol = 2))
row.names(Tab_psi_adj) <- c("Producing", "Scrounging")
colnames(Tab_psi_adj) <- c("Producing", "Scrounging")

# Psi in table
Tab_psi_adj[1,1] <- print(psi_mat_adj["1.1","Estimate"])
Tab_psi_adj[1,2] <- print(psi_mat_adj["1.2","Estimate"])
Tab_psi_adj[2,1] <- print(psi_mat_adj["2.1","Estimate"])
Tab_psi_adj[2,2] <- print(psi_mat_adj["2.2","Estimate"])

Tab_psi_adj[is.na(Tab_psi_adj)] <- "-"

Tab_psi_adj$. <- rownames(Tab_psi_adj)
Tab_psi_adj <- Tab_psi_adj %>% relocate(".")

tab_df(Tab_psi_adj, alternate.rows = T)
sjPlot::tab_df(Tab_psi_adj, alternate.rows = T, file = paste("Plots/Table_psi_adj_PS.doc"))



# Without accounting for groupsize, not reported in MS
# Multivariate psi (equation 14 Brodie & McGlothlin 2009)

psi_mat_posterior <- array(NA, dim = c(nrow(posterior_bivar_covs), 2, 2))

for(i in 1:nrow(posterior_bivar_covs)){
  psi_mat_posterior[i,,] <-  IIE_DIE_cov_mat[[i]] %*% solve(DIE_mat[[i]]) 
} 

# Summarize posterior results
psi_mat <- as.data.frame(apply(psi_mat_posterior, c(2,3), quantile,probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
psi_mat <- as.data.frame(t(psi_mat))
psi_mat <- psi_mat %>% mutate(across(c("2.5%", "50%", "97.5%"), ~ format(round(., 3), nsmall = 3)))
psi_mat <- psi_mat %>% unite(`2.5%`, `97.5%`, col = "Credible interval", sep = " ; ", remove = F)
psi_mat <- psi_mat %>% mutate(`Credible interval`= paste0("(",`Credible interval`, ")"))
psi_mat <- psi_mat %>% unite(`50%`, `Credible interval`, col = "Estimate", sep = " ", remove = F)
psi_mat <- psi_mat %>% relocate("Estimate")

# Make (capital letter) psi table
Tab_psi <- as.data.frame(matrix(nrow = 2, ncol = 2))
row.names(Tab_psi) <- c("Producing", "Scrounging")
colnames(Tab_psi) <- c("Producing", "Scrounging")

# Psi in table
Tab_psi[1,1] <- print(psi_mat["1.1","Estimate"])
Tab_psi[1,2] <- print(psi_mat["1.2","Estimate"])
Tab_psi[2,1] <- print(psi_mat["2.1","Estimate"])
Tab_psi[2,2] <- print(psi_mat["2.2","Estimate"])

Tab_psi[is.na(Tab_psi)] <- "-"

Tab_psi$. <- rownames(Tab_psi)
Tab_psi <- Tab_psi %>% relocate(".")

tab_df(Tab_psi, alternate.rows = T)
sjPlot::tab_df(Tab_psi, alternate.rows = T, file = paste("Plots/Table_psi_PS.doc"))

#===============================================================================
# Fig 2D: Total phenotpyic efect table 

# Read in BLUPS from model & wrangle
Bivar_full_P_eff <- as.data.frame(readRDS("Model_output/Sum_md_PS_bivar_full.rds"))
df <- read.csv("Data/Data_PS_DIE_IIE_cov.csv", header = TRUE)

# Wrangle data
Bivar_full_P_eff$Variable <- row.names(Bivar_full_P_eff)
Bivar_full_P_eff <- Bivar_full_P_eff %>% filter(stringr::str_detect(Variable, "P_eff"))
Bivar_full_P_eff <- Bivar_full_P_eff %>% mutate(Behaviour = ifelse(str_detect(Variable, "_2"), "Scrounging", "Producing")) 

# Back-transform variables
Bivar_full_P_eff <- Bivar_full_P_eff %>%
  mutate(across(c("2.5%", "50%", "97.5%"), 
                ~ ifelse(grepl("total_P_eff", Variable) & Behaviour == "Scrounging", . * var(df$Scrounging_events), 
                         ifelse(grepl("total_P_eff", Variable) & Behaviour == "Producing", . * var(df$Producing_events), .))))

Bivar_full_P_eff <- Bivar_full_P_eff %>% mutate(across(c("2.5%", "50%", "97.5%"), ~ format(round(., 3), nsmall = 3)))
Bivar_full_P_eff <- Bivar_full_P_eff %>% unite(`2.5%`, `97.5%`, col = "Credible interval", sep = " ; ", remove = F)
Bivar_full_P_eff <- Bivar_full_P_eff %>% mutate(`Credible interval`= paste0("(",`Credible interval`, ")"))
Bivar_full_P_eff <- Bivar_full_P_eff %>% unite(`50%`, `Credible interval`, col = "Estimate", sep = " ", remove = F)
Bivar_full_P_eff <- Bivar_full_P_eff %>% relocate("Estimate")

# Make (capital letter) psi table
Tab_P_eff <- as.data.frame(matrix(nrow = 2, ncol = 2))
row.names(Tab_P_eff) <- c("Variance", "R")
colnames(Tab_P_eff) <- c("Producing", "Scrounging")

# Psi in table
Tab_P_eff[1,1] <- print(Bivar_full_P_eff["total_P_eff_1","Estimate"])
Tab_P_eff[1,2] <- print(Bivar_full_P_eff["total_P_eff_2","Estimate"])
Tab_P_eff[2,1] <- print(Bivar_full_P_eff["rep_P_eff_1","Estimate"])
Tab_P_eff[2,2] <- print(Bivar_full_P_eff["rep_P_eff_2","Estimate"])

Tab_P_eff[is.na(Tab_P_eff)] <- "-"

Tab_P_eff$. <- rownames(Tab_P_eff)
Tab_P_eff <- Tab_P_eff %>% relocate(".")

tab_df(Tab_P_eff, alternate.rows = T)
sjPlot::tab_df(Tab_P_eff, alternate.rows = T, file = paste("Plots/Table_total_phenotypic_effect_PS.doc"))

#===============================================================================
# Fig 2E: Stacked barplot

Barchart_df <- Vars[c(1:5,7:11), ]# without the total variance
Barchart_df$Behaviour <- rep(c("Producing", "Scrounging"), each = 5)
Barchart_df <- Barchart_df %>% rename(Median = "50%")

# Remove redundant parts of variable names
Barchart_df$Variable <- gsub("_2","",as.character(Barchart_df$Variable))
Barchart_df$Variable <- gsub("var_","",as.character(Barchart_df$Variable))
# Rename variable names
Barchart_df$Variable <- gsub("res","Residual",as.character(Barchart_df$Variable))
Barchart_df$Variable <- gsub("trialID","Trial",as.character(Barchart_df$Variable))
Barchart_df$Variable <- gsub("groupID","Group",as.character(Barchart_df$Variable))

# Duplicate IIE variances for social partner 1 & 2
Barchart_df <- rbind(Barchart_df, Barchart_df[rep(2, 1), ])# Producing IIE
Barchart_df <- rbind(Barchart_df, Barchart_df[rep(7, 1), ])# Scrounging IIE
# Rename IIE variables
Barchart_df[2,7] <- "IIE partner1"
Barchart_df[7,7] <- "IIE partner1"
Barchart_df[11,7] <- "IIE partner2"
Barchart_df[12,7] <- "IIE partner2"

# Set order of variables in the plot
Barchart_df$Variable <- factor(Barchart_df$Variable, levels = c("Residual", "Trial", "Group", "IIE partner1", "IIE partner2", "DIE"))
level_order <- c('Producing', 'Scrounging') 

# Variance component
Stacked_bar_plot_PS <- ggplot(Barchart_df, aes(y = Median, x = factor(Behaviour, level = level_order), fill = Variable)) +
  geom_bar(position = "fill", stat = "identity", color = "black") + 
  theme_classic() + ylab("Variance component")+ xlab("") + 
  scale_fill_brewer(palette = "Greys", direction = -1) +
  theme_classic(base_size = 17) +
  theme(legend.title=element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"))

Stacked_bar_plot_PS

# Save plot(s)
png(file="Plots/Stacked_barplot_VC.png", width = 5.5*Res, height = 7*Res, res = Res)
Stacked_bar_plot_PS
dev.off()

#===============================================================================
# Fig3: Cross-year repeatability plot

Short_rep_prd <- as.data.frame(readRDS("Model_output/Sum_md_short_rep_prd.rds"))
Short_rep_scr <- as.data.frame(readRDS("Model_output/Sum_md_short_rep_scr.rds"))

Short_rep_prd$Variable <- paste0(row.names(Short_rep_prd),"_Prd")
Short_rep_scr$Variable <- paste0(row.names(Short_rep_scr),"_Scr")

Var1 <- rbind(Short_rep_prd,Short_rep_scr)

Rep_plot_data <- Var1 %>% filter(Variable %in% c("Intercept_rep_Prd", "Intercept_rep_opp_Prd", "Intercept_rep_Scr", "Intercept_rep_opp_Scr"))
Rep_plot_data <- Rep_plot_data %>% rename(Median = "50%",
                                          CI1 = "2.5%",
                                          CI2 = "97.5%")
Rep_plot_data <- Rep_plot_data %>% mutate(Behaviour = c("Producing","Prodicing","Scrounging", "Scrounging"),
                                          Variable = c("Focal", "Opponent", "Focal", "Opponent"))

# Half violin plot

Short_rep_prd_md <- readRDS("Model_output/md_short_rep_prd.rds")
Short_rep_scr_md <- readRDS("Model_output/md_short_rep_scr.rds")

posterior_prd <- as.data.frame(rstan::extract(Short_rep_prd_md, pars = c("Intercept_rep", "Intercept_rep_opp")))
posterior_prd <- posterior_prd %>%
  pivot_longer(cols = starts_with("Intercept"), names_to = "Variable", values_to = "Posterior")
posterior_prd$Behaviour <- "Producing"

posterior_scr <- as.data.frame(rstan::extract(Short_rep_scr_md, pars = c("Intercept_rep", "Intercept_rep_opp")))
posterior_scr <- posterior_scr %>%
  pivot_longer(cols = starts_with("Intercept"), names_to = "Variable", values_to = "Posterior")
posterior_scr$Behaviour <- "Scrounging"

dat_long <- rbind(posterior_prd, posterior_scr)

# Filter for 95% credible intervals
CI_df <- dat_long %>%
  group_by(Variable, Behaviour) %>%
  summarise(lower = quantile(Posterior, probs = 0.025),
            upper = quantile(Posterior, probs = 0.975))

CI_df <- dat_long %>%
  left_join(CI_df, by = c("Variable", "Behaviour")) %>%
  filter(Posterior >= lower & Posterior <= upper)

# Split violin plot
#devtools::install_github("psyteachr/introdataviz")

Violin_plot <- ggplot(CI_df, aes(x = Variable, y = Posterior, fill = Behaviour)) +
  introdataviz::geom_split_violin(alpha = 0.5, trim = TRUE) +
  geom_boxplot(width = 0.175, alpha = 0.5,show.legend = F, outlier.shape = NA) +
  geom_hline(yintercept = 1, linetype ="dashed", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype ="dashed", linewidth = 0.5) +
  scale_x_discrete(name = "", labels = c("DIE", "IIE")) +
  scale_y_continuous(name = "Cross-year repeatability",
                     breaks = seq(0, 1, 0.25), 
                     limits = c(0, 1.025)) +
  scale_fill_brewer(palette = "Greys", name = "") +
  theme_classic(base_size = 17) + 
  theme(legend.position = "bottom", legend.margin=margin(-20,0,0,0),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"))

Violin_plot

# Save plot
png(file="Plots/Cross_year_repeatability_violin_PS.png", width = Res*6, height = Res*5, res = Res)
Violin_plot
dev.off()

#===============================================================================
# Fig4: BLUP figures to show individual variation in traits and among individual correlations

# Read in BLUPS from model & wrangle
BLUPS <- as.data.frame(readRDS("Model_output/BLUPS_I_bivar_PS.rds"))
Corrs <- as.data.frame(readRDS("Model_output/Sum_md_PS_bivar_full.rds"))
df <- read.csv("Data/Data_PS_DIE_IIE_cov.csv", header = TRUE)

BLUPS$RowName <- row.names(BLUPS)
# Number of unique animals in the dataset
Na <- length(unique(df$RingNR)) 
BLUPS$ID <- rep(1:Na, each = 4)
BLUPS$Variable <- rep(paste(c("Producing_DIE", "Producing_IIE", "Scrounging_DIE", "Scrounging_IIE")), Na)
# Back-transform
BLUPS <- BLUPS %>% mutate(across(c("2.5%", "50%", "97.5%"), 
                                 ~ ifelse(grepl("Scrounging", .), 
                                          .* sd(df$Scrounging_events), 
                                          .* sd(df$Producing_events))))
BLUPS <- BLUPS %>% select("2.5%", "50%", "97.5%", "Variable", "ID")
BLUPS <- BLUPS %>% rename(CI1 = "2.5%", Median = "50%", CI2 = "97.5%")
BLUPS <- BLUPS %>% pivot_wider(names_from = Variable, values_from = c("CI1", "Median", "CI2"), id_cols = ID)

Corrs$Variable <- row.names(Corrs)
Corrs <- Corrs %>% filter(stringr::str_detect(Variable,"cor"))
Corrs <- Corrs %>% mutate(across(c("2.5%", "50%", "97.5%"), ~ format(round(., 3), nsmall = 3)))
Corrs <- Corrs %>% unite(`2.5%`, `97.5%`, col = "Credible interval", sep = " ; ", remove = F)
Corrs <- Corrs %>% mutate(`Credible interval`= paste0("(",`Credible interval`, ")"))
Corrs <- Corrs %>% unite(`50%`, `Credible interval`, col = "Estimate", sep = " ", remove = F)

# Set axes limits based on min and max values
Lim_min_P_D <- round(min(BLUPS$CI1_Producing_DIE),1) - 0.05
Lim_max_P_D <- round(max(BLUPS$CI2_Producing_DIE),1) + 0.05
Lim_min_P_I <- round(min(BLUPS$CI1_Producing_IIE),1) - 0.05
Lim_max_P_I <- round(max(BLUPS$CI2_Producing_IIE),1) + 0.05
Lim_min_S_D <- round(min(BLUPS$CI1_Scrounging_DIE),1) - 0.05
Lim_max_S_D <- round(max(BLUPS$CI2_Scrounging_DIE),1) + 0.05
Lim_min_S_I <- round(min(BLUPS$CI1_Scrounging_IIE),1) - 0.05
Lim_max_S_I <- round(max(BLUPS$CI2_Scrounging_IIE),1) + 0.05

Plot_cor1 <- ggplot(BLUPS, aes(y = Median_Producing_DIE, x = Median_Producing_IIE)) +
  geom_point(size = 2, alpha = 0.5) +  
  geom_errorbar(aes(ymin = CI1_Producing_DIE,ymax = CI2_Producing_DIE), alpha = 0.25) + 
  geom_errorbarh(aes(xmin = CI1_Producing_IIE,xmax = CI2_Producing_IIE), alpha = 0.25) +
  labs(y = "Propensity to produce (DIE)", x = "Producing elicited (IIE)") +
  scale_y_continuous(limits = c(Lim_min_P_D,Lim_max_P_D)) +
  scale_x_continuous(limits = c(Lim_min_P_I,Lim_max_P_I)) +
  annotate("text", x = mean(c(Lim_min_P_I,Lim_max_P_I)), y = Lim_max_P_D, label = paste0("r = ",Corrs[1,2]), size = 6, fontface = "bold") +
  theme_classic(base_size = 17) + 
  theme(legend.position = "none", 
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) 

Plot_cor2 <- ggplot(BLUPS, aes(y = Median_Producing_DIE, x = Median_Scrounging_DIE)) +
  geom_point(size = 2, alpha = 0.5) +  
  geom_errorbar(aes(ymin = CI1_Producing_DIE,ymax = CI2_Producing_DIE), alpha = 0.25) + 
  geom_errorbarh(aes(xmin = CI1_Scrounging_DIE,xmax = CI2_Scrounging_DIE), alpha = 0.25) +
  labs(y = "Propensity to produce (DIE)", x = "Propensity to scrounge (DIE)") +
  scale_y_continuous(limits = c(Lim_min_P_D,Lim_max_P_D)) +
  scale_x_continuous(limits = c(Lim_min_S_D,Lim_max_S_D)) +
  annotate("text", x = mean(c(Lim_min_S_D,Lim_max_S_D)), y = Lim_max_P_D, label = paste0("r = ",Corrs[2,2]), size = 6, fontface = "bold") +
  theme_classic(base_size = 17) + 
  theme(legend.position = "none", 
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) 

Plot_cor3 <- ggplot(BLUPS, aes(y = Median_Producing_DIE, x = Median_Scrounging_IIE)) +
  geom_point(size = 2, alpha = 0.5) +  
  geom_errorbar(aes(ymin = CI1_Producing_DIE,ymax = CI2_Producing_DIE), alpha = 0.25) + 
  geom_errorbarh(aes(xmin = CI1_Scrounging_IIE,xmax = CI2_Scrounging_IIE), alpha = 0.25) +
  labs(y = "Propensity to produce (DIE)", x = "Scrounging elicited (IIE)") +
  scale_y_continuous(limits = c(Lim_min_P_D,Lim_max_P_D)) +
  scale_x_continuous(limits = c(Lim_min_S_I,Lim_max_S_I)) +
  annotate("text", x = mean(c(Lim_min_S_I,Lim_max_S_I)), y = Lim_max_P_D, label = paste0("r = ",Corrs[3,2]), size = 6, fontface = "bold") +
  theme_classic(base_size = 17) + 
  theme(legend.position = "none", 
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) 

Plot_cor4 <- ggplot(BLUPS, aes(y = Median_Scrounging_DIE, x = Median_Producing_IIE)) +
  geom_point(size = 2, alpha = 0.5) +  
  geom_errorbar(aes(ymin = CI1_Scrounging_DIE,ymax = CI2_Scrounging_DIE), alpha = 0.25) + 
  geom_errorbarh(aes(xmin = CI1_Producing_IIE,xmax = CI2_Producing_IIE), alpha = 0.25) +
  labs(y = "Propensity to scrounge (DIE)", x = "Producing elicited (IIE)") +
  scale_y_continuous(limits = c(Lim_min_S_D,Lim_max_S_D)) +
  scale_x_continuous(limits = c(Lim_min_P_I,Lim_max_P_I)) +
  annotate("text", x = mean(c(Lim_min_P_I,Lim_max_P_I)), y = Lim_max_S_D, label = paste0("r = ",Corrs[4,2]), size = 6, fontface = "bold") +
  theme_classic(base_size = 17) + 
  theme(legend.position = "none", 
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) 

Plot_cor5 <- ggplot(BLUPS, aes(y = Median_Producing_IIE, x = Median_Scrounging_IIE)) +
  geom_point(size = 2, alpha = 0.5) +  
  geom_errorbar(aes(ymin = CI1_Producing_IIE,ymax = CI2_Producing_IIE), alpha = 0.25) + 
  geom_errorbarh(aes(xmin = CI1_Scrounging_IIE,xmax = CI2_Scrounging_IIE), alpha = 0.25) +
  labs(y = "Producing elicited (IIE)", x = "Scrounging elicited (IIE)") +
  scale_y_continuous(limits = c(Lim_min_P_I,Lim_max_P_I)) +
  scale_x_continuous(limits = c(Lim_min_S_I,Lim_max_S_I)) +
  annotate("text", x = mean(c(Lim_min_S_I,Lim_max_S_I)), y = Lim_max_P_I, label = paste0("r = ",Corrs[5,2]), size = 6, fontface = "bold") +
  theme_classic(base_size = 17) + 
  theme(legend.position = "none", 
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) 

Plot_cor6 <- ggplot(BLUPS, aes(y = Median_Scrounging_DIE, x = Median_Scrounging_IIE)) +
  geom_point(size = 2, alpha = 0.5) +  
  geom_errorbar(aes(ymin = CI1_Scrounging_DIE,ymax = CI2_Scrounging_DIE), alpha = 0.25) + 
  geom_errorbarh(aes(xmin = CI1_Scrounging_IIE,xmax = CI2_Scrounging_IIE), alpha = 0.25) +
  labs(y = "Propensity to scrounge (DIE)", x = "Scrounging elicited (IIE)") +
  scale_y_continuous(limits = c(Lim_min_S_D,Lim_max_S_D)) +
  scale_x_continuous(limits = c(Lim_min_S_I,Lim_max_S_I)) +
  annotate("text", x = mean(c(Lim_min_S_I,Lim_max_S_I)), y = Lim_max_S_D, label = paste0("r = ",Corrs[6,2]), size = 6, fontface = "bold") +
  theme_classic(base_size = 17) + 
  theme(legend.position = "none", 
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) 

grid_plot <- ggpubr::ggarrange(Plot_cor1, Plot_cor2, Plot_cor3, Plot_cor4, Plot_cor5, Plot_cor6,
                               nrow = 3, ncol = 2, labels = c("A", "B", "C", "D", "E", "F"),
                               font.label = list(color = "black",size = 17, face = "bold"),
                               align = "hv")

grid_plot
Plot_cor1
Plot_cor2
Plot_cor3
Plot_cor4
Plot_cor5
Plot_cor6

# Save plots with command
png(file="Plots/Grid_plot_P.png", width = 14*Res, height = 17*Res, res = Res)
grid_plot
dev.off()

png(file="Plots/Cor1_P.png", width = 5.5*Res, height = 7*Res, res = Res)
Plot_cor1
dev.off()

png(file="Plots/Cor2_P.png", width = 5.5*Res, height = 7*Res, res = Res)
Plot_cor2
dev.off()

png(file="Plots/Cor3_P.png", width = 5.5*Res, height = 7*Res, res = Res)
Plot_cor3
dev.off()

png(file="Plots/Cor4_P.png", width = 5.5*Res, height = 7*Res, res = Res)
Plot_cor4
dev.off()

png(file="Plots/Cor5_P.png", width = 5.5*Res, height = 7*Res, res = Res)
Plot_cor5
dev.off()

png(file="Plots/Cor6_P.png", width = 5.5*Res, height = 7*Res, res = Res)
Plot_cor6
dev.off()

#===============================================================================
# Supplemental tables & figures
#===============================================================================
# Suppl fig 5: Multinormal residual distribution

Bivar_full_md <- readRDS("Model_output/md_PS_bivar_full.rds")
df <- read.csv("Data/Data_PS_DIE_IIE_cov.csv", header = TRUE)

# Extract predicted model observations
Posterior_res <- as.data.frame(rstan::extract(Bivar_full_md, pars = c("z_exp", "z2_exp")))
rm(Bivar_full_md)

n_iter <- nrow(Posterior_res)

# Some data wrangling
Posterior_res <- Posterior_res %>%
  pivot_longer(col = starts_with("z"), names_to = c(".value", "ObservationID"), names_pattern = "(z2?|z)_(.*)")
Posterior_res$IterationID <- rep(1:nrow(df), each = n_iter)
Posterior_res$ObservationID <- str_remove(Posterior_res$ObservationID, pattern = "exp.")

# Insert the observed values
Posterior_res$Observed_prd <- rep(as.vector(scale(df$Producing_events)), times = n_iter)
Posterior_res$Observed_scr <- rep(as.vector(scale(df$Scrounging_events)), times = n_iter)

Posterior_res$Residual_prd <- Posterior_res$Observed_prd - Posterior_res$z 
Posterior_res$Residual_scr <- Posterior_res$Observed_scr - Posterior_res$z2 

# Randomly sample 100 iteration from the posterior
Posterior_res_random <- Posterior_res %>% filter(IterationID %in% sample(n_iter, 100, replace = F))

Mvn_residual_plot <- ggplot(Posterior_res_random, aes(x = Residual_prd, y = Residual_scr)) +
  geom_point(size = 2, alpha = 0.1) +
  geom_density_2d(color = "white") +
  labs(x = "Residual producing", y = "Residual scrounging") +
  scale_y_continuous(limits = c(-8,8)) +
  scale_x_continuous(limits = c(-8,8)) +
  theme_classic(base_size = 17)

Mvn_residual_plot <- ggMarginal(Mvn_residual_plot, type = "densigram", 
                                fill = "gray", color = "black",      
                                alpha = 0.75, linewidth = 1.05)

Mvn_residual_plot

# Save plots with command
# Needs to be saved manually in order for the histograms to be saved in the figure 

png(file="Plots/Mvn_residual_plot.png", width = 8*Res, height = 6*Res, res = Res)
Mvn_residual_plot
dev.off()

#===============================================================================
# Suppl tab2: Bivariate model tables for full model and adjusted repeatability model

# Read in BLUPS from model & wrangle
Bivar_full <- as.data.frame(readRDS("Model_output/Sum_md_PS_bivar_full.rds"))
Bivar_unadjusted <- as.data.frame(readRDS("Model_output/Sum_md_PS_bivar_full_unadjusted_rep.rds"))
df <- read.csv("Data/Data_PS_DIE_IIE_cov.csv", header = TRUE)

# Wrangle data
Bivar_full$Variable <- row.names(Bivar_full)
Bivar_unadjusted$Variable <- row.names(Bivar_unadjusted)

Bivar_full <- Bivar_full %>% filter(!stringr::str_detect(Variable, c("cor|cov")))
Bivar_unadjusted <- Bivar_unadjusted %>% filter(!stringr::str_detect(Variable, c("cor|cov")))

Bivar_full <- Bivar_full %>% mutate(Behaviour = ifelse(str_detect(Variable, "_2"), "Scrounging", "Producing")) 
Bivar_unadjusted<- Bivar_unadjusted %>% mutate(Behaviour = ifelse(str_detect(Variable, "_2"), "Scrounging", "Producing")) 

# Back-transform variables
Bivar_full <- Bivar_full %>%
  mutate(across(c("2.5%", "50%", "97.5%"), 
                ~ ifelse(grepl("B_|mu|int_", Variable) & Behaviour == "Scrounging", . * sd(df$Scrounging_events), 
                         ifelse(grepl("B_|mu|int_", Variable) & Behaviour == "Producing", . * sd(df$Producing_events), .))))

Bivar_full <- Bivar_full %>% mutate(across(c("2.5%", "50%", "97.5%"), 
                                           ~ ifelse(Variable == "mu" & Behaviour == "Producing", . + mean(df$Producing_events),
                                                    ifelse(Variable == "mu_2" & Behaviour == "Scrounging", . + mean(df$Scrounging_events), .))))

Bivar_full <- Bivar_full %>%
  mutate(across(c("2.5%", "50%", "97.5%"), 
                ~ ifelse(grepl("var_|total_P_eff", Variable) & Behaviour == "Scrounging", . * var(df$Scrounging_events), 
                         ifelse(grepl("var_|total_P_eff", Variable) & Behaviour == "Producing", . * var(df$Producing_events), .))))

Bivar_full <- Bivar_full %>% mutate(across(c("2.5%", "50%", "97.5%"), 
                                                 ~ ifelse(grepl("order", Variable), . / sd(df$TrialNR_Indiv), .)))

Bivar_unadjusted <- Bivar_unadjusted %>%
  mutate(across(c("2.5%", "50%", "97.5%"), 
                ~ ifelse(grepl("B_|mu|int_", Variable) & Behaviour == "Scrounging", . * sd(df$Scrounging_events), 
                         ifelse(grepl("B_|mu|int_", Variable) & Behaviour == "Producing", . * sd(df$Producing_events), .))))

Bivar_unadjusted <- Bivar_unadjusted %>% mutate(across(c("2.5%", "50%", "97.5%"), 
                                                       ~ ifelse(Variable == "mu" & Behaviour == "Producing", . + mean(df$Producing_events),
                                                                ifelse(Variable == "mu_2" & Behaviour == "Scrounging", . + mean(df$Scrounging_events), .))))

Bivar_unadjusted <- Bivar_unadjusted %>%
  mutate(across(c("2.5%", "50%", "97.5%"), 
                ~ ifelse(grepl("var_|total_P_eff", Variable) & Behaviour == "Scrounging", . * var(df$Scrounging_events), 
                         ifelse(grepl("var_|total_P_eff", Variable) & Behaviour == "Producing", . * var(df$Producing_events), .))))

Bivar_unadjusted <- Bivar_unadjusted %>% mutate(across(c("2.5%", "50%", "97.5%"), 
                                           ~ ifelse(grepl("order", Variable), . / sd(df$TrialNR_Indiv), .)))

# Wrangle further
Bivar_full <- Bivar_full %>% mutate(across(c("2.5%", "50%", "97.5%"), ~ format(round(., 3), nsmall = 3)))
Bivar_full <- Bivar_full %>% unite(`2.5%`, `97.5%`, col = "Credible interval", sep = " ; ", remove = F)
Bivar_full <- Bivar_full %>% mutate(`Credible interval`= paste0("(",`Credible interval`, ")"))
Bivar_full <- Bivar_full %>% unite(`50%`, `Credible interval`, col = "Estimate", sep = " ", remove = F)
Bivar_full$Variables1 <- c("Intercept", "Year", "Trial order", "Trial day", "Group trial before",
                           "Trial day * Group trial before", "Trial day * Trial order", "Trial order * Group trial before",  "Trial day * Trial order * Group trial before",
                           "Intercept", "Year", "Trial order", "Trial day", "Group trial before",
                           "Trial day * Group trial before", "Trial day * Trial order", "Trial order * Group trial before",  "Trial day * Trial order * Group trial before",
                           "Var DIE", "Var IIE","Var Trial", "Var Group",  "Residual variance", "Total variance",
                           "Var DIE", "Var IIE","Var Trial", "Var Group",  "Residual variance", "Total variance",
                           "R DIE", "R IIE", "R Trial", "R Group", "R Residual",  
                           "R DIE", "R IIE", "R Trial", "R Group", "R Residual", "Var total phenotypic effect", "Var total phenotypic effect","R total phenotypic effect", "R total phenotypic effect", "lp__") 

Prd_bivar <- Bivar_full %>% filter(Behaviour == "Producing") %>% dplyr::select("Variables1", "Estimate")
Scr_bivar <- Bivar_full %>% filter(Behaviour == "Scrounging") %>% dplyr::select("Variables1", "Estimate")

Bivar_unadjusted <- Bivar_unadjusted %>% mutate(across(c("2.5%", "50%", "97.5%"), ~ format(round(., 3), nsmall = 3)))
Bivar_unadjusted <- Bivar_unadjusted %>% unite(`2.5%`, `97.5%`, col = "Credible interval", sep = " ; ", remove = F)
Bivar_unadjusted <- Bivar_unadjusted %>% mutate(`Credible interval`= paste0("(",`Credible interval`, ")"))
Bivar_unadjusted <- Bivar_unadjusted %>% unite(`50%`, `Credible interval`, col = "Estimate", sep = " ", remove = F)
Bivar_unadjusted$Variables1 <- c("Intercept", "Intercept", 
                                 "Var DIE", "Var IIE","Var Trial", "Var Group",  "Residual variance", "Total variance",
                                 "Var DIE", "Var IIE","Var Trial", "Var Group",  "Residual variance", "Total variance",
                                 "R DIE", "R IIE", "R Trial", "R Group", "R Residual",  
                                 "R DIE", "R IIE", "R Trial", "R Group", "R Residual", "Var total phenotypic effect", "Var total phenotypic effect","R total phenotypic effect", "R total phenotypic effect", "lp__")

Prd_bivar_2 <- Bivar_unadjusted %>% filter(Behaviour == "Producing") %>% dplyr::select("Variables1", "Estimate")
Scr_bivar_2 <- Bivar_unadjusted %>% filter(Behaviour == "Scrounging") %>% dplyr::select("Variables1", "Estimate")


# Make table
bivar_table <-  c("Intercept", "Year", "Trial order", "Trial day", "Group trial before",
                  "Trial day * Group trial before", "Trial day * Trial order", "Trial order * Group trial before",  "Trial day * Trial order * Group trial before",
                  "Var DIE", "Var IIE","Var Trial", "Var Group", "Residual variance", "Var total phenotypic effect",  "Total variance",
                  "R DIE", "R IIE", "R Trial", "R Group", "R Residual", "R total phenotypic effect")

bivar_table <- as.data.frame(bivar_table)
bivar_table <- bivar_table %>% rename(Variables1 = bivar_table)
bivar_table$Order <- seq(bivar_table$Variables1) 

bivar_table <- merge(bivar_table, Prd_bivar, by = "Variables1", all.x = T)
bivar_table <- bivar_table %>% rename("Bivar producing" = Estimate)
bivar_table <- merge(bivar_table, Scr_bivar, by = "Variables1", all.x = T)
bivar_table <- bivar_table %>% rename("Bivar scrounging" = Estimate)
bivar_table <- merge(bivar_table, Prd_bivar_2, by = "Variables1", all.x = T)
bivar_table <- bivar_table %>% rename("Bivar producing unadjusted" = Estimate)
bivar_table <- merge(bivar_table, Scr_bivar_2, by = "Variables1", all.x = T)
bivar_table <- bivar_table %>% rename("Bivar scrounging unadjusted" = Estimate)

bivar_table[is.na(bivar_table)] <- "-"
bivar_table <- bivar_table %>% arrange(Order)
bivar_table <- bivar_table %>% dplyr::select(- "Order")
bivar_table <- bivar_table %>% rename(Variables = Variables1)

# Insert header rows:
# Fixed effects
Fixed_effects_row <- data.frame(matrix(c("Fixed effects", rep("\u03B2 (95% CI)", ncol(bivar_table) - 1)), nrow = 1))
colnames(Fixed_effects_row) <- colnames(bivar_table)
bivar_table <- add_row(bivar_table, Fixed_effects_row, .before = which(bivar_table$Variable == "Intercept"))

# Random effects
Random_effects_row <- data.frame(matrix(c("Random effects", rep("\u03C3\u00B2 (95% CI)", ncol(bivar_table) - 1)), nrow = 1))
colnames(Random_effects_row) <- colnames(bivar_table)
bivar_table <- add_row(bivar_table, Random_effects_row, .before = which(bivar_table$Variable == "Var DIE"))

# Repeatability
Repeatability_row <- data.frame(matrix(c("Repeatability", rep("R (95% CI)", ncol(bivar_table) - 1)), nrow = 1))
colnames(Repeatability_row) <- colnames(bivar_table)
bivar_table <- add_row(bivar_table, Repeatability_row, .before = which(bivar_table$Variable == "R DIE"))

# Remove identifying strings (e.g. Var, Cor)
bivar_table$Variables <- str_remove(bivar_table$Variable, pattern = str_c(c("Var ", "variance", "R "), collapse = "|"))

sjPlot::tab_df(bivar_table, alternate.rows = T)

sjPlot::tab_df(bivar_table, alternate.rows = T, file = paste("Plots/Table_bivar.doc"))

#===============================================================================
# Suppl tab3: Model tables for short & cross-year repeatability

# Read in BLUPS from model & wrangle
cy_rep_prd <- as.data.frame(readRDS("Model_output/Sum_md_short_rep_prd.rds"))
cy_rep_scr <- as.data.frame(readRDS("Model_output/Sum_md_short_rep_scr.rds"))

cy_rep_prd$Variable <- row.names(cy_rep_prd)
cy_rep_scr$Variable <- row.names(cy_rep_scr)

# Back-transform 
cy_rep_prd <- cy_rep_prd %>%
  mutate(across(c("2.5%", "50%", "97.5%"), 
                ~ ifelse(grepl("B_|int_|cov", Variable), . * sd(df$Producing_events), 
                         ifelse(grepl("Sigma2_", Variable), . * var(df$Producing_events), .))))

cy_rep_prd <- cy_rep_prd %>% mutate(across(c("2.5%", "50%", "97.5%"), 
                                           ~ ifelse(Variable == "B_0", . + mean(df$Producing_events), .)))

cy_rep_prd <- cy_rep_prd %>% mutate(across(c("2.5%", "50%", "97.5%"), 
                                                         ~ ifelse(grepl("order", Variable), . / sd(df$TrialNR_Indiv), .)))

cy_rep_scr <- cy_rep_scr %>%
  mutate(across(c("2.5%", "50%", "97.5%"), 
                ~ ifelse(grepl("B_|int_|cov", Variable), . * sd(df$Scrounging_events), 
                         ifelse(grepl("Sigma2_", Variable), . * var(df$Scrounging_events), .))))

cy_rep_scr <- cy_rep_scr %>% mutate(across(c("2.5%", "50%", "97.5%"), 
                                           ~ ifelse(Variable == "B_0", . + mean(df$Scrounging_events), .)))

cy_rep_scr <- cy_rep_scr %>% mutate(across(c("2.5%", "50%", "97.5%"), 
                                           ~ ifelse(grepl("order", Variable), . / sd(df$TrialNR_Indiv), .)))

# Wrangle data
cy_rep_prd <- cy_rep_prd %>% mutate(across(c("2.5%", "50%", "97.5%"), ~ format(round(., 3), nsmall = 3)))
Var8 <- cy_rep_prd %>% unite(`2.5%`, `97.5%`, col = "Credible interval", sep = " ; ", remove = F)
Var8 <- Var8 %>% mutate(`Credible interval`= paste0("(",`Credible interval`, ")"))
Var8 <- Var8 %>% unite(`50%`, `Credible interval`, col = "Estimate", sep = " ", remove = F)
Var8$Variables1 <- c("Intercept", "Year", "Trial order", "Trial day", "Group trial before",
                     "Trial day * Group trial before", "Trial day * Trial order", "Trial order * Group trial before",  "Trial day * Trial order * Group trial before",
                     "Var DIE", "Var IIE", "Var Group", "Var Trial", "Var Series DIE", "Var Series IIE",  "Residual variance", "Total variance",
                     "Cov DIE-IIE", "Cor DIE-IIE", "Cov series DIE-IIE", "Cor series DIE-IIE",
                     "R DIE", "R IIE", "R Group", "R Trial", "R Series DIE", "R Series IIE", "R Residual", 
                     "Cross-year rep DIE", "Short-term rep DIE", "Cross-year rep IIE","Short-term rep IIE","lp__")
cy_rep_prd2 <- Var8 %>% select("Variables1", "Estimate")

cy_rep_scr <- cy_rep_scr %>% mutate(across(c("2.5%", "50%", "97.5%"), ~ format(round(., 3), nsmall = 3)))
Var9 <- cy_rep_scr %>% unite(`2.5%`, `97.5%`, col = "Credible interval", sep = " ; ", remove = F)
Var9 <- Var9 %>% mutate(`Credible interval`= paste0("(",`Credible interval`, ")"))
Var9 <- Var9 %>% unite(`50%`, `Credible interval`, col = "Estimate", sep = " ", remove = F)
Var9$Variables1 <- c("Intercept", "Year", "Trial order", "Trial day", "Group trial before",
                     "Trial day * Group trial before", "Trial day * Trial order", "Trial order * Group trial before",  "Trial day * Trial order * Group trial before",
                     "Var DIE", "Var IIE", "Var Group", "Var Trial", "Var Series DIE", "Var Series IIE",  "Residual variance", "Total variance",
                     "Cov DIE-IIE", "Cor DIE-IIE", "Cov series DIE-IIE", "Cor series DIE-IIE",
                     "R DIE", "R IIE", "R Group", "R Trial", "R Series DIE", "R Series IIE", "R Residual", 
                     "Cross-year rep DIE", "Short-term rep DIE", "Cross-year rep IIE","Short-term rep IIE","lp__")
cy_rep_scr2 <- Var9 %>% select("Variables1", "Estimate")

# Make table
cy_rep_table <-  c("Intercept", "Year", "Trial order", "Trial day", "Group trial before",
                   "Trial day * Group trial before", "Trial day * Trial order", "Trial order * Group trial before",  "Trial day * Trial order * Group trial before",
                   "Var DIE", "Var IIE","Var Trial", "Var Group", "Var Series DIE", "Var Series IIE",  "Residual variance", "Total variance",
                   "Cov DIE-IIE",  "Cov series DIE-IIE", 
                   "R DIE", "R IIE", "R Trial", "R Group", "R Series DIE", "R Series IIE", "R Residual", 
                   "Cross-year rep DIE", "Cross-year rep IIE", 
                   "Cor DIE-IIE", "Cor series DIE-IIE")

cy_rep_table <- as.data.frame(cy_rep_table)
cy_rep_table <- cy_rep_table %>% rename(Variables1 = cy_rep_table)
cy_rep_table$Order <- seq(cy_rep_table$Variables1) 

cy_rep_table <- merge(cy_rep_table, cy_rep_prd2, by = "Variables1", all.x = T)
cy_rep_table <- cy_rep_table %>% rename("Cross-year producing" = Estimate)
cy_rep_table <- merge(cy_rep_table, cy_rep_scr2, by = "Variables1", all.x = T)
cy_rep_table <- cy_rep_table %>% rename("Cross-year scrounging" = Estimate)

cy_rep_table[is.na(cy_rep_table)] <- "-"
cy_rep_table <- cy_rep_table %>% arrange(Order)
cy_rep_table <- cy_rep_table %>% select(- "Order")
cy_rep_table <- cy_rep_table %>% rename(Variables = Variables1)

# Insert header rows:
# Fixed effects
Fixed_effects_row <- data.frame(matrix(c("Fixed effects", rep("\u03B2 (95% CI)", ncol(cy_rep_table) - 1)), nrow = 1))
colnames(Fixed_effects_row) <- colnames(cy_rep_table)
cy_rep_table <- add_row(cy_rep_table, Fixed_effects_row, .before = which(cy_rep_table$Variable == "Intercept"))

# Random effects
Random_effects_row <- data.frame(matrix(c("Random effects", rep("\u03C3\u00B2 (95% CI)", ncol(cy_rep_table) - 1)), nrow = 1))
colnames(Random_effects_row) <- colnames(cy_rep_table)
cy_rep_table <- add_row(cy_rep_table, Random_effects_row, .before = which(cy_rep_table$Variable == "Var DIE"))

# Covariance
Covariance_row <- data.frame(matrix(c("Covariance", rep("\u03C3 (95% CI)", ncol(cy_rep_table) - 1)), nrow = 1))
colnames(Covariance_row) <- colnames(cy_rep_table)
cy_rep_table <- add_row(cy_rep_table, Covariance_row, .before = which(cy_rep_table$Variable == "Cov DIE-IIE"))

# Repeatability
Repeatability_row <- data.frame(matrix(c("Repeatability", rep("R (95% CI)", ncol(cy_rep_table) - 1)), nrow = 1))
colnames(Repeatability_row) <- colnames(cy_rep_table)
cy_rep_table <- add_row(cy_rep_table, Repeatability_row, .before = which(cy_rep_table$Variable == "R DIE"))

# Correlations
Correlations_row <- data.frame(matrix(c("Correlations", rep("r (95% CI)", ncol(cy_rep_table) - 1)), nrow = 1))
colnames(Correlations_row) <- colnames(cy_rep_table)
cy_rep_table <- add_row(cy_rep_table, Correlations_row, .before = which(cy_rep_table$Variable == "Cor DIE-IIE"))

# Remove identifying strings (e.g. Var, Cor)
cy_rep_table$Variables <- str_remove(cy_rep_table$Variable, pattern = str_c(c("Var ", "variance", "Cov ", "R ", "Cor "), collapse = "|"))

sjPlot::tab_df(cy_rep_table, alternate.rows = T)

sjPlot::tab_df(cy_rep_table, alternate.rows = T, file = paste("Plots/Table_cy_rep.doc"))

#===============================================================================
# Suppl tab4 & 5: Table for all split models

# Load in the raw data for means and standard
df <- read.csv("Data/Data_PS_DIE_IIE_cov.csv", header = TRUE)

# Load in model output

# Producing
Sum_md_prd_22 <- readRDS("Model_output/Split_models/Sum_md_prd_22.rds")
Sum_md_prd_22_td1 <- readRDS("Model_output/Split_models/Sum_md_prd_22_td1.rds")
Sum_md_prd_22_td2 <- readRDS("Model_output/Split_models/Sum_md_prd_22_td2.rds")
Sum_md_prd_23 <- readRDS("Model_output/Split_models/Sum_md_prd_23.rds")
Sum_md_prd_23_td1 <- readRDS("Model_output/Split_models/Sum_md_prd_23_td1.rds")
Sum_md_prd_23_td2 <- readRDS("Model_output/Split_models/Sum_md_prd_23_td2.rds")
Sum_md_prd_full <- readRDS("Model_output/Split_models/Sum_md_prd_full.rds")

# Scrounging
Sum_md_scr_22 <- readRDS("Model_output/Split_models/Sum_md_scr_22.rds")
Sum_md_scr_22_td1 <- readRDS("Model_output/Split_models/Sum_md_scr_22_td1.rds")
Sum_md_scr_22_td2 <- readRDS("Model_output/Split_models/Sum_md_scr_22_td2.rds")
Sum_md_scr_23 <- readRDS("Model_output/Split_models/Sum_md_scr_23.rds")
Sum_md_scr_23_td1 <- readRDS("Model_output/Split_models/Sum_md_scr_23_td1.rds")
Sum_md_scr_23_td2 <- readRDS("Model_output/Split_models/Sum_md_scr_23_td2.rds")
Sum_md_scr_full <- readRDS("Model_output/Split_models/Sum_md_scr_full.rds")

# Split up datasets for means and sd

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

# Back transform

# scr_2022
Sum_md_scr_22 <- as.data.frame(Sum_md_scr_22)
Sum_md_scr_22$Variable <- row.names(Sum_md_scr_22)

Sum_md_scr_22 <- Sum_md_scr_22 %>%
  mutate(across(c("2.5%", "50%", "97.5%"), 
                ~ ifelse(grepl("B_|int_|cov", Variable), . * sd(df_22$Scrounging_events), 
                         ifelse(grepl("Sigma2_", Variable), . * var(df_22$Scrounging_events), .))))

Sum_md_scr_22 <- Sum_md_scr_22 %>% mutate(across(c("2.5%", "50%", "97.5%"), 
                                                 ~ ifelse(Variable == "B_0", . + mean(df_22$Scrounging_events), .)))

Sum_md_scr_22 <- Sum_md_scr_22 %>% mutate(across(c("2.5%", "50%", "97.5%"), 
                                                 ~ ifelse(grepl("order", Variable), . / sd(df_22$TrialNR_Indiv), .)))

# scr_2022_td1
Sum_md_scr_22_td1 <- as.data.frame(Sum_md_scr_22_td1)
Sum_md_scr_22_td1$Variable <- row.names(Sum_md_scr_22_td1)

Sum_md_scr_22_td1 <- Sum_md_scr_22_td1 %>%
  mutate(across(c("2.5%", "50%", "97.5%"), 
                ~ ifelse(grepl("B_|int_|cov", Variable), . * sd(df_22_td1$Scrounging_events), 
                         ifelse(grepl("Sigma2_", Variable), . * var(df_22_td1$Scrounging_events), .))))

Sum_md_scr_22_td1 <- Sum_md_scr_22_td1 %>% mutate(across(c("2.5%", "50%", "97.5%"), 
                                                         ~ ifelse(Variable == "B_0", . + mean(df_22_td1$Scrounging_events), .)))

Sum_md_scr_22_td1 <- Sum_md_scr_22_td1 %>% mutate(across(c("2.5%", "50%", "97.5%"), 
                                                         ~ ifelse(grepl("order", Variable), . / sd(df_22_td1$TrialNR_Indiv), .)))

# scr_2022_td2
Sum_md_scr_22_td2 <- as.data.frame(Sum_md_scr_22_td2)
Sum_md_scr_22_td2$Variable <- row.names(Sum_md_scr_22_td2)

Sum_md_scr_22_td2 <- Sum_md_scr_22_td2 %>%
  mutate(across(c("2.5%", "50%", "97.5%"), 
                ~ ifelse(grepl("B_|int_|cov", Variable), . * sd(df_22_td2$Scrounging_events), 
                         ifelse(grepl("Sigma2_", Variable), . * var(df_22_td2$Scrounging_events), .))))

Sum_md_scr_22_td2 <- Sum_md_scr_22_td2 %>% mutate(across(c("2.5%", "50%", "97.5%"), 
                                                         ~ ifelse(Variable == "B_0", . + mean(df_22_td2$Scrounging_events), .)))

Sum_md_scr_22_td2 <- Sum_md_scr_22_td2 %>% mutate(across(c("2.5%", "50%", "97.5%"), 
                                                         ~ ifelse(grepl("order", Variable), . / sd(df_22_td2$TrialNR_Indiv), .)))

# scr_2023
Sum_md_scr_23 <- as.data.frame(Sum_md_scr_23)
Sum_md_scr_23$Variable <- row.names(Sum_md_scr_23)

Sum_md_scr_23 <- Sum_md_scr_23 %>%
  mutate(across(c("2.5%", "50%", "97.5%"), 
                ~ ifelse(grepl("B_|int_|cov", Variable), . * sd(df_23$Scrounging_events), 
                         ifelse(grepl("Sigma2_", Variable), . * var(df_23$Scrounging_events), .))))

Sum_md_scr_23 <- Sum_md_scr_23 %>% mutate(across(c("2.5%", "50%", "97.5%"), 
                                                 ~ ifelse(Variable == "B_0", . + mean(df_23$Scrounging_events), .)))

Sum_md_scr_23 <- Sum_md_scr_23 %>% mutate(across(c("2.5%", "50%", "97.5%"), 
                                                 ~ ifelse(grepl("order", Variable), . / sd(df_23$TrialNR_Indiv), .)))

# scr_2023_td1
Sum_md_scr_23_td1 <- as.data.frame(Sum_md_scr_23_td1)
Sum_md_scr_23_td1$Variable <- row.names(Sum_md_scr_23_td1)

Sum_md_scr_23_td1 <- Sum_md_scr_23_td1 %>%
  mutate(across(c("2.5%", "50%", "97.5%"), 
                ~ ifelse(grepl("B_|int_|cov", Variable), . * sd(df_23_td1$Scrounging_events), 
                         ifelse(grepl("Sigma2_", Variable), . * var(df_23_td1$Scrounging_events), .))))

Sum_md_scr_23_td1 <- Sum_md_scr_23_td1 %>% mutate(across(c("2.5%", "50%", "97.5%"), 
                                                         ~ ifelse(Variable == "B_0", . + mean(df_23_td1$Scrounging_events), .)))

Sum_md_scr_23_td1 <- Sum_md_scr_23_td1 %>% mutate(across(c("2.5%", "50%", "97.5%"), 
                                                         ~ ifelse(grepl("order", Variable), . / sd(df_23_td1$TrialNR_Indiv), .)))

# scr_2023_td2
Sum_md_scr_23_td2 <- as.data.frame(Sum_md_scr_23_td2)
Sum_md_scr_23_td2$Variable <- row.names(Sum_md_scr_23_td2)

Sum_md_scr_23_td2 <- Sum_md_scr_23_td2 %>%
  mutate(across(c("2.5%", "50%", "97.5%"), 
                ~ ifelse(grepl("B_|int_|cov", Variable), . * sd(df_23_td2$Scrounging_events), 
                         ifelse(grepl("Sigma2_", Variable), . * var(df_23_td2$Scrounging_events), .))))

Sum_md_scr_23_td2 <- Sum_md_scr_23_td2 %>% mutate(across(c("2.5%", "50%", "97.5%"), 
                                                         ~ ifelse(Variable == "B_0", . + mean(df_23_td2$Scrounging_events), .)))

Sum_md_scr_23_td2 <- Sum_md_scr_23_td2 %>% mutate(across(c("2.5%", "50%", "97.5%"), 
                                                         ~ ifelse(grepl("order", Variable), . / sd(df_23_td2$TrialNR_Indiv), .)))

# scr_full
Sum_md_scr_full <- as.data.frame(Sum_md_scr_full)
Sum_md_scr_full$Variable <- row.names(Sum_md_scr_full)

Sum_md_scr_full <- Sum_md_scr_full %>%
  mutate(across(c("2.5%", "50%", "97.5%"), 
                ~ ifelse(grepl("B_|int_|cov", Variable), . * sd(df$Scrounging_events), 
                         ifelse(grepl("Sigma2_", Variable), . * var(df$Scrounging_events), .))))

Sum_md_scr_full <- Sum_md_scr_full %>% mutate(across(c("2.5%", "50%", "97.5%"), 
                                                     ~ ifelse(Variable == "B_0", . + mean(df$Scrounging_events), .)))

Sum_md_scr_full <- Sum_md_scr_full %>% mutate(across(c("2.5%", "50%", "97.5%"), 
                                                     ~ ifelse(grepl("order", Variable), . / sd(df$TrialNR_Indiv), .)))


# prd_2022
Sum_md_prd_22 <- as.data.frame(Sum_md_prd_22)
Sum_md_prd_22$Variable <- row.names(Sum_md_prd_22)

Sum_md_prd_22 <- Sum_md_prd_22 %>%
  mutate(across(c("2.5%", "50%", "97.5%"), 
                ~ ifelse(grepl("B_|int_|cov", Variable), . * sd(df_22$Producing_events), 
                         ifelse(grepl("Sigma2_", Variable), . * var(df_22$Producing_events), .))))

Sum_md_prd_22 <- Sum_md_prd_22 %>% mutate(across(c("2.5%", "50%", "97.5%"), 
                                           ~ ifelse(Variable == "B_0", . + mean(df_22$Producing_events), .)))

Sum_md_prd_22 <- Sum_md_prd_22 %>% mutate(across(c("2.5%", "50%", "97.5%"), 
                                                 ~ ifelse(grepl("order", Variable), . / sd(df_22$TrialNR_Indiv), .)))

# prd_2022_td1
Sum_md_prd_22_td1 <- as.data.frame(Sum_md_prd_22_td1)
Sum_md_prd_22_td1$Variable <- row.names(Sum_md_prd_22_td1)

Sum_md_prd_22_td1 <- Sum_md_prd_22_td1 %>%
  mutate(across(c("2.5%", "50%", "97.5%"), 
                ~ ifelse(grepl("B_|int_|cov", Variable), . * sd(df_22_td1$Producing_events), 
                         ifelse(grepl("Sigma2_", Variable), . * var(df_22_td1$Producing_events), .))))

Sum_md_prd_22_td1 <- Sum_md_prd_22_td1 %>% mutate(across(c("2.5%", "50%", "97.5%"), 
                                                 ~ ifelse(Variable == "B_0", . + mean(df_22_td1$Producing_events), .)))

Sum_md_prd_22_td1 <- Sum_md_prd_22_td1 %>% mutate(across(c("2.5%", "50%", "97.5%"), 
                                                 ~ ifelse(grepl("order", Variable), . / sd(df_22_td1$TrialNR_Indiv), .)))

# prd_2022_td2
Sum_md_prd_22_td2 <- as.data.frame(Sum_md_prd_22_td2)
Sum_md_prd_22_td2$Variable <- row.names(Sum_md_prd_22_td2)

Sum_md_prd_22_td2 <- Sum_md_prd_22_td2 %>%
  mutate(across(c("2.5%", "50%", "97.5%"), 
                ~ ifelse(grepl("B_|int_|cov", Variable), . * sd(df_22_td2$Producing_events), 
                         ifelse(grepl("Sigma2_", Variable), . * var(df_22_td2$Producing_events), .))))

Sum_md_prd_22_td2 <- Sum_md_prd_22_td2 %>% mutate(across(c("2.5%", "50%", "97.5%"), 
                                                         ~ ifelse(Variable == "B_0", . + mean(df_22_td2$Producing_events), .)))

Sum_md_prd_22_td2 <- Sum_md_prd_22_td2 %>% mutate(across(c("2.5%", "50%", "97.5%"), 
                                                         ~ ifelse(grepl("order", Variable), . / sd(df_22_td2$TrialNR_Indiv), .)))

# prd_2023
Sum_md_prd_23 <- as.data.frame(Sum_md_prd_23)
Sum_md_prd_23$Variable <- row.names(Sum_md_prd_23)

Sum_md_prd_23 <- Sum_md_prd_23 %>%
  mutate(across(c("2.5%", "50%", "97.5%"), 
                ~ ifelse(grepl("B_|int_|cov", Variable), . * sd(df_23$Producing_events), 
                         ifelse(grepl("Sigma2_", Variable), . * var(df_23$Producing_events), .))))

Sum_md_prd_23 <- Sum_md_prd_23 %>% mutate(across(c("2.5%", "50%", "97.5%"), 
                                                 ~ ifelse(Variable == "B_0", . + mean(df_23$Producing_events), .)))

Sum_md_prd_23 <- Sum_md_prd_23 %>% mutate(across(c("2.5%", "50%", "97.5%"), 
                                                 ~ ifelse(grepl("order", Variable), . / sd(df_23$TrialNR_Indiv), .)))

# prd_2023_td1
Sum_md_prd_23_td1 <- as.data.frame(Sum_md_prd_23_td1)
Sum_md_prd_23_td1$Variable <- row.names(Sum_md_prd_23_td1)

Sum_md_prd_23_td1 <- Sum_md_prd_23_td1 %>%
  mutate(across(c("2.5%", "50%", "97.5%"), 
                ~ ifelse(grepl("B_|int_|cov", Variable), . * sd(df_23_td1$Producing_events), 
                         ifelse(grepl("Sigma2_", Variable), . * var(df_23_td1$Producing_events), .))))

Sum_md_prd_23_td1 <- Sum_md_prd_23_td1 %>% mutate(across(c("2.5%", "50%", "97.5%"), 
                                                         ~ ifelse(Variable == "B_0", . + mean(df_23_td1$Producing_events), .)))

Sum_md_prd_23_td1 <- Sum_md_prd_23_td1 %>% mutate(across(c("2.5%", "50%", "97.5%"), 
                                                         ~ ifelse(grepl("order", Variable), . / sd(df_23_td1$TrialNR_Indiv), .)))

# prd_2023_td2
Sum_md_prd_23_td2 <- as.data.frame(Sum_md_prd_23_td2)
Sum_md_prd_23_td2$Variable <- row.names(Sum_md_prd_23_td2)

Sum_md_prd_23_td2 <- Sum_md_prd_23_td2 %>%
  mutate(across(c("2.5%", "50%", "97.5%"), 
                ~ ifelse(grepl("B_|int_|cov", Variable), . * sd(df_23_td2$Producing_events), 
                         ifelse(grepl("Sigma2_", Variable), . * var(df_23_td2$Producing_events), .))))

Sum_md_prd_23_td2 <- Sum_md_prd_23_td2 %>% mutate(across(c("2.5%", "50%", "97.5%"), 
                                                         ~ ifelse(Variable == "B_0", . + mean(df_23_td2$Producing_events), .)))

Sum_md_prd_23_td2 <- Sum_md_prd_23_td2 %>% mutate(across(c("2.5%", "50%", "97.5%"), 
                                                         ~ ifelse(grepl("order", Variable), . / sd(df_23_td2$TrialNR_Indiv), .)))

# prd_full
Sum_md_prd_full <- as.data.frame(Sum_md_prd_full)
Sum_md_prd_full$Variable <- row.names(Sum_md_prd_full)

Sum_md_prd_full <- Sum_md_prd_full %>%
  mutate(across(c("2.5%", "50%", "97.5%"), 
                ~ ifelse(grepl("B_|int_|cov", Variable), . * sd(df$Producing_events), 
                         ifelse(grepl("Sigma2_", Variable), . * var(df$Producing_events), .))))

Sum_md_prd_full <- Sum_md_prd_full %>% mutate(across(c("2.5%", "50%", "97.5%"), 
                                                 ~ ifelse(Variable == "B_0", . + mean(df$Producing_events), .)))

Sum_md_prd_full <- Sum_md_prd_full %>% mutate(across(c("2.5%", "50%", "97.5%"), 
                                                 ~ ifelse(grepl("order", Variable), . / sd(df$TrialNR_Indiv), .)))


# Wrangle data for table
Var1 <- as.data.frame(Sum_md_scr_22)
Var1 <- Var1 %>% mutate(across(c("2.5%", "50%", "97.5%"), ~ format(round(., 3), nsmall = 3)))
Var1 <- Var1 %>% unite(`2.5%`, `97.5%`, col = "Credible interval", sep = " ; ", remove = F)
Var1 <- Var1 %>% mutate(`Credible interval`= paste0("(",`Credible interval`, ")"))
Var1 <- Var1 %>% unite(`50%`, `Credible interval`, col = "Estimate", sep = " ", remove = F)
Var1$Variables1 <- c("Intercept scrounging", "BMR", "Trial day", "Trial order", "Trial day * Trial order", "Var DIE", "Var IIE", "Var Trial", "Var Group",  "Residual variance", "Total variance", "Cov DIE-IIE", "Cor DIE-IIE",
                     "R DIE", "R IIE", "R Group", "R Trial",  "R Residual", "lp__")
scr_mod1 <- Var1 %>% select("Variables1", "Estimate")

Var2 <- as.data.frame(Sum_md_scr_22_td1)
Var2 <- Var2 %>% mutate(across(c("2.5%", "50%", "97.5%"), ~ format(round(., 3), nsmall = 3)))
Var2 <- Var2 %>% unite(`2.5%`, `97.5%`, col = "Credible interval", sep = " ; ", remove = F)
Var2 <- Var2 %>% mutate(`Credible interval`= paste0("(",`Credible interval`, ")"))
Var2 <- Var2 %>% unite(`50%`, `Credible interval`, col = "Estimate", sep = " ", remove = F)
Var2$Variables1 <- c("Intercept scrounging", "BMR", "Trial order", "Var DIE", "Var IIE", "Var Trial", "Var Group", "Residual variance", "Total variance", "Cov DIE-IIE", "Cor DIE-IIE",
                      "R DIE", "R IIE", "R Group", "R Trial",  "R Residual", "lp__")
scr_mod2 <- Var2 %>% select("Variables1", "Estimate")

Var3 <- as.data.frame(Sum_md_scr_22_td2)
Var3 <- Var3 %>% mutate(across(c("2.5%", "50%", "97.5%"), ~ format(round(., 3), nsmall = 3)))
Var3 <- Var3 %>% unite(`2.5%`, `97.5%`, col = "Credible interval", sep = " ; ", remove = F)
Var3 <- Var3 %>% mutate(`Credible interval`= paste0("(",`Credible interval`, ")"))
Var3 <- Var3 %>% unite(`50%`, `Credible interval`, col = "Estimate", sep = " ", remove = F)
Var3$Variables1 <- c("Intercept scrounging", "BMR", "Trial order", "Var DIE", "Var IIE", "Var Trial", "Var Group", "Residual variance", "Total variance", "Cov DIE-IIE", "Cor DIE-IIE",
                      "R DIE", "R IIE", "R Group", "R Trial",  "R Residual", "lp__")
scr_mod3 <- Var3 %>% select("Variables1", "Estimate")

Var4 <- as.data.frame(Sum_md_scr_23)
Var4 <- Var4 %>% mutate(across(c("2.5%", "50%", "97.5%"), ~ format(round(., 3), nsmall = 3)))
Var4 <- Var4 %>% unite(`2.5%`, `97.5%`, col = "Credible interval", sep = " ; ", remove = F)
Var4 <- Var4 %>% mutate(`Credible interval`= paste0("(",`Credible interval`, ")"))
Var4 <- Var4 %>% unite(`50%`, `Credible interval`, col = "Estimate", sep = " ", remove = F)
Var4$Variables1 <- c("Intercept scrounging", "Trial order", "Trial day", "Group trial before", "Trial day * Trial order", "Trial day * Group trial before", "Trial order * Group trial before", "Trial day * Trial order * Group trial before",  "Var DIE", "Var IIE", "Var Trial", "Var Group",  "Residual variance", "Total variance", "Cov DIE-IIE", "Cor DIE-IIE",
                      "R DIE", "R IIE", "R Group", "R Trial",  "R Residual", "lp__")
scr_mod4 <- Var4 %>% select("Variables1", "Estimate")

Var5 <- as.data.frame(Sum_md_scr_23_td1)
Var5 <- Var5 %>% mutate(across(c("2.5%", "50%", "97.5%"), ~ format(round(., 3), nsmall = 3)))
Var5 <- Var5 %>% unite(`2.5%`, `97.5%`, col = "Credible interval", sep = " ; ", remove = F)
Var5 <- Var5 %>% mutate(`Credible interval`= paste0("(",`Credible interval`, ")"))
Var5 <- Var5 %>% unite(`50%`, `Credible interval`, col = "Estimate", sep = " ", remove = F)
Var5$Variables1 <- c("Intercept scrounging", "Trial order", "Group trial before", "Trial order * Group trial before", "Var DIE", "Var IIE", "Var Trial", "Var Group",  "Residual variance", "Total variance", "Cov DIE-IIE", "Cor DIE-IIE",
                      "R DIE", "R IIE", "R Group", "R Trial",  "R Residual", "lp__")
scr_mod5 <- Var5 %>% select("Variables1", "Estimate")

Var6 <- as.data.frame(Sum_md_scr_23_td2)
Var6 <- Var6 %>% mutate(across(c("2.5%", "50%", "97.5%"), ~ format(round(., 3), nsmall = 3)))
Var6 <- Var6 %>% unite(`2.5%`, `97.5%`, col = "Credible interval", sep = " ; ", remove = F)
Var6 <- Var6 %>% mutate(`Credible interval`= paste0("(",`Credible interval`, ")"))
Var6 <- Var6 %>% unite(`50%`, `Credible interval`, col = "Estimate", sep = " ", remove = F)
Var6$Variables1 <- c("Intercept scrounging", "Trial order", "Group trial before", "Trial order * Group trial before", "Var DIE", "Var IIE", "Var Trial", "Var Group",  "Residual variance", "Total variance", "Cov DIE-IIE", "Cor DIE-IIE",
                      "R DIE", "R IIE", "R Group", "R Trial",  "R Residual", "lp__")
scr_mod6 <- Var6 %>% select("Variables1", "Estimate")

Var7 <- as.data.frame(Sum_md_scr_full)
Var7 <- Var7 %>% mutate(across(c("2.5%", "50%", "97.5%"), ~ format(round(., 3), nsmall = 3)))
Var7 <- Var7 %>% unite(`2.5%`, `97.5%`, col = "Credible interval", sep = " ; ", remove = F)
Var7 <- Var7 %>% mutate(`Credible interval`= paste0("(",`Credible interval`, ")"))
Var7 <- Var7 %>% unite(`50%`, `Credible interval`, col = "Estimate", sep = " ", remove = F)
Var7$Variables1 <- c("Intercept scrounging", "Year", "Trial order", "Trial day", "Group trial before", "Trial day * Trial order", "Trial day * Group trial before", "Trial order * Group trial before", "Trial day * Trial order * Group trial before", "Var DIE", "Var IIE", "Var Trial", "Var Group",  "Residual variance", "Total variance", "Cov DIE-IIE", "Cor DIE-IIE",
                      "R DIE", "R IIE", "R Group", "R Trial",  "R Residual", "lp__")
scr_mod7 <- Var7 %>% select("Variables1", "Estimate")

#==================================
# Make table
Scr_table <-  c("Intercept scrounging", "Year", "BMR", "Trial order", "Trial day", "Group trial before",
                "Trial day * Group trial before","Trial order * Group trial before", "Trial day * Trial order", "Trial day * Trial order * Group trial before",
                "Var DIE", "Var IIE", "Var Group", "Var Trial", "Residual variance", "Total variance",
                "Cov DIE-IIE", 
                "R DIE", "R IIE", "R Group", "R Trial", "R Residual",
                "Cor DIE-IIE")

Scr_table <- as.data.frame(Scr_table)
Scr_table <- Scr_table %>% rename(Variables1 = Scr_table)
Scr_table$Order <- seq(Scr_table$Variables1) 

Scr_table <- merge(Scr_table, scr_mod1, by = "Variables1", all.x = T)
Scr_table <- Scr_table %>% rename("mod 2022" = Estimate)
Scr_table <- merge(Scr_table, scr_mod2, by = "Variables1", all.x = T)
Scr_table <- Scr_table %>% rename("mod 2022 day1" = Estimate)
Scr_table <- merge(Scr_table, scr_mod3, by = "Variables1", all.x = T)
Scr_table <- Scr_table %>% rename("mod 2022 day2" = Estimate)
Scr_table <- merge(Scr_table, scr_mod4, by = "Variables1", all.x = T)
Scr_table <- Scr_table %>% rename("mod 2023" = Estimate)
Scr_table <- merge(Scr_table, scr_mod5, by = "Variables1", all.x = T)
Scr_table <- Scr_table %>% rename("mod 2023 day1" = Estimate)
Scr_table <- merge(Scr_table, scr_mod6, by = "Variables1", all.x = T)
Scr_table <- Scr_table %>% rename("mod 2023 day2" = Estimate)
Scr_table <- merge(Scr_table, scr_mod7, by = "Variables1", all.x = T)
Scr_table <- Scr_table %>% rename("mod 2022/2023" = Estimate)

Scr_table[is.na(Scr_table)] <- "-"
Scr_table <- Scr_table %>% arrange(Order)
Scr_table <- Scr_table %>% select(- "Order")
Scr_table <- Scr_table %>% rename(Variables = Variables1)

# Insert header rows:
# Fixed effects
Fixed_effects_row <- data.frame(matrix(c("Fixed effects", rep("\u03B2 (95% CI)", ncol(Scr_table) - 1)), nrow = 1))
colnames(Fixed_effects_row) <- colnames(Scr_table)
Scr_table <- add_row(Scr_table, Fixed_effects_row, .before = which(Scr_table$Variable == "Intercept scrounging"))

# Random effects
Random_effects_row <- data.frame(matrix(c("Random effects", rep("\u03C3\u00B2 (95% CI)", ncol(Scr_table) - 1)), nrow = 1))
colnames(Random_effects_row) <- colnames(Scr_table)
Scr_table <- add_row(Scr_table, Random_effects_row, .before = which(Scr_table$Variable == "Var DIE"))

# Covariance
Covariance_row <- data.frame(matrix(c("Covariance", rep("\u03C3 (95% CI)", ncol(Scr_table) - 1)), nrow = 1))
colnames(Covariance_row) <- colnames(Scr_table)
Scr_table <- add_row(Scr_table, Covariance_row, .before = which(Scr_table$Variable == "Cov DIE-IIE"))

# Repeatability
Repeatability_row <- data.frame(matrix(c("Repeatability", rep("R (95% CI)", ncol(Scr_table) - 1)), nrow = 1))
colnames(Repeatability_row) <- colnames(Scr_table)
Scr_table <- add_row(Scr_table, Repeatability_row, .before = which(Scr_table$Variable == "R DIE"))

# Correlations
Correlations_row <- data.frame(matrix(c("Correlations", rep("r (95% CI)", ncol(Scr_table) - 1)), nrow = 1))
colnames(Correlations_row) <- colnames(Scr_table)
Scr_table <- add_row(Scr_table, Correlations_row, .before = which(Scr_table$Variable == "Cor DIE-IIE"))

# Remove identifying strings (e.g. Var, Cor)
Scr_table$Variables <- str_remove(Scr_table$Variable, pattern = str_c(c("Var ", "variance", "Cov ", "R ", "Cor "), collapse = "|"))

sjPlot::tab_df(Scr_table, alternate.rows = T)

sjPlot::tab_df(Scr_table, alternate.rows = T, file = paste("Plots/Table_scr_split.doc"))

#=======================================
# For producing
# Wrangle data for table
Var1 <- as.data.frame(Sum_md_prd_22)
Var1 <- Var1 %>% mutate(across(c("2.5%", "50%", "97.5%"), ~ format(round(., 3), nsmall = 3)))
Var1 <- Var1 %>% unite(`2.5%`, `97.5%`, col = "Credible interval", sep = " ; ", remove = F)
Var1 <- Var1 %>% mutate(`Credible interval`= paste0("(",`Credible interval`, ")"))
Var1 <- Var1 %>% unite(`50%`, `Credible interval`, col = "Estimate", sep = " ", remove = F)
Var1$Variables1 <- c("Intercept producing", "BMR", "Trial day", "Trial order", "Trial day * Trial order", "Var DIE", "Var IIE", "Var Trial", "Var Group",  "Residual variance", "Total variance", "Cov DIE-IIE", "Cor DIE-IIE",
                     "R DIE", "R IIE", "R Group", "R Trial",  "R Residual", "lp__")
prd_mod1 <- Var1 %>% select("Variables1", "Estimate")

Var2 <- as.data.frame(Sum_md_prd_22_td1)
Var2 <- Var2 %>% mutate(across(c("2.5%", "50%", "97.5%"), ~ format(round(., 3), nsmall = 3)))
Var2 <- Var2 %>% unite(`2.5%`, `97.5%`, col = "Credible interval", sep = " ; ", remove = F)
Var2 <- Var2 %>% mutate(`Credible interval`= paste0("(",`Credible interval`, ")"))
Var2 <- Var2 %>% unite(`50%`, `Credible interval`, col = "Estimate", sep = " ", remove = F)
Var2$Variables1 <- c("Intercept producing", "BMR", "Trial order", "Var DIE", "Var IIE", "Var Trial", "Var Group", "Residual variance", "Total variance", "Cov DIE-IIE", "Cor DIE-IIE",
                     "R DIE", "R IIE", "R Group", "R Trial",  "R Residual", "lp__")
prd_mod2 <- Var2 %>% select("Variables1", "Estimate")

Var3 <- as.data.frame(Sum_md_prd_22_td2)
Var3 <- Var3 %>% mutate(across(c("2.5%", "50%", "97.5%"), ~ format(round(., 3), nsmall = 3)))
Var3 <- Var3 %>% unite(`2.5%`, `97.5%`, col = "Credible interval", sep = " ; ", remove = F)
Var3 <- Var3 %>% mutate(`Credible interval`= paste0("(",`Credible interval`, ")"))
Var3 <- Var3 %>% unite(`50%`, `Credible interval`, col = "Estimate", sep = " ", remove = F)
Var3$Variables1 <- c("Intercept producing", "BMR", "Trial order", "Var DIE", "Var IIE", "Var Trial", "Var Group", "Residual variance", "Total variance", "Cov DIE-IIE", "Cor DIE-IIE",
                     "R DIE", "R IIE", "R Group", "R Trial",  "R Residual", "lp__")
prd_mod3 <- Var3 %>% select("Variables1", "Estimate")

Var4 <- as.data.frame(Sum_md_prd_23)
Var4 <- Var4 %>% mutate(across(c("2.5%", "50%", "97.5%"), ~ format(round(., 3), nsmall = 3)))
Var4 <- Var4 %>% unite(`2.5%`, `97.5%`, col = "Credible interval", sep = " ; ", remove = F)
Var4 <- Var4 %>% mutate(`Credible interval`= paste0("(",`Credible interval`, ")"))
Var4 <- Var4 %>% unite(`50%`, `Credible interval`, col = "Estimate", sep = " ", remove = F)
Var4$Variables1 <- c("Intercept producing", "Trial order", "Trial day", "Group trial before", "Trial day * Trial order", "Trial day * Group trial before", "Trial order * Group trial before", "Trial day * Trial order * Group trial before",  "Var DIE", "Var IIE", "Var Trial", "Var Group",  "Residual variance", "Total variance", "Cov DIE-IIE", "Cor DIE-IIE",
                     "R DIE", "R IIE", "R Group", "R Trial",  "R Residual", "lp__")
prd_mod4 <- Var4 %>% select("Variables1", "Estimate")

Var5 <- as.data.frame(Sum_md_prd_23_td1)
Var5 <- Var5 %>% mutate(across(c("2.5%", "50%", "97.5%"), ~ format(round(., 3), nsmall = 3)))
Var5 <- Var5 %>% unite(`2.5%`, `97.5%`, col = "Credible interval", sep = " ; ", remove = F)
Var5 <- Var5 %>% mutate(`Credible interval`= paste0("(",`Credible interval`, ")"))
Var5 <- Var5 %>% unite(`50%`, `Credible interval`, col = "Estimate", sep = " ", remove = F)
Var5$Variables1 <- c("Intercept producing", "Trial order", "Group trial before", "Trial order * Group trial before", "Var DIE", "Var IIE", "Var Trial", "Var Group",  "Residual variance", "Total variance", "Cov DIE-IIE", "Cor DIE-IIE",
                     "R DIE", "R IIE", "R Group", "R Trial",  "R Residual", "lp__")
prd_mod5 <- Var5 %>% select("Variables1", "Estimate")

Var6 <- as.data.frame(Sum_md_prd_23_td2)
Var6 <- Var6 %>% mutate(across(c("2.5%", "50%", "97.5%"), ~ format(round(., 3), nsmall = 3)))
Var6 <- Var6 %>% unite(`2.5%`, `97.5%`, col = "Credible interval", sep = " ; ", remove = F)
Var6 <- Var6 %>% mutate(`Credible interval`= paste0("(",`Credible interval`, ")"))
Var6 <- Var6 %>% unite(`50%`, `Credible interval`, col = "Estimate", sep = " ", remove = F)
Var6$Variables1 <- c("Intercept producing", "Trial order", "Group trial before", "Trial order * Group trial before", "Var DIE", "Var IIE", "Var Trial", "Var Group",  "Residual variance", "Total variance", "Cov DIE-IIE", "Cor DIE-IIE",
                     "R DIE", "R IIE", "R Group", "R Trial",  "R Residual", "lp__")
prd_mod6 <- Var6 %>% select("Variables1", "Estimate")

Var7 <- as.data.frame(Sum_md_prd_full)
Var7 <- Var7 %>% mutate(across(c("2.5%", "50%", "97.5%"), ~ format(round(., 3), nsmall = 3)))
Var7 <- Var7 %>% unite(`2.5%`, `97.5%`, col = "Credible interval", sep = " ; ", remove = F)
Var7 <- Var7 %>% mutate(`Credible interval`= paste0("(",`Credible interval`, ")"))
Var7 <- Var7 %>% unite(`50%`, `Credible interval`, col = "Estimate", sep = " ", remove = F)
Var7$Variables1 <- c("Intercept producing", "Year", "Trial order", "Trial day", "Group trial before", "Trial day * Trial order", "Trial day * Group trial before", "Trial order * Group trial before", "Trial day * Trial order * Group trial before", "Var DIE", "Var IIE", "Var Trial", "Var Group",  "Residual variance", "Total variance", "Cov DIE-IIE", "Cor DIE-IIE",
                     "R DIE", "R IIE", "R Group", "R Trial",  "R Residual", "lp__")
prd_mod7 <- Var7 %>% select("Variables1", "Estimate")

#==================================
# Make table
prd_table <-  c("Intercept producing", "Year", "BMR", "Trial order", "Trial day", "Group trial before",
                "Trial day * Group trial before","Trial order * Group trial before", "Trial day * Trial order", "Trial day * Trial order * Group trial before",
                "Var DIE", "Var IIE", "Var Group", "Var Trial", "Residual variance", "Total variance",
                "Cov DIE-IIE", 
                "R DIE", "R IIE", "R Group", "R Trial", "R Residual",
                "Cor DIE-IIE")

prd_table <- as.data.frame(prd_table)
prd_table <- prd_table %>% rename(Variables1 = prd_table)
prd_table$Order <- seq(prd_table$Variables1) 

prd_table <- merge(prd_table, prd_mod1, by = "Variables1", all.x = T)
prd_table <- prd_table %>% rename("mod 2022" = Estimate)
prd_table <- merge(prd_table, prd_mod2, by = "Variables1", all.x = T)
prd_table <- prd_table %>% rename("mod 2022 day1" = Estimate)
prd_table <- merge(prd_table, prd_mod3, by = "Variables1", all.x = T)
prd_table <- prd_table %>% rename("mod 2022 day2" = Estimate)
prd_table <- merge(prd_table, prd_mod4, by = "Variables1", all.x = T)
prd_table <- prd_table %>% rename("mod 2023" = Estimate)
prd_table <- merge(prd_table, prd_mod5, by = "Variables1", all.x = T)
prd_table <- prd_table %>% rename("mod 2023 day1" = Estimate)
prd_table <- merge(prd_table, prd_mod6, by = "Variables1", all.x = T)
prd_table <- prd_table %>% rename("mod 2023 day2" = Estimate)
prd_table <- merge(prd_table, prd_mod7, by = "Variables1", all.x = T)
prd_table <- prd_table %>% rename("mod 2022/2023" = Estimate)

prd_table[is.na(prd_table)] <- "-"
prd_table <- prd_table %>% arrange(Order)
prd_table <- prd_table %>% select(- "Order")
prd_table <- prd_table %>% rename(Variables = Variables1)

# Insert header rows:
# Fixed effects
Fixed_effects_row <- data.frame(matrix(c("Fixed effects", rep("\u03B2 (95% CI)", ncol(prd_table) - 1)), nrow = 1))
colnames(Fixed_effects_row) <- colnames(prd_table)
prd_table <- add_row(prd_table, Fixed_effects_row, .before = which(prd_table$Variable == "Intercept producing"))

# Random effects
Random_effects_row <- data.frame(matrix(c("Random effects", rep("\u03C3\u00B2 (95% CI)", ncol(prd_table) - 1)), nrow = 1))
colnames(Random_effects_row) <- colnames(prd_table)
prd_table <- add_row(prd_table, Random_effects_row, .before = which(prd_table$Variable == "Var DIE"))

# Covariance
Covariance_row <- data.frame(matrix(c("Covariance", rep("\u03C3 (95% CI)", ncol(prd_table) - 1)), nrow = 1))
colnames(Covariance_row) <- colnames(prd_table)
prd_table <- add_row(prd_table, Covariance_row, .before = which(prd_table$Variable == "Cov DIE-IIE"))

# Repeatability
Repeatability_row <- data.frame(matrix(c("Repeatability", rep("R (95% CI)", ncol(prd_table) - 1)), nrow = 1))
colnames(Repeatability_row) <- colnames(prd_table)
prd_table <- add_row(prd_table, Repeatability_row, .before = which(prd_table$Variable == "R DIE"))

# Correlations
Correlations_row <- data.frame(matrix(c("Correlations", rep("r (95% CI)", ncol(prd_table) - 1)), nrow = 1))
colnames(Correlations_row) <- colnames(prd_table)
prd_table <- add_row(prd_table, Correlations_row, .before = which(prd_table$Variable == "Cor DIE-IIE"))

# Remove identifying strings (e.g. Var, Cor)
prd_table$Variables <- str_remove(prd_table$Variable, pattern = str_c(c("Var ", "variance", "Cov ", "R ", "Cor "), collapse = "|"))

sjPlot::tab_df(prd_table, alternate.rows = T)

sjPlot::tab_df(prd_table, alternate.rows = T, file = paste("Plots/Table_prd_split.doc"))

#===============================================================================
# End of script #