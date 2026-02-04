
#############
#########
####
##
#

#----------- Archipelagos Models

#--- Past and present variables

# Marine - Region level; Terrestrial - Island level
#rm (list =ls())

library(glmmTMB)
library(ggplot2)
library(dplyr)
library(tidyr)
library(broom.mixed)
library(sjPlot)


if (!require("pacman")) install.packages("pacman")

pacman::p_load(tidyverse, BiocManager, rstan, brms, tidybayes, bayestestR, bayesplot, shinystan, marginaleffects, cowplot,
               tidymodels, xgboost, caret, doSNOW, modelr, paletteer, forcats, patchwork, parameters, ggeffects, ggridges, gganimate, gifski, transformr, glmmTMB, emmeans
)



#

#options(mc.cores = parallel::detectCores())

#

# Parallel 

total_cores <- 100L
chains <- 4L
threads_per_chain <- floor(total_cores / chains)  # = 25
message("using chains=", chains, " threads_per_chain=", threads_per_chain,
        " total threads=", chains * threads_per_chain)


#Beta diversity and variables data ----

#-- MARINE

# Corals taxonomic -----
model_corals_taxonomic <- read.csv("model_data_coral_taxonomic_arch.csv", header = TRUE) %>% 
  mutate_if(is.character, as.factor) %>% glimpse

# Corals functional -----
model_corals_functional <- read.csv("model_data_coral_functional_arch.csv", header = TRUE)%>% 
  mutate_if(is.character, as.factor) %>% glimpse

# Fish taxonomic -----
model_fish_taxonomic <- read.csv("model_data_fish_taxonomic_arch.csv", header = TRUE)%>% 
  mutate_if(is.character, as.factor) %>% glimpse

# Fish functional -----
model_fish_functional <- read.csv("model_data_fish_functional_arch_VERSION2.csv", header = TRUE)%>% 
  mutate_if(is.character, as.factor) %>% glimpse



#-- TERRESTRIAL

# Plants taxonomic -----
model_plants_taxonomic <- read.csv("plant_taxonomic_only_native_final.csv", header = TRUE) %>% 
  mutate_if(is.character, as.factor) %>% glimpse

# Plants functional -----
model_plants_functional <- read.csv("plant_functional_only_native_final.csv", header = TRUE) %>% 
  mutate_if(is.character, as.factor) %>% glimpse

# Birds taxonomic -----
model_birds_taxonomic <- read.csv("bird_taxonomic_only_native.csv", header = TRUE) %>% 
  mutate_if(is.character, as.factor) %>% glimpse

# Birds functional -----
model_birds_functional <- read.csv("bird_functional_only_native.csv", header = TRUE) %>% 
  mutate_if(is.character, as.factor) %>% glimpse



#-------------------------------------------------------------------------------


#1. Marine - CORALS BETA TAXONOMIC ----

## log for arch distance and arch quaternary isolation

# model_corals_taxonomic <- model_corals_taxonomic %>%
#   mutate(
#     log_present_isolation = log1p(diff_isolation),
#     log_past_isolation    = log1p(diff_past_isolation),
#     log_distance          = log1p(Distance)
#   )



## Scaling variables - difference between marine regions ----
corals_taxonomic <- model_corals_taxonomic %>%
  mutate(
    beta_sor_adj = beta_sor - 0.000001,  
    diff_isolation_sc = scale(diff_isolation)[,1], #quaternary distance # log
    diff_reef_area_sc = scale(diff_reef_area)[,1],
    diff_past_isolation_sc = scale(diff_past_isolation)[,1], #quaternary distance # log
    diff_past_reef_area_sc = scale(diff_past_reef_area)[,1],
    diff_age_sc = scale(diff_age)[,1],
    dist_sc = scale(Distance)[,1], # log
    diff_present_sst = scale(diff_sst)[,1],
    diff_past_sst = scale(diff_past_sst)
  )

corals_taxonomic %>% glimpse

#------------------------- Models

# Null

# model_beta_corals_taxonomic_present_null_gp3 <- brm(
#   beta_sor_adj ~ 1 + (1|Group1) + (1|Group2) + (1|key),
#   data = corals_taxonomic,
#   family = beta_family(),
#   chains = 4,
#   cores = 4,
#   iter = 4000,
#   warmup = 3000,
#   control = list(adapt_delta = 0.99, max_treedepth = 14),
#   save_pars = save_pars(all = T)
# )

# model_beta_corals_taxonomic_present_null_gp3 <- add_criterion(model_beta_corals_taxonomic_present_null_gp3, "loo", save_psis = TRUE, reloo = T)


# model_beta_corals_taxonomic_present_null <- brm(
#   beta_sor_adj ~ 1,
#   data = corals_taxonomic,
#   family = beta_family(),
#   chains = 4,
#   cores = 4,
#   iter = 4000,
#   warmup = 3000,
#   control = list(adapt_delta = 0.99, max_treedepth = 14),
#   save_pars = save_pars(all = T)
# )

# model_beta_corals_taxonomic_present_null <- add_criterion(model_beta_corals_taxonomic_present_null, "loo", save_psis = TRUE, reloo = T)
# 
# loo_compare(model_beta_corals_taxonomic_present_null, model_beta_corals_taxonomic_present_null_gp2, model_beta_corals_taxonomic_present_null_gp1, model_beta_corals_taxonomic_present_null_gp3)

model_beta_corals_taxonomic_present_null_gp2 <- brm(
  beta_sor_adj ~ 1 + 
    (1+dist_sc|Group1) + 
    (1+dist_sc|Group2),
  data = corals_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  control = list(adapt_delta = 0.99, max_treedepth = 14),
  save_pars = save_pars(all = T)
)

model_beta_corals_taxonomic_present_null_gp2 <- add_criterion(model_beta_corals_taxonomic_present_null_gp2, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_corals_taxonomic_present_null_gp2,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_corals_taxonomic_present_null_gp2.rds"
)


# Assess model validity using LOO
loo(model_beta_corals_taxonomic_present_null_gp2)

# Present
model_beta_corals_taxonomic_present_bayes_all <- brm(
  beta_sor_adj ~ diff_isolation_sc + diff_reef_area_sc + diff_age_sc + dist_sc + diff_present_sst +  
    (1+dist_sc|Group1) + (1+dist_sc|Group2),
  data = corals_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_corals_taxonomic_present_bayes_all <- add_criterion(model_beta_corals_taxonomic_present_bayes_all, "loo", save_psis = TRUE, reloo = T)

model_beta_corals_taxonomic_present_bayes_all

saveRDS(
  model_beta_corals_taxonomic_present_bayes_all,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_corals_taxonomic_present_bayes_all.rds"
)

# Assess model validity using LOO
loo(model_beta_corals_taxonomic_present_bayes_all)

#Look at model summary its on 'logit' scale, because of bernouli family link function
#model_beta_corals_taxonomic_present_bayes_2

#loo_compare(model_beta_corals_taxonomic_present_bayes_all, model_beta_corals_taxonomic_present_null_gp2)
## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_corals_taxonomic_present_bayes_all , ndraws = 50) 
pp_check(model_beta_corals_taxonomic_present_bayes_all, type = "ecdf_overlay", ndraws = 50) 

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_corals_taxonomic_present_bayes_all$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_corals_taxonomic_present_bayes_all, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_corals_taxonomic_present_bayes_all, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_corals_taxonomic_present_bayes_all, type = 'rhat') 


#back transform model summary with point estimates - Now Odds ratios (3:1 odds of bleaching)
model_parameters(model_beta_corals_taxonomic_present_bayes_all, exponentiate = T, ci = 0.95)



#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####

#Plotting isolation effects
# 
p1 <- corals_taxonomic %>% 
  modelr::data_grid(
    diff_isolation_sc = seq_range(diff_isolation_sc, n = 10),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_present_sst = mean(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_taxonomic_present_bayes_all,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_isolation_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in isolation", y = "Beta diversity")

#
p2 <- corals_taxonomic %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = seq_range(diff_reef_area_sc, n = 10),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_present_sst = mean(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_taxonomic_present_bayes_all,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_reef_area_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in reef area", y = "Beta diversity")

# 
p3 <- corals_taxonomic %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    diff_age_sc = seq_range(diff_age_sc, n = 10),
    dist_sc = mean(dist_sc),
    diff_present_sst = mean(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_taxonomic_present_bayes_all,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_age_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in age", y = "Beta diversity")

# 
p4 <- corals_taxonomic %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    dist_sc = seq_range(dist_sc, n = 10),
    diff_age_sc = mean(diff_age_sc),
    diff_present_sst = mean(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_taxonomic_present_bayes_all,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = dist_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch distance", y = "Beta diversity")

p5 <- corals_taxonomic %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    dist_sc = mean(dist_sc),
    diff_age_sc = mean(diff_age_sc),
    diff_present_sst = seq_range(diff_present_sst, n = 10)
  ) %>% 
  add_epred_draws(model_beta_corals_taxonomic_present_bayes_all,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_present_sst)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch SST", y = "Beta diversity")


library(cowplot)

p1 / p2 / p3 / p4 / p5






#2. Marine - CORALS BETA FUNCTIONAL ----


## log for arch distance and arch quaternary isolation

# model_corals_functional <- model_corals_functional %>%
#   mutate(
#     log_present_isolation = log1p(diff_isolation),
#     log_past_isolation    = log1p(diff_past_isolation),
#     log_distance          = log1p(Distance)
#   )


## Scaling variables - difference between marine regions ----
corals_functional <- model_corals_functional %>%
  mutate(
    beta_sor_adj = beta_sor - 0.000001,  
    diff_isolation_sc = scale(diff_isolation)[,1], #quaternary distance # log
    diff_reef_area_sc = scale(diff_reef_area)[,1],
    diff_past_isolation_sc = scale(diff_past_isolation)[,1], #quaternary distance # log
    diff_past_reef_area_sc = scale(diff_past_reef_area)[,1],
    diff_age_sc = scale(diff_age)[,1],
    dist_sc = scale(Distance)[,1], # log
    diff_present_sst = scale(diff_sst)[,1],
    diff_past_sst = scale(diff_past_sst)[,1]
  )

corals_functional %>% glimpse

#------------------------- Models


# Bayes models
# Null

model_beta_corals_functional_present_null_gp2 <- brm(
  beta_sor_adj ~ 1 + 
    (1+dist_sc|Group1) + (1+dist_sc|Group2),
  data = corals_functional,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  control = list(adapt_delta = 0.99, max_treedepth = 14),
  save_pars = save_pars(all = T)
)

model_beta_corals_functional_present_null_gp2 <- add_criterion(model_beta_corals_functional_present_null_gp2, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_corals_functional_present_null_gp2,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_corals_functional_present_null_gp2.rds"
)

# Assess model validity using LOO
loo(model_beta_corals_functional_present_null_gp2)

# Present
model_beta_corals_functional_present_bayes_all <- brm(
  beta_sor_adj ~ diff_isolation_sc + diff_reef_area_sc + diff_age_sc + dist_sc + diff_present_sst + 
   (1+dist_sc|Group1) + (1+dist_sc|Group2),
  data = corals_functional,
  family = beta_family(),
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_corals_functional_present_bayes_all <- add_criterion(model_beta_corals_functional_present_bayes_all, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_corals_functional_present_bayes_all,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_corals_functional_present_bayes_all.rds"
)

# Assess model validity using LOO
loo(model_beta_corals_functional_present_bayes_all)

#Look at model summary its on 'logit' scale, because of bernouli family link function
model_beta_corals_functional_present_bayes_all

#loo_compare(model_beta_corals_functional_present_bayes_all, model_beta_corals_functional_present_null_gp2)

## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_corals_functional_present_null_gp2 , ndraws = 50) 
pp_check(model_beta_corals_functional_present_bayes_all, type = "ecdf_overlay", ndraws = 50) 

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_corals_functional_present_bayes_all$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_corals_functional_present_bayes_all, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_corals_functional_present_bayes_all, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_corals_functional_present_bayes_all, type = 'rhat') 


#back transform model summary with point estimates - Now Odds ratios (3:1 odds of bleaching)
model_parameters(model_beta_corals_functional_present_bayes_all, exponentiate = T, ci = 0.95)


#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####

#Plotting isolation effects
p1 <- corals_functional %>% 
  modelr::data_grid(
    diff_isolation_sc = seq_range(diff_isolation_sc, n = 10),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_present_sst = mean(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_functional_present_bayes_all,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_isolation_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in isolation", y = "Beta diversity")

#
p2 <- corals_functional %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = seq_range(diff_reef_area_sc, n = 10),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_present_sst = mean(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_functional_present_bayes_all,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_reef_area_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in reef area", y = "Beta diversity")

# 
p3 <- corals_functional %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    diff_age_sc = seq_range(diff_age_sc, n = 10),
    dist_sc = mean(dist_sc),
    diff_present_sst = mean(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_functional_present_bayes_all,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_age_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in age", y = "Beta diversity")

# 
p4 <- corals_functional %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    dist_sc = seq_range(dist_sc, n = 10),
    diff_age_sc = mean(diff_age_sc),
    diff_present_sst = mean(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_functional_present_bayes_all,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = dist_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch distance", y = "Beta diversity")


p5 <- corals_functional %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    dist_sc = mean(dist_sc),
    diff_age_sc = mean(diff_age_sc),
    diff_present_sst = seq_range(diff_present_sst, n = 10)
  ) %>% 
  add_epred_draws(model_beta_corals_functional_present_bayes_all,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_present_sst)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch SST", y = "Beta diversity")


library(cowplot)

p1 / p2 / p3 / p4 / p5





#3. Marine - FISH BETA TAXONOMIC ----


## Scaling variables - difference between marine regions ----
## log for arch distance and arch quaternary isolation

# model_fish_taxonomic <- model_fish_taxonomic %>%
#   mutate(
#     log_present_isolation = log1p(diff_isolation),
#     log_past_isolation    = log1p(diff_past_isolation),
#     log_distance          = log1p(Distance)
#   )


## Scaling variables - difference between marine regions ----
fish_taxonomic <- model_fish_taxonomic %>%
  mutate(
    beta_sor_adj = beta_sor - 0.000001,  
    diff_isolation_sc = scale(diff_isolation)[,1], #quaternary distance # log
    diff_reef_area_sc = scale(diff_reef_area)[,1],
    diff_past_isolation_sc = scale(diff_past_isolation)[,1], #quaternary distance # log
    diff_past_reef_area_sc = scale(diff_past_reef_area)[,1],
    diff_age_sc = scale(diff_age)[,1],
    dist_sc = scale(Distance)[,1], # log
    diff_present_sst = scale(diff_sst)[,1],
    diff_past_sst = scale(diff_past_sst)
  )

fish_taxonomic %>% glimpse



#------------------------- Models

# Null

# model_beta_fish_taxonomic_present_null_gp3 <- add_criterion(model_beta_fish_taxonomic_present_null_gp3, "loo", save_psis = TRUE, reloo = T)
# model_beta_fish_taxonomic_present_null_gp2 <- brm(
#   beta_sor_adj ~ 1 + (1|Group1) + (1|Group2) ,
#   data = fish_taxonomic,
#   family = beta_family(),
#   chains = 4,
#   cores = 4,
#   iter = 4000,
#   warmup = 3000,
#   control = list(adapt_delta = 0.99, max_treedepth = 14),
#   save_pars = save_pars(all = T)
# )

model_beta_fish_taxonomic_present_null_gp2 <- brm(
  beta_sor_adj ~ 1 + (1+dist_sc|Group1) + (1+dist_sc|Group2),
  data = fish_taxonomic,
  family = beta_family(),
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 3000,
  control = list(adapt_delta = 0.99, max_treedepth = 14),
  save_pars = save_pars(all = T)
)

model_beta_fish_taxonomic_present_null_gp2 <- add_criterion(model_beta_fish_taxonomic_present_null_gp2, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_fish_taxonomic_present_null_gp2,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_fish_taxonomic_present_null_gp2.rds"
)

# Assess model validity using LOO
loo(model_beta_fish_taxonomic_present_null_gp2)

#Look at model summary its on 'logit' scale, because of bernouli family link function
model_beta_fish_taxonomic_present_null_gp2


# Present
model_beta_fish_taxonomic_present_bayes_all <- brm(
  beta_sor_adj ~ diff_isolation_sc + diff_reef_area_sc + diff_age_sc + dist_sc + diff_present_sst + 
  (1+dist_sc|Group1) + (1+dist_sc|Group2),
  data = fish_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_fish_taxonomic_present_bayes_all <- add_criterion(model_beta_fish_taxonomic_present_bayes_all, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_fish_taxonomic_present_bayes_all,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_fish_taxonomic_present_bayes_all.rds"
)

model_beta_fish_taxonomic_present_bayes_all

# Assess model validity using LOO
loo(model_beta_fish_taxonomic_present_bayes_all)

#loo_compare(model_beta_fish_taxonomic_present_bayes_all, model_beta_fish_taxonomic_present_null_gp2)
## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_fish_taxonomic_present_bayes_all , ndraws = 50) 
pp_check(model_beta_fish_taxonomic_present_bayes_all, type = "ecdf_overlay", ndraws = 50) 

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_fish_taxonomic_present_bayes_all$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_fish_taxonomic_present_bayes_all, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_fish_taxonomic_present_bayes_all, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_fish_taxonomic_present_bayes_all, type = 'rhat') 


#back transform model summary with point estimates - Now Odds ratios (3:1 odds of bleaching)
model_parameters(model_beta_fish_taxonomic_present_bayes_all, exponentiate = T, ci = 0.95)


#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####

#Plotting isolation effects
p1 <- fish_taxonomic %>% 
  modelr::data_grid(
    diff_isolation_sc = seq_range(diff_isolation_sc, n = 10),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_present_sst = mean(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_taxonomic_present_bayes_all,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_isolation_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in isolation", y = "Beta diversity")

#
p2 <- fish_taxonomic %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = seq_range(diff_reef_area_sc, n = 10),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_present_sst = mean(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_taxonomic_present_bayes_all,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_reef_area_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in reef area", y = "Beta diversity")

# 
p3 <- fish_taxonomic %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    diff_age_sc = seq_range(diff_age_sc, n = 10),
    dist_sc = mean(dist_sc),
    diff_present_sst = mean(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_taxonomic_present_bayes_all,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_age_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in age", y = "Beta diversity")

# 
p4 <- fish_taxonomic %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    dist_sc = seq_range(dist_sc, n = 10),
    diff_age_sc = mean(diff_age_sc),
    diff_present_sst = mean(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_taxonomic_present_bayes_all,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = dist_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch distance", y = "Beta diversity")

p5 <- fish_taxonomic %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    dist_sc = mean(dist_sc),
    diff_age_sc = mean(diff_age_sc),
    diff_present_sst = seq_range(diff_present_sst, n = 10)
  ) %>% 
  add_epred_draws(model_beta_fish_taxonomic_present_bayes_all,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_present_sst)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch SST", y = "Beta diversity")


library(cowplot)

p1 / p2 / p3 / p4 / p5






#4. Marine - FISH BETA FUNCTIONAL ----


# ## Scaling variables - difference between marine regions ----
# model_fish_functional <- model_fish_functional %>%
#   mutate(
#     log_present_isolation = log1p(diff_isolation),
#     log_past_isolation    = log1p(diff_past_isolation),
#     log_distance          = log1p(Distance)
#   )


## Scaling variables - difference between marine regions ----
fish_functional <- model_fish_functional %>%
  mutate(
    beta_sor_adj = beta_sor - 0.000001,  
    diff_isolation_sc = scale(diff_isolation)[,1], #quaternary distance # log
    diff_reef_area_sc = scale(diff_reef_area)[,1],
    diff_past_isolation_sc = scale(diff_past_isolation)[,1], #quaternary distance # log
    diff_past_reef_area_sc = scale(diff_past_reef_area)[,1],
    diff_age_sc = scale(diff_age)[,1],
    dist_sc = scale(Distance)[,1], # log
    diff_present_sst = scale(diff_sst)[,1],
    diff_past_sst = scale(diff_past_sst)[,1]
  )

fish_functional %>% glimpse

#------------------------- Models


# Bayes models
# Null
model_beta_fish_functional_present_null_gp2 <- brm(
  beta_sor_adj ~ 1 + (1+dist_sc|Group1) + (1+dist_sc|Group2) ,
  data = fish_functional,
  family = beta_family(),
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 3000,
  control = list(adapt_delta = 0.99, max_treedepth = 14),
  save_pars = save_pars(all = T)
)

model_beta_fish_functional_present_null_gp2 <- add_criterion(model_beta_fish_functional_present_null_gp2, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_fish_functional_present_null_gp2,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_fish_functional_present_null_gp2.rds"
)

# Assess model validity using LOO
loo(model_beta_fish_functional_present_null_gp2)


# Present
model_beta_fish_functional_present_bayes_all <- brm(
  beta_sor_adj ~  1 + diff_isolation_sc + diff_reef_area_sc + diff_age_sc + dist_sc + diff_present_sst +
  (1+dist_sc|Group1) + (1+dist_sc|Group2),
  data = fish_functional,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)



# Add Leave One Out cross validation criterion
model_beta_fish_functional_present_bayes_all <- add_criterion(model_beta_fish_functional_present_bayes_all, "loo", save_psis = TRUE, reloo = T, moment_match = TRUE)

saveRDS(
  model_beta_fish_functional_present_bayes_all,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_fish_functional_present_bayes_all.rds"
)


# Assess model validity using LOO
loo(model_beta_fish_functional_present_bayes_all)

#loo_compare(model_beta_fish_functional_present_bayes_all, model_beta_fish_functional_present_null_gp2)

## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_fish_functional_present_bayes_all, ndraws = 50) 
pp_check(model_beta_fish_functional_present_bayes_all, type = "ecdf_overlay", ndraws = 50) 
pp_check(model_beta_fish_functional_present_bayes_all, ndraws = 50, type = "dens_overlay_grouped", group = "Group1")
pp_check(model_beta_fish_functional_present_bayes_all, ndraws = 50, type = "scatter_avg_grouped", group = "Group1")
pp_check(model_beta_fish_functional_present_bayes_all, ndraws = 50, type = "error_scatter_avg_grouped", group = "Group1")
pp_check(model_beta_fish_functional_present_bayes_all, type='pit_ecdf_grouped', group = "Group1", prob = 0.95, plot_diff = F)
pp_check(model_beta_fish_functional_present_bayes_all, type='pit_ecdf_grouped', group = "Group2", prob = 0.95, plot_diff = F)

#pp_check(model_beta_fish_functional_present_bayes_all, ndraws = 50, type = "scatter_avg_grouped", group = "Group2")

#pairs(model_beta_fish_functional_present_bayes_all, variable = variables(model_beta_fish_functional_present_bayes_all_g)[1:5])
#fish_functional %>% select(Group1, Group2) %>% unique %>% view

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_fish_functional_present_bayes_all$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_fish_functional_present_bayes_all, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_fish_functional_present_bayes_all, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_fish_functional_present_bayes_all, type = 'rhat') 


#back transform model summary with point estimates - Now Odds ratios (3:1 odds of bleaching)
model_parameters(model_beta_fish_functional_present_bayes_all, exponentiate = T, ci = 0.95)



#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####

#Plotting isolation effects
p1 <- fish_functional %>% 
  modelr::data_grid(
    diff_isolation_sc = seq_range(diff_isolation_sc, n = 10),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_present_sst = mean(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_functional_present_bayes_all,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_isolation_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in isolation", y = "Beta diversity")

#
p2 <- fish_functional %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = seq_range(diff_reef_area_sc, n = 10),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_present_sst = mean(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_functional_present_bayes_all,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_reef_area_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in reef area", y = "Beta diversity")

# 
p3 <- fish_functional %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    diff_age_sc = seq_range(diff_age_sc, n = 10),
    dist_sc = mean(dist_sc),
    diff_present_sst = mean(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_functional_present_bayes_all,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_age_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in age", y = "Beta diversity")

# 
p4 <- fish_functional %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    dist_sc = seq_range(dist_sc, n = 10),
    diff_age_sc = mean(diff_age_sc),
    diff_present_sst = mean(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_functional_present_bayes_all,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = dist_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch distance", y = "Beta diversity")

p5 <- fish_functional %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    dist_sc = mean(dist_sc),
    diff_age_sc = mean(diff_age_sc),
    diff_present_sst = seq_range(diff_present_sst, n = 10)
  ) %>% 
  add_epred_draws(model_beta_fish_functional_present_bayes_all,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_present_sst)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch SST", y = "Beta diversity")


library(cowplot)

p1 / p2 / p3 / p4 / p5






#5. Terrestrial - PLANTS BETA TAXONOMIC ----


model_plants_taxonomic <- model_plants_taxonomic %>%
  filter(!Island1 %in% "Farquhar_Atoll",
         !Island2 %in% "Farquhar_Atoll")


# Scaling variables - difference between archipelagos islands
plants_taxonomic <- model_plants_taxonomic %>%
  mutate(
    beta_sor_adj = beta_sor - 0.000001,
    beta_sor_adj = ifelse(beta_sor_adj <= 0, 1e-6, beta_sor_adj),
    d_past_sc = scale(d_past.y), #emerse distance past
    diff_emerse_past_sc = scale(diff_emerse_past), #emerse island area present
    diff_emerse_pres_sc = scale(diff_emerse_pres), #emerse island area past
    diff_island_age_sc = scale(diff_island_age),
    distance_km_sc = scale(Distance_km), #emerse distance present
    diff_present_air = scale(diif_present_air),
    diff_past_air = scale(diff_past_air)
  )

# scaled_diff_emerse_pres <- plants_taxonomic$diff_emerse_pres - mean(plants_taxonomic$diff_emerse_pres)/sd(plants_taxonomic$diff_emerse_pres)
# 
# diff_emerse_pres <- scaled_diff_emerse_pres * sd(plants_taxonomic$diff_emerse_pres) + mean(plants_taxonomic$diff_emerse_pres)
# 
# diff_emerse_pres_sc <-  scale(plants_taxonomic$diff_island_age_sc)
# diff_emerse_pres_sc %>% glimpse


#------------------------- Models

# Null

model_beta_plants_taxonomic_present_null <- brm(
  beta_sor_adj ~ 1 +
    (1 | Archipelago),
  data = plants_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 5000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


model_beta_plants_taxonomic_present_null <- add_criterion(model_beta_plants_taxonomic_present_null, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_plants_taxonomic_present_null,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_plants_taxonomic_present_null.rds"
)


# Assess model validity using LOO
loo(model_beta_plants_taxonomic_present_null)



# Present

model_beta_plants_taxonomic_present <- brm(
  beta_sor_adj ~  diff_emerse_pres_sc + diff_island_age_sc + distance_km_sc + diff_present_air +
    (1 | Archipelago),
  data = plants_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_plants_taxonomic_present <- add_criterion(model_beta_plants_taxonomic_present, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_plants_taxonomic_present,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_plants_taxonomic_present.rds"
)


loo(model_beta_plants_taxonomic_present)

#loo_compare(model_beta_plants_taxonomic_present_null, model_beta_plants_taxonomic_present)

#Look at model summary its on 'logit' scale, because of Beta family link function
model_beta_plants_taxonomic_present

## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_plants_taxonomic_present, ndraws = 50) 
pp_check(model_beta_plants_taxonomic_present, type = "ecdf_overlay", ndraws = 50) 

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_plants_taxonomic_present$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_plants_taxonomic_present, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_plants_taxonomic_present, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_plants_taxonomic_present, type = 'rhat') 


#back transform model summary with point estimates (Ratios)
model_parameters(model_beta_plants_taxonomic_present, exponentiate = T, ci = 0.95)


#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####


#Plotting isolation effects
p1 <- plants_taxonomic %>% 
  modelr::data_grid(distance_km_sc = seq_range(distance_km_sc, n = 100) %>% as.numeric,
                    diff_emerse_pres_sc = median(diff_emerse_pres_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_present_air = median(diff_present_air, na.rm = T))%>% 
  add_epred_draws(model_beta_plants_taxonomic_present, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = distance_km_sc)) +
  # geom_point(data = corals_taxonomic$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Distance between islands", y = "Beta diversity")

p2 <- plants_taxonomic %>% 
  modelr::data_grid(distance_km_sc = median(distance_km_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_pres_sc = seq_range(diff_emerse_pres_sc, n = 10),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_present_air = median(diff_present_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_plants_taxonomic_present, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_emerse_pres_sc)) +
  # geom_point(data = corals_taxonomic$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Emersed Island Area", y = "Beta diversity")

p3 <- plants_taxonomic %>% 
  modelr::data_grid(distance_km_sc = median(distance_km_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_pres_sc = median(diff_emerse_pres_sc, na.rm = T),
                    diff_island_age_sc = seq_range(diff_island_age_sc, n = 10),
                    diff_present_air = median(diff_present_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_plants_taxonomic_present, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_island_age_sc)) +
  # geom_point(data = corals_taxonomic$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Island Age", y = "Beta diversity")


p4 <- plants_taxonomic %>% 
  modelr::data_grid(distance_km_sc = median(distance_km_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_pres_sc = median(diff_emerse_pres_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_present_air = seq_range(diff_present_air, n = 10)) %>% 
  add_epred_draws(model_beta_plants_taxonomic_present, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_present_air)) +
  # geom_point(data = corals_taxonomic$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Air Temperature", y = "Beta diversity")

library(cowplot)

p1 / p2 / p3 / p4






#6. Terrestrial - PLANTS BETA FUNCTIONAL ----


model_plants_functional <- model_plants_functional %>%
  filter(!Island1 %in% "Farquhar_Atoll",
         !Island2 %in% "Farquhar_Atoll")



# Scaling variables - difference between archipelagos islands
plants_functional <- model_plants_functional %>%
  mutate(
    beta_sor_adj = beta_sor - 0.000001, 
    beta_sor_adj = ifelse(beta_sor_adj <= 0, 1e-6, beta_sor_adj),
    d_past_sc = scale(d_past.y), #emerse distance past
    diff_emerse_past_sc = scale(diff_emerse_past), #emerse island area present
    diff_emerse_pres_sc = scale(diff_emerse_pres), #emerse island area past
    diff_island_age_sc = scale(diff_island_age),
    distance_km_sc = scale(Distance_km), #emerse distance present
    diff_present_air = scale(diif_present_air),
    diff_past_air = scale(diff_past_air)
  )


#------------------------- Models

# Null

model_beta_plants_functional_present_null <- brm(
  beta_sor_adj ~ 1 +
    (1 | Archipelago),
  data = plants_functional,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


model_beta_plants_functional_present_null<- add_criterion(model_beta_plants_functional_present_null, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_plants_functional_present_null,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_plants_functional_present_null.rds"
)

# Assess model validity using LOO
loo(model_beta_plants_functional_present_null)



# Present

model_beta_plants_functional_present <- brm(
  beta_sor_adj ~  diff_emerse_pres_sc + diff_island_age_sc + distance_km_sc + diff_present_air +
    (1 | Archipelago),
  data = plants_functional,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_plants_functional_present <- add_criterion(model_beta_plants_functional_present, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_plants_functional_present,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_plants_functional_present.rds"
)

loo(model_beta_plants_functional_present)

#loo_compare(model_beta_plants_functional_present_null, model_beta_plants_functional_present)

#Look at model summary its on 'logit' scale, because of Beta family link function
model_beta_plants_functional_present

## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_plants_functional_present, ndraws = 50) 
pp_check(model_beta_plants_functional_present, type = "ecdf_overlay", ndraws = 50) 

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_plants_functional_present$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_plants_functional_present, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_plants_functional_present, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_plants_functional_present, type = 'rhat') 


#back transform model summary with point estimates (Ratios)
model_parameters(model_beta_plants_functional_present, exponentiate = T, ci = 0.95)


#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####


#Plotting isolation effects
p1 <- plants_functional %>% 
  modelr::data_grid(distance_km_sc = seq_range(distance_km_sc, n = 100) %>% as.numeric,
                    diff_emerse_pres_sc = median(diff_emerse_pres_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_present_air = median(diff_present_air, na.rm = T))%>% 
  add_epred_draws(model_beta_plants_functional_present, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = distance_km_sc)) +
  # geom_point(data = corals_functional$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Distance between islands", y = "Beta diversity")

p2 <- plants_functional %>% 
  modelr::data_grid(distance_km_sc = median(distance_km_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_pres_sc = seq_range(diff_emerse_pres_sc, n = 10),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_present_air = median(diff_present_air, na.rm = T))%>%  
  add_epred_draws(model_beta_plants_functional_present, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_emerse_pres_sc)) +
  # geom_point(data = corals_functional$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Emersed Island Area", y = "Beta diversity")

p3 <- plants_functional %>% 
  modelr::data_grid(distance_km_sc = median(distance_km_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_pres_sc = median(diff_emerse_pres_sc, na.rm = T),
                    diff_island_age_sc = seq_range(diff_island_age_sc, n = 10),
                    diff_present_air = median(diff_present_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_plants_functional_present, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_island_age_sc)) +
  # geom_point(data = corals_functional$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Island Age", y = "Beta diversity")


p4 <- plants_taxonomic %>% 
  modelr::data_grid(distance_km_sc = median(distance_km_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_pres_sc = median(diff_emerse_pres_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_present_air = seq_range(diff_present_air, n = 10)) %>% 
  add_epred_draws(model_beta_plants_taxonomic_present, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_present_air)) +
  # geom_point(data = corals_taxonomic$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Air Temperature", y = "Beta diversity")

p4 <- plants_taxonomic %>% 
  modelr::data_grid(distance_km_sc = median(distance_km_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_pres_sc = median(diff_emerse_pres_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_present_air = seq_range(diff_present_air, n = 10)) %>% 
  add_epred_draws(model_beta_plants_taxonomic_present, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_present_air)) +
  # geom_point(data = corals_taxonomic$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Air Temperature", y = "Beta diversity")


library(cowplot)

p1 / p2 / p3 / p4 




#7. Terrestrial - BIRDS BETA TAXONOMIC ----


model_birds_taxonomic <- model_birds_taxonomic %>%
  filter(!Island1 %in% "Farquhar_Atoll",
         !Island2 %in% "Farquhar_Atoll")



# Scaling variables - difference between archipelagos islands
birds_taxonomic <- model_birds_taxonomic %>%
  mutate(
    beta_sor_adj = beta_sor - 0.000001,
    beta_sor_adj = ifelse(beta_sor_adj <= 0, 1e-6, beta_sor_adj),
    d_past_sc = scale(d_past.y), #emerse distance past
    diff_emerse_past_sc = scale(diff_emerse_past), #emerse island area present
    diff_emerse_pres_sc = scale(diff_emerse_pres), #emerse island area past
    diff_island_age_sc = scale(diff_island_age),
    distance_km_sc = scale(Distance_km), #emerse distance present
    diff_present_air = scale(diif_present_air),
    diff_past_air = scale(diff_past_air)
  )




#------------------------- Models

# Null

model_beta_birds_taxonomic_present_null <- brm(
  beta_sor_adj ~ 1 +
    (1 | Archipelago),
  data = birds_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 5000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


model_beta_birds_taxonomic_present_null <- add_criterion(model_beta_birds_taxonomic_present_null, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_birds_taxonomic_present_null,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_birds_taxonomic_present_null.rds"
)


# Assess model validity using LOO
loo(model_beta_birds_taxonomic_present_null)



# Present

model_beta_birds_taxonomic_present <- brm(
  beta_sor_adj ~  diff_emerse_pres_sc + diff_island_age_sc + distance_km_sc + diff_present_air +
    (1 | Archipelago),
  data = birds_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_birds_taxonomic_present <- add_criterion(model_beta_birds_taxonomic_present, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_birds_taxonomic_present,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_birds_taxonomic_present.rds"
)

loo(model_beta_birds_taxonomic_present)

#loo_compare(model_beta_birds_taxonomic_present_null, model_beta_birds_taxonomic_present)

#Look at model summary its on 'logit' scale, because of Beta family link function
model_beta_birds_taxonomic_present

## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_birds_taxonomic_present, ndraws = 50) 
pp_check(model_beta_birds_taxonomic_present, type = "ecdf_overlay", ndraws = 50) 

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_birds_taxonomic_present$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_birds_taxonomic_present, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_birds_taxonomic_present, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_birds_taxonomic_present, type = 'rhat') 


#back transform model summary with point estimates (Ratios)
model_parameters(model_beta_birds_taxonomic_present, exponentiate = T, ci = 0.95)


#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####


#Plotting isolation effects
p1 <- birds_taxonomic %>% 
  modelr::data_grid(distance_km_sc = seq_range(distance_km_sc, n = 100) %>% as.numeric,
                    diff_emerse_pres_sc = median(diff_emerse_pres_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_present_air = median(diff_present_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_birds_taxonomic_present, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = distance_km_sc)) +
  # geom_point(data = corals_taxonomic$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Distance between islands", y = "Beta diversity")

p2 <- birds_taxonomic %>% 
  modelr::data_grid(distance_km_sc = median(distance_km_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_pres_sc = seq_range(diff_emerse_pres_sc, n = 10),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_present_air = median(diff_present_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_birds_taxonomic_present, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_emerse_pres_sc)) +
  # geom_point(data = corals_taxonomic$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Emersed Island Area", y = "Beta diversity")

p3 <- birds_taxonomic %>% 
  modelr::data_grid(distance_km_sc = median(distance_km_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_pres_sc = median(diff_emerse_pres_sc, na.rm = T),
                    diff_island_age_sc = seq_range(diff_island_age_sc, n = 10),
                    diff_present_air = median(diff_present_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_birds_taxonomic_present, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_island_age_sc)) +
  # geom_point(data = corals_taxonomic$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Island Age", y = "Beta diversity")


p4 <- birds_taxonomic %>% 
  modelr::data_grid(distance_km_sc = median(distance_km_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_pres_sc = median(diff_emerse_pres_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_present_air = seq_range(diff_present_air, n = 10)) %>% 
  add_epred_draws(model_beta_birds_taxonomic_present, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_present_air)) +
  # geom_point(data = corals_taxonomic$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Air Temperature", y = "Beta diversity")

library(cowplot)

p1 / p2 / p3 / p4





#8. Terrestrial - BIRDS BETA FUNCTIONAL ----


model_birds_functional <- model_birds_functional %>%
  filter(!Island1 %in% "Farquhar_Atoll",
         !Island2 %in% "Farquhar_Atoll")



# Scaling variables - difference between archipelagos islands
birds_functional <- model_birds_functional %>%
  mutate(
    beta_sor_adj = beta_sor - 0.000001,
    beta_sor_adj = ifelse(beta_sor_adj <= 0, 1e-6, beta_sor_adj),
    d_past_sc = scale(d_past.y), #emerse distance past
    diff_emerse_past_sc = scale(diff_emerse_past), #emerse island area present
    diff_emerse_pres_sc = scale(diff_emerse_pres), #emerse island area past
    diff_island_age_sc = scale(diff_island_age),
    distance_km_sc = scale(Distance_km), #emerse distance present
    diff_present_air = scale(diif_present_air),
    diff_past_air = scale(diff_past_air)
  )




#------------------------- Models

# Null

model_beta_birds_functional_present_null <- brm(
  beta_sor_adj ~ 1 +
    (1 | Archipelago),
  data = birds_functional,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 5000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.98, max_treedepth = 14),
)


model_beta_birds_functional_present_null <- add_criterion(model_beta_birds_functional_present_null, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_birds_functional_present_null,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_birds_functional_present_null.rds"
)


# Assess model validity using LOO
loo(model_beta_birds_functional_present_null)



# Present

model_beta_birds_functional_present <- brm(
  beta_sor_adj ~  diff_emerse_pres_sc + diff_island_age_sc + distance_km_sc + diff_present_air +
    (1 | Archipelago),
  data = birds_functional,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_birds_functional_present <- add_criterion(model_beta_birds_functional_present, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_birds_functional_present,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_birds_functional_present.rds"
)

loo(model_beta_birds_functional_present)


#loo_compare(model_beta_birds_functional_present_null, model_beta_birds_functional_present)

#Look at model summary its on 'logit' scale, because of Beta family link function
model_beta_birds_functional_present

## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_birds_functional_present, ndraws = 50) 
pp_check(model_beta_birds_functional_present, type = "ecdf_overlay", ndraws = 50) 

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_birds_functional_present$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_birds_functional_present, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_birds_functional_present, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_birds_functional_present, type = 'rhat') 


#back transform model summary with point estimates (Ratios)
model_parameters(model_beta_birds_functional_present, exponentiate = T, ci = 0.95)


#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####


#Plotting isolation effects
p1 <- birds_functional %>% 
  modelr::data_grid(distance_km_sc = seq_range(distance_km_sc, n = 100) %>% as.numeric,
                    diff_emerse_pres_sc = median(diff_emerse_pres_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_present_air = median(diff_present_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_birds_functional_present, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = distance_km_sc)) +
  # geom_point(data = corals_functional$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Distance between islands", y = "Beta diversity")

p2 <- birds_functional %>% 
  modelr::data_grid(distance_km_sc = median(distance_km_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_pres_sc = seq_range(diff_emerse_pres_sc, n = 10),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_present_air = median(diff_present_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_birds_functional_present, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_emerse_pres_sc)) +
  # geom_point(data = corals_functional$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Emersed Island Area", y = "Beta diversity")

p3 <- birds_functional %>% 
  modelr::data_grid(distance_km_sc = median(distance_km_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_pres_sc = median(diff_emerse_pres_sc, na.rm = T),
                    diff_island_age_sc = seq_range(diff_island_age_sc, n = 10),
                    diff_present_air = median(diff_present_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_birds_functional_present, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_island_age_sc)) +
  # geom_point(data = corals_functional$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Island Age", y = "Beta diversity")


p4 <- birds_functional %>% 
  modelr::data_grid(distance_km_sc = median(distance_km_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_pres_sc = median(diff_emerse_pres_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_present_air = seq_range(diff_present_air, n = 10)) %>% 
  add_epred_draws(model_beta_birds_functional_present, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_present_air)) +
  # geom_point(data = corals_taxonomic$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Air Temperature", y = "Beta diversity")

library(cowplot)

p1 / p2 / p3 / p4


# PAST ----

#-------------------------------------------------------------------------------


#1. Marine - CORALS BETA TAXONOMIC ----


## Scaling variables - difference between marine regions ----
corals_taxonomic <- model_corals_taxonomic %>%
  mutate(
    beta_sor_adj = beta_sor - 0.000001,  
    diff_isolation_sc = scale(diff_isolation)[,1], #quaternary distance # log
    diff_reef_area_sc = scale(diff_reef_area)[,1],
    diff_past_isolation_sc = scale(diff_past_isolation)[,1], #quaternary distance # log
    diff_past_reef_area_sc = scale(diff_past_reef_area)[,1],
    diff_age_sc = scale(diff_age)[,1],
    dist_sc = scale(Distance)[,1], # log
    diff_present_sst = scale(diff_sst)[,1],
    diff_past_sst = scale(diff_past_sst)[,1]
  )

corals_taxonomic %>% glimpse

#------------------------- Models

# past


# Bayes models
# Null

# model_beta_corals_taxonomic_past_null_gp3 <- brm(
#   beta_sor_adj ~ 1 + (1|Group1) + (1|Group2) + (1|key),
#   data = corals_taxonomic,
#   family = beta_family(),
#   chains = 4,
#   cores = 4,
#   iter = 4000,
#   warmup = 3000,
#   control = list(adapt_delta = 0.99, max_treedepth = 14),
#   save_pars = save_pars(all = T)
# )

# model_beta_corals_taxonomic_past_null_gp3 <- add_criterion(model_beta_corals_taxonomic_past_null_gp3, "loo", save_psis = TRUE, reloo = T)


# model_beta_corals_taxonomic_past_null <- brm(
#   beta_sor_adj ~ 1,
#   data = corals_taxonomic,
#   family = beta_family(),
#   chains = 4,
#   cores = 4,
#   iter = 4000,
#   warmup = 3000,
#   control = list(adapt_delta = 0.99, max_treedepth = 14),
#   save_pars = save_pars(all = T)
# )

# model_beta_corals_taxonomic_past_null <- add_criterion(model_beta_corals_taxonomic_past_null, "loo", save_psis = TRUE, reloo = T)
# 
# loo_compare(model_beta_corals_taxonomic_past_null, model_beta_corals_taxonomic_past_null_gp2, model_beta_corals_taxonomic_past_null_gp1, model_beta_corals_taxonomic_past_null_gp3)

model_beta_corals_taxonomic_past_null_gp2 <- brm(
  beta_sor_adj ~ 1 + 
    (1+dist_sc|Group1) + (1+dist_sc|Group2),
  data = corals_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  control = list(adapt_delta = 0.99, max_treedepth = 14),
  save_pars = save_pars(all = T)
)

model_beta_corals_taxonomic_past_null_gp2 <- add_criterion(model_beta_corals_taxonomic_past_null_gp2, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_corals_taxonomic_past_null_gp2,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_corals_taxonomic_past_null_gp2.rds"
)

# Assess model validity using LOO
loo(model_beta_corals_taxonomic_past_null_gp2)

# past
model_beta_corals_taxonomic_past_bayes_all <- brm(
  beta_sor_adj ~ diff_past_isolation_sc + diff_past_reef_area_sc + diff_age_sc + dist_sc + diff_past_sst +
    (1+dist_sc|Group1) + (1+dist_sc|Group2),
  data = corals_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_corals_taxonomic_past_bayes_all <- add_criterion(model_beta_corals_taxonomic_past_bayes_all, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_corals_taxonomic_past_bayes_all,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_corals_taxonomic_past_bayes_all.rds"
)

#Look at model summary its on 'logit' scale, because of bernouli family link function
model_beta_corals_taxonomic_past_bayes_all

# Assess model validity using LOO
loo(model_beta_corals_taxonomic_past_bayes_all)


#loo_compare(model_beta_corals_taxonomic_past_bayes_all, model_beta_corals_taxonomic_past_null_gp2)
## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_corals_taxonomic_past_bayes_all , ndraws = 50) 
pp_check(model_beta_corals_taxonomic_past_bayes_all, type = "ecdf_overlay", ndraws = 50) 

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_corals_taxonomic_past_bayes_all$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_corals_taxonomic_past_bayes_all, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_corals_taxonomic_past_bayes_all, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_corals_taxonomic_past_bayes_all, type = 'rhat') 


#back transform model summary with point estimates - Now Odds ratios (3:1 odds of bleaching)
model_parameters(model_beta_corals_taxonomic_past_bayes_all, exponentiate = T, ci = 0.95)



#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####

#Plotting isolation effects
# 
p1 <- corals_taxonomic %>% 
  modelr::data_grid(
    diff_past_isolation_sc = seq_range(diff_isolation_sc, n = 10),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_past_sst = mean(diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_taxonomic_past_bayes_all,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_past_isolation_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in isolation", y = "Beta diversity")

#
p2 <- corals_taxonomic %>% 
  modelr::data_grid(
    diff_past_isolation_sc = seq_range(diff_past_isolation_sc, n=10),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_past_sst = mean(diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_taxonomic_past_bayes_all,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_past_isolation_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in area", y = "Beta diversity")

# 
p3 <- corals_taxonomic %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_past_isolation_sc),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    diff_age_sc = seq_range(diff_age_sc, n = 10),
    dist_sc = mean(dist_sc),
    diff_past_sst = mean(diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_taxonomic_past_bayes_all,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_age_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in age", y = "Beta diversity")

# 
p4 <- corals_taxonomic %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_past_isolation_sc),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    dist_sc = seq_range(dist_sc, n = 10),
    diff_age_sc = mean(diff_age_sc),
    diff_past_sst = mean(diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_taxonomic_past_bayes_all,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = dist_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch distance", y = "Beta diversity")

#

p5 <- corals_taxonomic %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_past_isolation_sc),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    dist_sc = mean(dist_sc),
    diff_age_sc = mean(diff_age_sc),
    diff_past_sst = seq_range(diff_past_sst, n = 10)
  ) %>% 
  add_epred_draws(model_beta_corals_taxonomic_past_bayes_all,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_past_sst)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch SST", y = "Beta diversity")


library(cowplot)

p1 / p2 / p3 / p4 / p5





#2. Marine - CORALS BETA FUNCTIONAL ----


## Scaling variables - difference between marine regions ----
corals_functional <- model_corals_functional %>%
  mutate(
    beta_sor_adj = beta_sor - 0.000001,  
    diff_isolation_sc = scale(diff_isolation)[,1],
    diff_reef_area_sc = scale(diff_reef_area)[,1],
    diff_past_isolation_sc = scale(diff_past_isolation)[,1], #quaternary distance
    diff_past_reef_area_sc = scale(diff_past_reef_area)[,1],
    diff_age_sc = scale(diff_age)[,1],
    dist_sc = scale(Distance)[,1],
    diff_present_sst = scale(diff_sst)[,1],
    diff_past_sst = scale(diff_past_sst)[,1]
  )

corals_functional %>% glimpse

#------------------------- Models


# Bayes models
# Null

model_beta_corals_functional_past_null_gp2 <- brm(
  beta_sor_adj ~ 1 + 
    (1+dist_sc|Group1) + (1+dist_sc|Group2),
  data = corals_functional,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  control = list(adapt_delta = 0.99, max_treedepth = 14),
  save_pars = save_pars(all = T)
)

model_beta_corals_functional_past_null_gp2 <- add_criterion(model_beta_corals_functional_past_null_gp2, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_corals_functional_past_null_gp2,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_corals_functional_past_null_gp2.rds"
)

# Assess model validity using LOO
loo(model_beta_corals_functional_past_null_gp2)

# past
model_beta_corals_functional_past_bayes_all <- brm(
  beta_sor_adj ~ diff_past_isolation_sc + diff_past_reef_area_sc + diff_age_sc + dist_sc + diff_past_sst +
  (1+dist_sc|Group1) + (1+dist_sc|Group2),
  data = corals_functional,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_corals_functional_past_bayes_all <- add_criterion(model_beta_corals_functional_past_bayes_all, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_corals_functional_past_bayes_all,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_corals_functional_past_bayes_all.rds"
)

# Assess model validity using LOO
loo(model_beta_corals_functional_past_bayes_all)

#Look at model summary its on 'logit' scale, because of bernouli family link function
model_beta_corals_functional_past_bayes_all

#loo_compare(model_beta_corals_functional_past_bayes_all, model_beta_corals_functional_past_null_gp2)

## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_corals_functional_past_null_gp2 , ndraws = 50) 
pp_check(model_beta_corals_functional_past_bayes_all, type = "ecdf_overlay", ndraws = 50) 

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_corals_functional_past_bayes_all$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_corals_functional_past_bayes_all, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_corals_functional_past_bayes_all, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_corals_functional_past_bayes_all, type = 'rhat') 


#back transform model summary with point estimates - Now Odds ratios (3:1 odds of bleaching)
model_parameters(model_beta_corals_functional_past_bayes_all, exponentiate = T, ci = 0.95)


#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####

#Plotting isolation effects
p1 <- corals_functional %>% 
  modelr::data_grid(
    diff_past_isolation_sc = seq_range(diff_past_isolation_sc, n = 10),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_past_sst = mean(diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_functional_past_bayes_all,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_past_isolation_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in isolation", y = "Beta diversity")

#
p2 <- corals_functional %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_past_isolation_sc),
    diff_past_reef_area_sc = seq_range(diff_past_reef_area_sc, n = 10),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_past_sst = mean(diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_functional_past_bayes_all,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_past_reef_area_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in reef area", y = "Beta diversity")

# 
p3 <- corals_functional %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_past_isolation_sc),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    diff_age_sc = seq_range(diff_age_sc, n = 10),
    dist_sc = mean(dist_sc),
    diff_past_sst = mean(diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_functional_past_bayes_all,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_age_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in age", y = "Beta diversity")

# 
p4 <- corals_functional %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_past_isolation_sc),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    dist_sc = seq_range(dist_sc, n = 10),
    diff_age_sc = mean(diff_age_sc),
    diff_past_sst = mean(diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_functional_past_bayes_all,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = dist_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch distance", y = "Beta diversity")

#

p5 <- corals_functional %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_past_isolation_sc),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    dist_sc = mean(dist_sc),
    diff_age_sc = mean(diff_age_sc),
    diff_past_sst = seq_range(diff_past_sst, n = 10)
  ) %>% 
  add_epred_draws(model_beta_corals_functional_past_bayes_all,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_past_sst)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch SST", y = "Beta diversity")


library(cowplot)

p1 / p2 / p3 / p4 / p5



#3. Marine - FISH BETA TAXONOMIC ----


## Scaling variables - difference between marine regions ----
fish_taxonomic <- model_fish_taxonomic %>%
  mutate(
    beta_sor_adj = beta_sor - 0.000001,  
    diff_isolation_sc = scale(diff_isolation)[,1],
    diff_reef_area_sc = scale(diff_reef_area)[,1],
    diff_past_isolation_sc = scale(diff_past_isolation)[,1], #quaternary distance
    diff_past_reef_area_sc = scale(diff_past_reef_area)[,1],
    diff_age_sc = scale(diff_age)[,1],
    dist_sc = scale(Distance)[,1],
    diff_present_sst = scale(diff_sst)[,1],
    diff_past_sst = scale(diff_past_sst)[,1]
  )

fish_taxonomic %>% glimpse



#------------------------- Models


# Bayes models
# Null


# model_beta_fish_taxonomic_past_null_gp3 <- add_criterion(model_beta_fish_taxonomic_past_null_gp3, "loo", save_psis = TRUE, reloo = T)
# model_beta_fish_taxonomic_past_null_gp2 <- brm(
#   beta_sor_adj ~ 1 + (1|Group1) + (1|Group2) ,
#   data = fish_taxonomic,
#   family = beta_family(),
#   chains = 4,
#   cores = 4,
#   iter = 4000,
#   warmup = 3000,
#   control = list(adapt_delta = 0.99, max_treedepth = 14),
#   save_pars = save_pars(all = T)
# )

model_beta_fish_taxonomic_past_null_gp2 <- brm(
  beta_sor_adj ~ 1 + (1+dist_sc|Group1) + (1+dist_sc|Group2),
  data = fish_taxonomic,
  family = beta_family(),
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 3000,
  control = list(adapt_delta = 0.99, max_treedepth = 14),
  save_pars = save_pars(all = T)
)

model_beta_fish_taxonomic_past_null_gp2 <- add_criterion(model_beta_fish_taxonomic_past_null_gp2, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_fish_taxonomic_past_null_gp2,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_fish_taxonomic_past_null_gp2.rds"
)

# Assess model validity using LOO
loo(model_beta_fish_taxonomic_past_null_gp2)

#Look at model summary its on 'logit' scale, because of bernouli family link function
model_beta_fish_taxonomic_past_null_gp2


# past
model_beta_fish_taxonomic_past_bayes_all <- brm(
  beta_sor_adj ~ diff_past_isolation_sc + diff_past_reef_area_sc + diff_age_sc + dist_sc + diff_past_sst +
  (1+dist_sc|Group1) + (1+dist_sc|Group2),
  data = fish_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_fish_taxonomic_past_bayes_all <- add_criterion(model_beta_fish_taxonomic_past_bayes_all, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_fish_taxonomic_past_bayes_all,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_fish_taxonomic_past_bayes_all.rds"
)

model_beta_fish_taxonomic_past_bayes_all

# Assess model validity using LOO
loo(model_beta_fish_taxonomic_past_bayes_all)

#loo_compare(model_beta_fish_taxonomic_past_bayes_all, model_beta_fish_taxonomic_past_null_gp2)
## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_fish_taxonomic_past_bayes_all , ndraws = 50) 
pp_check(model_beta_fish_taxonomic_past_bayes_all, type = "ecdf_overlay", ndraws = 50) 

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_fish_taxonomic_past_bayes_all$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_fish_taxonomic_past_bayes_all, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_fish_taxonomic_past_bayes_all, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_fish_taxonomic_past_bayes_all, type = 'rhat') 


#back transform model summary with point estimates - Now Odds ratios (3:1 odds of bleaching)
model_parameters(model_beta_fish_taxonomic_past_bayes_all, exponentiate = T, ci = 0.95)


#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####

#Plotting isolation effects
p1 <- fish_taxonomic %>% 
  modelr::data_grid(
    diff_past_isolation_sc = seq_range(diff_past_isolation_sc, n = 10),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_past_sst = mean(diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_taxonomic_past_bayes_all,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_past_isolation_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in isolation", y = "Beta diversity")

#
p2 <- fish_taxonomic %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_past_isolation_sc),
    diff_past_reef_area_sc = seq_range(diff_past_reef_area_sc, n = 10),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_past_sst = mean(diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_taxonomic_past_bayes_all,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_past_reef_area_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in reef area", y = "Beta diversity")

# 
p3 <- fish_taxonomic %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_past_isolation_sc),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    diff_age_sc = seq_range(diff_age_sc, n = 10),
    dist_sc = mean(dist_sc),
    diff_past_sst = mean(diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_taxonomic_past_bayes_all,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_age_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in age", y = "Beta diversity")

# 
p4 <- fish_taxonomic %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_past_isolation_sc),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    dist_sc = seq_range(dist_sc, n = 10),
    diff_age_sc = mean(diff_age_sc),
    diff_past_sst = mean(diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_taxonomic_past_bayes_all,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = dist_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch distance", y = "Beta diversity")

#

p5 <- fish_taxonomic %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_past_isolation_sc),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    dist_sc = mean(dist_sc),
    diff_age_sc = mean(diff_age_sc),
    diff_past_sst = seq_range(diff_past_sst, n = 10)
  ) %>% 
  add_epred_draws(model_beta_fish_taxonomic_past_bayes_all,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_past_sst)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch SST", y = "Beta diversity")


library(cowplot)

p1 / p2 / p3 / p4 / p5





#4. Marine - FISH BETA FUNCTIONAL ----


## Scaling variables - difference between marine regions ----
fish_functional <- model_fish_functional %>%
  mutate(
    beta_sor_adj = beta_sor - 0.000001,  
    diff_isolation_sc = scale(diff_isolation)[,1],
    diff_reef_area_sc = scale(diff_reef_area)[,1],
    diff_past_isolation_sc = scale(diff_past_isolation)[,1], #quaternary distance
    diff_past_reef_area_sc = scale(diff_past_reef_area)[,1],
    diff_age_sc = scale(diff_age)[,1],
    dist_sc = scale(Distance)[,1],
    diff_present_sst = scale(diff_sst)[,1],
    diff_past_sst = scale(diff_past_sst)[,1]
  ) 

fish_functional %>% glimpse

#------------------------- Models


# Bayes models
# Null
model_beta_fish_functional_past_null_gp2 <- brm(
  beta_sor_adj ~ 1 + (1+dist_sc|Group1) + (1+dist_sc|Group2) ,
  data = fish_functional,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  control = list(adapt_delta = 0.99, max_treedepth = 14),
  save_pars = save_pars(all = T)
)

model_beta_fish_functional_past_null_gp2 <- add_criterion(model_beta_fish_functional_past_null_gp2, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_fish_functional_past_null_gp2,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_fish_functional_past_null_gp2.rds"
)


# Assess model validity using LOO
loo(model_beta_fish_functional_past_null_gp2)



# past
model_beta_fish_functional_past_bayes_all <- brm(
  beta_sor_adj ~  1 + diff_past_isolation_sc + diff_past_reef_area_sc + diff_age_sc + dist_sc + diff_past_sst + 
  (1+dist_sc|Group1) + (1+dist_sc|Group2),
  data = fish_functional,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_fish_functional_past_bayes_all <- add_criterion(model_beta_fish_functional_past_bayes_all, "loo", save_psis = TRUE, reloo = T, moment_match = TRUE)

saveRDS(
  model_beta_fish_functional_past_bayes_all,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_fish_functional_past_bayes_all.rds"
)

# Assess model validity using LOO
loo(model_beta_fish_functional_past_bayes_all)

#loo_compare(model_beta_fish_functional_past_bayes_all, model_beta_fish_functional_past_null_gp2)

## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_fish_functional_past_bayes_all , ndraws = 50) 
pp_check(model_beta_fish_functional_past_bayes_all, type = "ecdf_overlay", ndraws = 50) 
pp_check(model_beta_fish_functional_past_bayes_all , ndraws = 50, type = "dens_overlay_grouped", group = "Group1")
pp_check(model_beta_fish_functional_past_bayes_all , ndraws = 50, type = "scatter_avg_grouped", group = "Group1")
pp_check(model_beta_fish_functional_past_bayes_all , ndraws = 50, type = "error_scatter_avg_grouped", group = "Group1")
pp_check(model_beta_fish_functional_past_bayes_all, type='pit_ecdf_grouped', group = "Group1", prob = 0.95, plot_diff = F)
pp_check(model_beta_fish_functional_past_bayes_all, type='pit_ecdf_grouped', group = "Group2", prob = 0.95, plot_diff = F)

pp_check(model_beta_fish_functional_past_bayes_all , ndraws = 50, type = "scatter_avg_grouped", group = "Group2")

pairs(model_beta_fish_functional_past_bayes_all, variable = variables(model_beta_fish_functional_past_bayes_all)[1:5])
fish_functional %>% select(Group1, Group2) %>% unique %>% view

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_fish_functional_past_bayes_all$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_fish_functional_past_bayes_all, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_fish_functional_past_bayes_all, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_fish_functional_past_bayes_all, type = 'rhat') 


#back transform model summary with point estimates - Now Odds ratios (3:1 odds of bleaching)
model_parameters(model_beta_fish_functional_past_bayes_all, exponentiate = T, ci = 0.95)



#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####

#Plotting isolation effects
p1 <- fish_functional %>% 
  modelr::data_grid(
    diff_past_isolation_sc = seq_range(diff_past_isolation_sc, n = 10),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_past_sst = mean(diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_functional_past_bayes_all,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_past_isolation_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in isolation", y = "Beta diversity")

#
p2 <- fish_functional %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_past_isolation_sc),
    diff_past_reef_area_sc = seq_range(diff_past_reef_area_sc, n = 10),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_past_sst = mean(diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_functional_past_bayes_all,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_past_reef_area_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in reef area", y = "Beta diversity")

# 
p3 <- fish_functional %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_past_isolation_sc),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    diff_age_sc = seq_range(diff_age_sc, n = 10),
    dist_sc = mean(dist_sc),
    diff_past_sst = mean(diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_functional_past_bayes_all,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_age_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in age", y = "Beta diversity")

# 
p4 <- fish_functional %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_past_isolation_sc),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    dist_sc = seq_range(dist_sc, n = 10),
    diff_age_sc = mean(diff_age_sc),
    diff_past_sst = mean(diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_functional_past_bayes_all,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = dist_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch distance", y = "Beta diversity")

#
p5 <- fish_functional %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_past_isolation_sc),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    dist_sc = mean(dist_sc),
    diff_age_sc = mean(diff_age_sc),
    diff_past_sst = seq_range(diff_past_sst, n = 10)
  ) %>% 
  add_epred_draws(model_beta_fish_functional_past_bayes_all,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_past_sst)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch SST", y = "Beta diversity")


library(cowplot)

p1 / p2 / p3 / p4 / p5






#5. Terrestrial - PLANTS BETA TAXONOMIC ----


model_plants_taxonomic <- model_plants_taxonomic %>%
  filter(!Island1 %in% "Farquhar_Atoll",
         !Island2 %in% "Farquhar_Atoll")


# Scaling variables - difference between archipelagos islands
plants_taxonomic <- model_plants_taxonomic %>%
  mutate(
    beta_sor_adj = beta_sor - 0.000001,
    beta_sor_adj = ifelse(beta_sor_adj <= 0, 1e-6, beta_sor_adj),
    d_past_sc = scale(d_past.y), #emerse distance past
    diff_emerse_past_sc = scale(diff_emerse_past), #emerse island area past
    diff_emerse_pres_sc = scale(diff_emerse_pres), #emerse island area past
    diff_island_age_sc = scale(diff_island_age),
    distance_km_sc = scale(Distance_km), #emerse distance past
    diff_present_air = scale(diif_present_air),
    diff_past_air = scale(diff_past_air)
  )

# scaled_diff_emerse_pres <- plants_taxonomic$diff_emerse_pres - mean(plants_taxonomic$diff_emerse_pres)/sd(plants_taxonomic$diff_emerse_pres)
# 
# diff_emerse_pres <- scaled_diff_emerse_pres * sd(plants_taxonomic$diff_emerse_pres) + mean(plants_taxonomic$diff_emerse_pres)
# 
# diff_emerse_pres_sc <-  scale(plants_taxonomic$diff_island_age_sc)
# diff_emerse_pres_sc %>% glimpse


#------------------------- Models

# Null

model_beta_plants_taxonomic_past_null<- brm(
  beta_sor_adj ~ 1 +
    (1 | Archipelago),
  data = plants_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 5000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.98, max_treedepth = 14),
)


model_beta_plants_taxonomic_past_null<- add_criterion(model_beta_plants_taxonomic_past_null, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_plants_taxonomic_past_null,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_plants_taxonomic_past_null.rds"
)

# Assess model validity using LOO
loo(model_beta_plants_taxonomic_past_null)



# past

model_beta_plants_taxonomic_past <- brm(
  beta_sor_adj ~  diff_emerse_past_sc + diff_island_age_sc + d_past_sc + diff_past_air + 
    (1 | Archipelago),
  data = plants_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_plants_taxonomic_past <- add_criterion(model_beta_plants_taxonomic_past, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_plants_taxonomic_past,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_plants_taxonomic_past.rds"
)

loo(model_beta_plants_taxonomic_past)

#loo_compare(model_beta_plants_taxonomic_past_null, model_beta_plants_taxonomic_past)

#Look at model summary its on 'logit' scale, because of Beta family link function
model_beta_plants_taxonomic_past

## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_plants_taxonomic_past, ndraws = 50) 
pp_check(model_beta_plants_taxonomic_past, type = "ecdf_overlay", ndraws = 50) 

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_plants_taxonomic_past$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_plants_taxonomic_past, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_plants_taxonomic_past, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_plants_taxonomic_past, type = 'rhat') 


#back transform model summary with point estimates (Ratios)
model_parameters(model_beta_plants_taxonomic_past, exponentiate = T, ci = 0.95)


#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####


#Plotting isolation effects
p1 <- plants_taxonomic %>% 
  modelr::data_grid(d_past_sc = seq_range(d_past_sc, n = 100) %>% as.numeric,
                    diff_emerse_past_sc = median(diff_emerse_past_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_past_air = median(diff_past_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_plants_taxonomic_past, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = d_past_sc)) +
  # geom_point(data = corals_taxonomic$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Distance between islands", y = "Beta diversity")

p2 <- plants_taxonomic %>% 
  modelr::data_grid(d_past_sc = median(d_past_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_past_sc = seq_range(diff_emerse_past_sc, n = 10),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_past_air = median(diff_past_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_plants_taxonomic_past, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_emerse_past_sc)) +
  # geom_point(data = corals_taxonomic$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Emersed Island Area", y = "Beta diversity")

p3 <- plants_taxonomic %>% 
  modelr::data_grid(d_past_sc = median(d_past_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_past_sc = median(diff_emerse_past_sc, na.rm = T),
                    diff_island_age_sc = seq_range(diff_island_age_sc, n = 10),
                    diff_past_air = median(diff_past_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_plants_taxonomic_past, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_island_age_sc)) +
  # geom_point(data = corals_taxonomic$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Island Age", y = "Beta diversity")


p4 <- plants_taxonomic %>% 
  modelr::data_grid(d_past_sc = median(d_past_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_past_sc = median(diff_emerse_past_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_past_air = seq_range(diff_past_air, n = 10)) %>% 
  add_epred_draws(model_beta_plants_taxonomic_past, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_past_air)) +
  # geom_point(data = corals_taxonomic$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Air Temperature", y = "Beta diversity")

library(cowplot)

p1 / p2 / p3 / p4





#6. Terrestrial - PLANTS BETA FUNCTIONAL ----


model_plants_functional <- model_plants_functional %>%
  filter(!Island1 %in% "Farquhar_Atoll",
         !Island2 %in% "Farquhar_Atoll")



# Scaling variables - difference between archipelagos islands
plants_functional <- model_plants_functional %>%
  mutate(
    beta_sor_adj = beta_sor - 0.000001, 
    beta_sor_adj = ifelse(beta_sor_adj <= 0, 1e-6, beta_sor_adj),
    d_past_sc = scale(d_past.y), #emerse distance past
    diff_emerse_past_sc = scale(diff_emerse_past), #emerse island area past
    diff_emerse_pres_sc = scale(diff_emerse_pres), #emerse island area past
    diff_island_age_sc = scale(diff_island_age),
    distance_km_sc = scale(Distance_km), #emerse distance past
    diff_present_air = scale(diif_present_air),
    diff_past_air = scale(diff_past_air)
  )


#------------------------- Models

# Null

model_beta_plants_functional_past_null<- brm(
  beta_sor_adj ~ 1 +
    (1 | Archipelago),
  data = plants_functional,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 5000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.98, max_treedepth = 14),
)


model_beta_plants_functional_past_null<- add_criterion(model_beta_plants_functional_past_null, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_plants_functional_past_null,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_plants_functional_past_null.rds"
)

# Assess model validity using LOO
loo(model_beta_plants_functional_past_null)



# past

model_beta_plants_functional_past <- brm(
  beta_sor_adj ~  diff_emerse_past_sc + diff_island_age_sc + d_past_sc + diff_past_air +
    (1 | Archipelago),
  data = plants_functional,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_plants_functional_past <- add_criterion(model_beta_plants_functional_past, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_plants_functional_past,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_plants_functional_past.rds"
)

loo(model_beta_plants_functional_past)

#loo_compare(model_beta_plants_functional_past_null, model_beta_plants_functional_past)

#Look at model summary its on 'logit' scale, because of Beta family link function
model_beta_plants_functional_past

## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_plants_functional_past, ndraws = 50) 
pp_check(model_beta_plants_functional_past, type = "ecdf_overlay", ndraws = 50) 

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_plants_functional_past$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_plants_functional_past, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_plants_functional_past, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_plants_functional_past, type = 'rhat') 


#back transform model summary with point estimates (Ratios)
model_parameters(model_beta_plants_functional_past, exponentiate = T, ci = 0.95)


#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####


#Plotting isolation effects
p1 <- plants_functional %>% 
  modelr::data_grid(d_past_sc = seq_range(d_past_sc, n = 100) %>% as.numeric,
                    diff_emerse_past_sc = median(diff_emerse_past_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_past_air = median(diff_past_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_plants_functional_past, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = d_past_sc)) +
  # geom_point(data = corals_functional$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Distance between islands", y = "Beta diversity")

p2 <- plants_functional %>% 
  modelr::data_grid(d_past_sc = median(d_past_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_past_sc = seq_range(diff_emerse_past_sc, n = 10),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_past_air = median(diff_past_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_plants_functional_past, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_emerse_past_sc)) +
  # geom_point(data = corals_functional$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Emersed Island Area", y = "Beta diversity")

p3 <- plants_functional %>% 
  modelr::data_grid(d_past_sc = median(d_past_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_past_sc = median(diff_emerse_past_sc, na.rm = T),
                    diff_island_age_sc = seq_range(diff_island_age_sc, n = 10),
                    diff_past_air = median(diff_past_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_plants_functional_past, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_island_age_sc)) +
  # geom_point(data = corals_functional$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Island Age", y = "Beta diversity")


p4 <- plants_functional %>% 
  modelr::data_grid(d_past_sc = median(d_past_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_past_sc = median(diff_emerse_past_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_past_air = seq_range(diff_past_air, n = 10)) %>% 
  add_epred_draws(model_beta_plants_functional_past, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_past_air)) +
  # geom_point(data = corals_taxonomic$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Air Temperature", y = "Beta diversity")

library(cowplot)

p1 / p2 / p3 / p4


#7. Terrestrial - BIRDS BETA TAXONOMIC ----


model_birds_taxonomic <- model_birds_taxonomic %>%
  filter(!Island1 %in% "Farquhar_Atoll",
         !Island2 %in% "Farquhar_Atoll")



# Scaling variables - difference between archipelagos islands
birds_taxonomic <- model_birds_taxonomic %>%
  mutate(
    beta_sor_adj = beta_sor - 0.000001,
    beta_sor_adj = ifelse(beta_sor_adj <= 0, 1e-6, beta_sor_adj),
    d_past_sc = scale(d_past.y), #emerse distance past
    diff_emerse_past_sc = scale(diff_emerse_past), #emerse island area past
    diff_emerse_pres_sc = scale(diff_emerse_pres), #emerse island area past
    diff_island_age_sc = scale(diff_island_age),
    distance_km_sc = scale(Distance_km), #emerse distance past
    diff_present_air = scale(diif_present_air),
    diff_past_air = scale(diff_past_air)
  )




#------------------------- Models

# Null

model_beta_birds_taxonomic_past_null <- brm(
  beta_sor_adj ~ 1 +
    (1 | Archipelago),
  data = birds_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 5000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.98, max_treedepth = 14),
)


model_beta_birds_taxonomic_past_null <- add_criterion(model_beta_birds_taxonomic_past_null, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_birds_taxonomic_past_null,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_birds_taxonomic_past_null.rds"
)

# Assess model validity using LOO
loo(model_beta_birds_taxonomic_past_null)


# past

model_beta_birds_taxonomic_past <- brm(
  beta_sor_adj ~  diff_emerse_past_sc + diff_island_age_sc + d_past_sc + diff_past_air +
    (1 | Archipelago),
  data = birds_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_birds_taxonomic_past <- add_criterion(model_beta_birds_taxonomic_past, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_birds_taxonomic_past,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_birds_taxonomic_past.rds"
)


loo(model_beta_birds_taxonomic_past)

#loo_compare(model_beta_birds_taxonomic_past_null, model_beta_birds_taxonomic_past)

#Look at model summary its on 'logit' scale, because of Beta family link function
model_beta_birds_taxonomic_past

## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_birds_taxonomic_past, ndraws = 50) 

pp_check(model_beta_birds_taxonomic_past, type = "ecdf_overlay", ndraws = 50) 

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_birds_taxonomic_past$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_birds_taxonomic_past, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_birds_taxonomic_past, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_birds_taxonomic_past, type = 'rhat') 


#back transform model summary with point estimates (Ratios)
model_parameters(model_beta_birds_taxonomic_past, exponentiate = T, ci = 0.95)


#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####


#Plotting isolation effects
p1 <- birds_taxonomic %>% 
  modelr::data_grid(d_past_sc = seq_range(d_past_sc, n = 100) %>% as.numeric,
                    diff_emerse_past_sc = median(diff_emerse_past_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_past_air = median(diff_past_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_birds_taxonomic_past, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = d_past_sc)) +
  # geom_point(data = corals_taxonomic$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Distance between islands", y = "Beta diversity")

p2 <- birds_taxonomic %>% 
  modelr::data_grid(d_past_sc = median(d_past_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_past_sc = seq_range(diff_emerse_past_sc, n = 10),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_past_air = median(diff_past_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_birds_taxonomic_past, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_emerse_past_sc)) +
  # geom_point(data = corals_taxonomic$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Emersed Island Area", y = "Beta diversity")

p3 <- birds_taxonomic %>% 
  modelr::data_grid(d_past_sc = median(d_past_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_past_sc = median(diff_emerse_past_sc, na.rm = T),
                    diff_island_age_sc = seq_range(diff_island_age_sc, n = 10),
                    diff_past_air = median(diff_past_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_birds_taxonomic_past, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_island_age_sc)) +
  # geom_point(data = corals_taxonomic$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Island Age", y = "Beta diversity")

p4 <- birds_taxonomic %>% 
  modelr::data_grid(d_past_sc = median(d_past_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_past_sc = median(diff_emerse_past_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_past_air = seq_range(diff_past_air, n = 10)) %>% 
  add_epred_draws(model_beta_birds_taxonomic_past, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_past_air)) +
  # geom_point(data = corals_taxonomic$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Air Temperature", y = "Beta diversity")

library(cowplot)

p1 / p2 / p3 / p4





#8. Terrestrial - BIRDS BETA FUNCTIONAL ----


model_birds_functional <- model_birds_functional %>%
  filter(!Island1 %in% "Farquhar_Atoll",
         !Island2 %in% "Farquhar_Atoll")



# Scaling variables - difference between archipelagos islands
birds_functional <- model_birds_functional %>%
  mutate(
    beta_sor_adj = beta_sor - 0.000001,
    beta_sor_adj = ifelse(beta_sor_adj <= 0, 1e-6, beta_sor_adj),
    d_past_sc = scale(d_past.y), #emerse distance past
    diff_emerse_past_sc = scale(diff_emerse_past), #emerse island area past
    diff_emerse_pres_sc = scale(diff_emerse_pres), #emerse island area past
    diff_island_age_sc = scale(diff_island_age),
    distance_km_sc = scale(Distance_km), #emerse distance past
    diff_present_air = scale(diif_present_air),
    diff_past_air = scale(diff_past_air)
  )




#------------------------- Models

# Null

model_beta_birds_functional_past_null <- brm(
  beta_sor_adj ~ 1 +
    (1 | Archipelago),
  data = birds_functional,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 5000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.98, max_treedepth = 14),
)


model_beta_birds_functional_past_null <- add_criterion(model_beta_birds_functional_past_null, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_birds_functional_past_null,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_birds_functional_past_null.rds"
)

# Assess model validity using LOO
loo(model_beta_birds_functional_past_null)


# past

model_beta_birds_functional_past <- brm(
  beta_sor_adj ~  diff_emerse_past_sc + diff_island_age_sc + d_past_sc + diff_past_air +
    (1 | Archipelago),
  data = birds_functional,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_birds_functional_past <- add_criterion(model_beta_birds_functional_past, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_birds_functional_past,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_birds_functional_past.rds"
)

loo(model_beta_birds_functional_past)

#loo_compare(model_beta_birds_functional_past_null, model_beta_birds_functional_past)

#Look at model summary its on 'logit' scale, because of Beta family link function
model_beta_birds_functional_past

## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_birds_functional_past, ndraws = 50) 
pp_check(model_beta_birds_functional_past, type = "ecdf_overlay", ndraws = 50) 

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_birds_functional_past$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_birds_functional_past, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_birds_functional_past, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_birds_functional_past, type = 'rhat') 


#back transform model summary with point estimates (Ratios)
model_parameters(model_beta_birds_functional_past, exponentiate = T, ci = 0.95)


#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####


#Plotting isolation effects
p1 <- birds_functional %>% 
  modelr::data_grid(d_past_sc = seq_range(d_past_sc, n = 100) %>% as.numeric,
                    diff_emerse_past_sc = median(diff_emerse_past_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_past_air = median(diff_past_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_birds_functional_past, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = d_past_sc)) +
  # geom_point(data = corals_functional$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Distance between islands", y = "Beta diversity")

p2 <- birds_functional %>% 
  modelr::data_grid(d_past_sc = median(d_past_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_past_sc = seq_range(diff_emerse_past_sc, n = 10),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_past_air = median(diff_past_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_birds_functional_past, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_emerse_past_sc)) +
  # geom_point(data = corals_functional$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Emersed Island Area", y = "Beta diversity")

p3 <- birds_functional %>% 
  modelr::data_grid(d_past_sc = median(d_past_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_past_sc = median(diff_emerse_past_sc, na.rm = T),
                    diff_island_age_sc = seq_range(diff_island_age_sc, n = 10),
                    diff_past_air = median(diff_past_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_birds_functional_past, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_island_age_sc)) +
  # geom_point(data = corals_functional$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Island Age", y = "Beta diversity")

p4 <- birds_functional %>% 
  modelr::data_grid(d_past_sc = median(d_past_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_past_sc = median(diff_emerse_past_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_past_air = seq_range(diff_past_air, n = 10)) %>% 
  add_epred_draws(model_beta_birds_functional_past, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_past_air)) +
  # geom_point(data = corals_taxonomic$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Air Temperature", y = "Beta diversity")

library(cowplot)

p1 / p2 / p3 / p4







######################################################################################################################################################################


# NESTEDNESS ----

#-------------------------------------------------------------------------------


#1. Marine - CORALS BETA TAXONOMIC ----


## Scaling variables - difference between marine regions ----
corals_taxonomic <- model_corals_taxonomic %>%
  mutate(
    beta_sne_adj = beta_sne - 0.000001, 
    beta_sne_adj = ifelse(beta_sne_adj <= 0, 1e-6, beta_sne_adj),
    diff_isolation_sc = scale(diff_isolation)[,1],
    diff_reef_area_sc = scale(diff_reef_area)[,1],
    diff_past_isolation_sc = scale(diff_past_isolation)[,1], #quaternary distance
    diff_past_reef_area_sc = scale(diff_past_reef_area)[,1],
    diff_age_sc = scale(diff_age)[,1],
    dist_sc = scale(Distance)[,1],
    diff_present_sst = scale(diff_sst)[,1],
    diff_past_sst = scale(diff_past_sst)[,1]
  )

corals_taxonomic %>% glimpse

#------------------------- Models

# Present


# Bayes models
# Null

# model_beta_corals_taxonomic_present_null_nest_nestgp3 <- brm(
#   beta_sne_adj ~ 1 + (1|Group1) + (1|Group2) + (1|key),
#   data = corals_taxonomic,
#   family = beta_family(),
#   chains = 4,
#   cores = 4,
#   iter = 4000,
#   warmup = 3000,
#   control = list(adapt_delta = 0.99, max_treedepth = 14),
#   save_pars = save_pars(all = T)
# )

# model_beta_corals_taxonomic_present_null_nest_nest_gp3 <- add_criterion(model_beta_corals_taxonomic_present_null_nest_nest_gp3, "loo", save_psis = TRUE, reloo = T)


# model_beta_corals_taxonomic_present_null <- brm(
#   beta_sne_adj ~ 1,
#   data = corals_taxonomic,
#   family = beta_family(),
#   chains = 4,
#   cores = 4,
#   iter = 4000,
#   warmup = 3000,
#   control = list(adapt_delta = 0.99, max_treedepth = 14),
#   save_pars = save_pars(all = T)
# )

# model_beta_corals_taxonomic_present_null <- add_criterion(model_beta_corals_taxonomic_present_null, "loo", save_psis = TRUE, reloo = T)
# 
# loo_compare(model_beta_corals_taxonomic_present_null, model_beta_corals_taxonomic_present_null_nest_nest_gp2, model_beta_corals_taxonomic_present_null_nest_nest_gp1, model_beta_corals_taxonomic_present_null_nest_nest_gp3)

model_beta_corals_taxonomic_present_null_gp2 <- brm(
  beta_sne_adj ~ 1 + 
    (1+dist_sc|Group1) + (1+dist_sc|Group2),
  data = corals_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  control = list(adapt_delta = 0.99, max_treedepth = 14),
  save_pars = save_pars(all = T)
)

model_beta_corals_taxonomic_present_null_gp2 <- add_criterion(model_beta_corals_taxonomic_present_null_gp2, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_corals_taxonomic_present_null_gp2,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_corals_taxonomic_present_null_gp2.rds"
)


# Assess model validity using LOO
loo(model_beta_corals_taxonomic_present_null_gp2)

# Present
model_beta_corals_taxonomic_present_bayes_nest <- brm(
  beta_sne_adj ~ diff_isolation_sc + diff_reef_area_sc + diff_age_sc + dist_sc + diff_present_sst +
    (1+dist_sc|Group1) + (1+dist_sc|Group2),
  data = corals_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_corals_taxonomic_present_bayes_nest <- add_criterion(model_beta_corals_taxonomic_present_bayes_nest, "loo", save_psis = TRUE, reloo = T)

model_beta_corals_taxonomic_present_bayes_nest

saveRDS(
  model_beta_corals_taxonomic_present_bayes_nest,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_corals_taxonomic_present_bayes_nest.rds"
)

# Assess model validity using LOO
loo(model_beta_corals_taxonomic_present_bayes_nest)

#Look at model summary its on 'logit' scale, because of bernouli family link function
#model_beta_corals_taxonomic_present_bayes_2

#loo_compare(model_beta_corals_taxonomic_present_bayes_nest, model_beta_corals_taxonomic_present_null_gp2)
## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_corals_taxonomic_present_bayes_nest , ndraws = 50) 
pp_check(model_beta_corals_taxonomic_present_bayes_nest, type = "ecdf_overlay", ndraws = 50) 

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_corals_taxonomic_present_bayes_nest$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_corals_taxonomic_present_bayes_nest, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_corals_taxonomic_present_bayes_nest, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_corals_taxonomic_present_bayes_nest, type = 'rhat') 


#back transform model summary with point estimates - Now Odds ratios (3:1 odds of bleaching)
model_parameters(model_beta_corals_taxonomic_present_bayes_nest, exponentiate = T, ci = 0.95)



#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####

#Plotting isolation effects
# 
p1 <- corals_taxonomic %>% 
  modelr::data_grid(
    diff_isolation_sc = seq_range(diff_isolation_sc, n = 10),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_present_sst = mean(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_taxonomic_present_bayes_nest,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_isolation_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in isolation", y = "Beta diversity")

#
p2 <- corals_taxonomic %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = seq_range(diff_reef_area_sc, n = 10),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_present_sst = mean(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_taxonomic_present_bayes_nest,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_reef_area_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in reef area", y = "Beta diversity")

# 
p3 <- corals_taxonomic %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    diff_age_sc = seq_range(diff_age_sc, n = 10),
    dist_sc = mean(dist_sc),
    diff_present_sst = mean(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_taxonomic_present_bayes_nest,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_age_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in age", y = "Beta diversity")

# 
p4 <- corals_taxonomic %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    dist_sc = seq_range(dist_sc, n = 10),
    diff_age_sc = mean(diff_age_sc),
    diff_present_sst = mean(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_taxonomic_present_bayes_nest,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = dist_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch distance", y = "Beta diversity")

#
p5 <- corals_taxonomic %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    dist_sc = mean(dist_sc),
    diff_age_sc = mean(diff_age_sc),
    diff_present_sst = seq_range(diff_present_sst, n = 10)
  ) %>% 
  add_epred_draws(model_beta_corals_taxonomic_present_bayes_nest,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_present_sst)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch SST", y = "Beta diversity")


library(cowplot)

p1 / p2 / p3 / p4 / p5






#2. Marine - CORALS BETA FUNCTIONAL ----


## Scaling variables - difference between marine regions ----
corals_functional <- model_corals_functional %>%
  mutate(
    beta_sne_adj = beta_sne - 0.000001, 
    beta_sne_adj = ifelse(beta_sne_adj <= 0, 1e-6, beta_sne_adj),
    diff_isolation_sc = scale(diff_isolation)[,1],
    diff_reef_area_sc = scale(diff_reef_area)[,1],
    diff_past_isolation_sc = scale(diff_past_isolation)[,1], #quaternary distance
    diff_past_reef_area_sc = scale(diff_past_reef_area)[,1],
    diff_age_sc = scale(diff_age)[,1],
    dist_sc = scale(Distance)[,1],
    diff_present_sst = scale(diff_sst)[,1],
    diff_past_sst = scale(diff_sst)[,1]
  )

corals_functional %>% glimpse

#------------------------- Models


# Bayes models
# Null

model_beta_corals_functional_present_null_nest_gp2 <- brm(
  beta_sne_adj ~ 1 + 
    (1+dist_sc|Group1) + (1+dist_sc|Group2),
  data = corals_functional,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  control = list(adapt_delta = 0.99, max_treedepth = 14),
  save_pars = save_pars(all = T)
)

model_beta_corals_functional_present_null_nest_gp2 <- add_criterion(model_beta_corals_functional_present_null_nest_gp2, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_corals_functional_present_null_nest_gp2,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_corals_functional_present_null_nest_gp2.rds"
)

# Assess model validity using LOO
loo(model_beta_corals_functional_present_null_nest_gp2)

# Present
model_beta_corals_functional_present_bayes_nest <- brm(
  beta_sne_adj ~ diff_isolation_sc + diff_reef_area_sc + diff_age_sc + dist_sc + diff_present_sst + 
  (1+dist_sc|Group1) + (1+dist_sc|Group2),
  data = corals_functional,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_corals_functional_present_bayes_nest <- add_criterion(model_beta_corals_functional_present_bayes_nest, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_corals_functional_present_bayes_nest,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_corals_functional_present_bayes_nest.rds"
)

# Assess model validity using LOO
loo(model_beta_corals_functional_present_bayes_nest)

#Look at model summary its on 'logit' scale, because of bernouli family link function
model_beta_corals_functional_present_bayes_nest

#loo_compare(model_beta_corals_functional_present_bayes_nest, model_beta_corals_functional_present_null_nest_gp2)

## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_corals_functional_present_null_nest_gp2 , ndraws = 50) 
pp_check(model_beta_corals_functional_present_bayes_nest, type = "ecdf_overlay", ndraws = 50) 

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_corals_functional_present_bayes_nest$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_corals_functional_present_bayes_nest, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_corals_functional_present_bayes_nest, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_corals_functional_present_bayes_nest, type = 'rhat') 


#back transform model summary with point estimates - Now Odds ratios (3:1 odds of bleaching)
model_parameters(model_beta_corals_functional_present_bayes_nest, exponentiate = T, ci = 0.95)


#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####

#Plotting isolation effects
p1 <- corals_functional %>% 
  modelr::data_grid(
    diff_isolation_sc = seq_range(diff_isolation_sc, n = 10),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_present_sst = mean(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_functional_present_bayes_nest,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_isolation_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in isolation", y = "Beta diversity")

#
p2 <- corals_functional %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = seq_range(diff_reef_area_sc, n = 10),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_present_sst = mean(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_functional_present_bayes_nest,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_reef_area_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in reef area", y = "Beta diversity")

# 
p3 <- corals_functional %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    diff_age_sc = seq_range(diff_age_sc, n = 10),
    dist_sc = mean(dist_sc),
    diff_present_sst = mean(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_functional_present_bayes_nest,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_age_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in age", y = "Beta diversity")

# 
p4 <- corals_taxonomic %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    dist_sc = seq_range(dist_sc, n = 10),
    diff_age_sc = mean(diff_age_sc),
    diff_present_sst = mean(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_functional_present_bayes_nest,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = dist_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch distance", y = "Beta diversity")

#
p5 <- corals_taxonomic %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    dist_sc = mean(dist_sc),
    diff_age_sc = mean(diff_age_sc),
    diff_present_sst = seq_range(diff_present_sst, n = 10)
  ) %>% 
  add_epred_draws(model_beta_corals_functional_present_bayes_nest,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_present_sst)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch SST", y = "Beta diversity")


library(cowplot)

p1 / p2 / p3 / p4 / p5




#3. Marine - FISH BETA TAXONOMIC ----


## Scaling variables - difference between marine regions ----
fish_taxonomic <- model_fish_taxonomic %>%
  mutate(
    beta_sne_adj = beta_sne - 0.000001,
    beta_sne_adj = ifelse(beta_sne_adj <= 0, 1e-6, beta_sne_adj),
    diff_isolation_sc = scale(diff_isolation)[,1],
    diff_reef_area_sc = scale(diff_reef_area)[,1],
    diff_past_isolation_sc = scale(diff_past_isolation)[,1], #quaternary distance
    diff_past_reef_area_sc = scale(diff_past_reef_area)[,1],
    diff_age_sc = scale(diff_age)[,1],
    dist_sc = scale(Distance)[,1],
    diff_present_sst = scale(diff_sst)[,1],
    diff_past_sst = scale(diff_past_sst)[,1]
  )

fish_taxonomic %>% glimpse



#------------------------- Models


# Bayes models
# Null


# model_beta_fish_taxonomic_present_null_nest_gp3 <- add_criterion(model_beta_fish_taxonomic_present_null_nest_gp3, "loo", save_psis = TRUE, reloo = T)
# model_beta_fish_taxonomic_present_null_nest_gp2 <- brm(
#   beta_sne_adj ~ 1 + (1|Group1) + (1|Group2) ,
#   data = fish_taxonomic,
#   family = beta_family(),
#   chains = 4,
#   cores = 4,
#   iter = 4000,
#   warmup = 3000,
#   control = list(adapt_delta = 0.99, max_treedepth = 14),
#   save_pars = save_pars(all = T)
# )

model_beta_fish_taxonomic_present_null_nest_gp2 <- brm(
  beta_sne_adj ~ 1 + (1+dist_sc|Group1) + (1+dist_sc|Group2),
  data = fish_taxonomic,
  family = beta_family(),
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 3000,
  control = list(adapt_delta = 0.99, max_treedepth = 14),
  save_pars = save_pars(all = T)
)

model_beta_fish_taxonomic_present_null_nest_gp2 <- add_criterion(model_beta_fish_taxonomic_present_null_nest_gp2, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_fish_taxonomic_present_null_nest_gp2,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_fish_taxonomic_present_null_nest_gp2.rds"
)

# Assess model validity using LOO
loo(model_beta_fish_taxonomic_present_null_nest_gp2)

#Look at model summary its on 'logit' scale, because of bernouli family link function
model_beta_fish_taxonomic_present_null_nest_gp2


# Present
model_beta_fish_taxonomic_present_bayes_nest <- brm(
  beta_sne_adj ~ diff_isolation_sc + diff_reef_area_sc + diff_age_sc + dist_sc + diff_present_sst + 
  (1+dist_sc|Group1) + (1+dist_sc|Group2),
  data = fish_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_fish_taxonomic_present_bayes_nest <- add_criterion(model_beta_fish_taxonomic_present_bayes_nest, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_fish_taxonomic_present_bayes_nest,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_fish_taxonomic_present_bayes_nest.rds"
)

model_beta_fish_taxonomic_present_bayes_nest

# Assess model validity using LOO
loo(model_beta_fish_taxonomic_present_bayes_nest)

#loo_compare(model_beta_fish_taxonomic_present_bayes_nest, model_beta_fish_taxonomic_present_null_nest_gp2)
## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_fish_taxonomic_present_bayes_nest , ndraws = 50) 
pp_check(model_beta_fish_taxonomic_present_bayes_nest, type = "ecdf_overlay", ndraws = 50) 

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_fish_taxonomic_present_bayes_nest$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_fish_taxonomic_present_bayes_nest, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_fish_taxonomic_present_bayes_nest, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_fish_taxonomic_present_bayes_nest, type = 'rhat') 


#back transform model summary with point estimates - Now Odds ratios (3:1 odds of bleaching)
model_parameters(model_beta_fish_taxonomic_present_bayes_nest, exponentiate = T, ci = 0.95)


#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####

#Plotting isolation effects
p1 <- fish_taxonomic %>% 
  modelr::data_grid(
    diff_isolation_sc = seq_range(diff_isolation_sc, n = 10),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_present_sst = mean(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_taxonomic_present_bayes_nest,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_isolation_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in isolation", y = "Beta diversity")

#
p2 <- fish_taxonomic %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = seq_range(diff_reef_area_sc, n = 10),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_present_sst = mean(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_taxonomic_present_bayes_nest,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_reef_area_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in reef area", y = "Beta diversity")

# 
p3 <- fish_taxonomic %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    diff_age_sc = seq_range(diff_age_sc, n = 10),
    dist_sc = mean(dist_sc),
    diff_present_sst = mean(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_taxonomic_present_bayes_nest,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_age_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in age", y = "Beta diversity")

# 
p4 <- fish_taxonomic %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    dist_sc = seq_range(dist_sc, n = 10),
    diff_age_sc = mean(diff_age_sc),
    diff_present_sst = mean(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_taxonomic_present_bayes_nest,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = dist_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch distance", y = "Beta diversity")

#
p5 <- fish_taxonomic %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    dist_sc = mean(dist_sc),
    diff_age_sc = mean(diff_age_sc),
    diff_present_sst = seq_range(diff_present_sst, n = 10)
  ) %>% 
  add_epred_draws(model_beta_fish_taxonomic_present_bayes_nest,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_present_sst)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch SST", y = "Beta diversity")


library(cowplot)

p1 / p2 / p3 / p4 / p5




#4. Marine - FISH BETA FUNCTIONAL ----


## Scaling variables - difference between marine regions ----
fish_functional <- model_fish_functional %>%
  mutate(
    beta_sne_adj = beta_sne - 0.000001,  
    beta_sne_adj = ifelse(beta_sne_adj <= 0, 1e-6, beta_sne_adj),
    diff_isolation_sc = scale(diff_isolation)[,1],
    diff_reef_area_sc = scale(diff_reef_area)[,1],
    diff_past_isolation_sc = scale(diff_past_isolation)[,1], #quaternary distance
    diff_past_reef_area_sc = scale(diff_past_reef_area)[,1],
    diff_age_sc = scale(diff_age)[,1],
    dist_sc = scale(Distance)[,1],
    diff_present_sst = scale(diff_sst)[,1],
    diff_past_sst = scale(diff_past_sst)[,1]
  ) 

fish_functional %>% glimpse

#------------------------- Models


# Bayes models
# Null
model_beta_fish_functional_present_null_nest_gp2 <- brm(
  beta_sne_adj ~ 1 + (1+dist_sc|Group1) + (1+dist_sc|Group2) ,
  data = fish_functional,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  control = list(adapt_delta = 0.99, max_treedepth = 14),
  save_pars = save_pars(all = T)
)

model_beta_fish_functional_present_null_nest_gp2 <- add_criterion(model_beta_fish_functional_present_null_nest_gp2, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_fish_functional_present_null_nest_gp2,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_fish_functional_present_null_nest_gp2.rds"
)

# Assess model validity using LOO
loo(model_beta_fish_functional_present_null_nest_gp2)



# Present
model_beta_fish_functional_present_bayes_nest <- brm(
  beta_sne_adj ~  1 + diff_isolation_sc + diff_reef_area_sc + diff_age_sc + dist_sc + diff_present_sst + 
  (1+dist_sc|Group1) + (1+dist_sc|Group2),
  data = fish_functional,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)



# Add Leave One Out cross validation criterion
model_beta_fish_functional_present_bayes_nest <- add_criterion(model_beta_fish_functional_present_bayes_nest, "loo", save_psis = TRUE, reloo = T, moment_match = TRUE)

saveRDS(
  model_beta_fish_functional_present_bayes_nest,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_fish_functional_present_bayes_nest.rds"
)

# Assess model validity using LOO
loo(model_beta_fish_functional_present_bayes_nest)

#loo_compare(model_beta_fish_functional_present_bayes_nest, model_beta_fish_functional_present_null_nest_gp2)

## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_fish_functional_present_bayes_nest , ndraws = 50) 
pp_check(model_beta_fish_functional_present_bayes_nest, type = "ecdf_overlay", ndraws = 50) 
pp_check(model_beta_fish_functional_present_bayes_nest , ndraws = 50, type = "dens_overlay_grouped", group = "Group1")
# pp_check(model_beta_fish_functional_present_bayes_nest_g , ndraws = 50, type = "scatter_avg_grouped", group = "Group1")
# pp_check(model_beta_fish_functional_present_bayes_nest_g , ndraws = 50, type = "error_scatter_avg_grouped", group = "Group1")
# pp_check(model_beta_fish_functional_present_bayes_nest_g, type='pit_ecdf_grouped', group = "Group1", prob = 0.95, plot_diff = F)
# pp_check(model_beta_fish_functional_present_bayes_nest_g, type='pit_ecdf_grouped', group = "Group2", prob = 0.95, plot_diff = F)
# 
# pp_check(model_beta_fish_functional_present_bayes_nest_g , ndraws = 50, type = "scatter_avg_grouped", group = "Group2")
# 
# pairs(model_beta_fish_functional_present_bayes_nest_g, variable = variables(model_beta_fish_functional_present_bayes_nest_g)[1:5])
# fish_functional %>% select(Group1, Group2) %>% unique %>% view

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_fish_functional_present_bayes_nest$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_fish_functional_present_bayes_nest, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_fish_functional_present_bayes_nest, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_fish_functional_present_bayes_nest, type = 'rhat') 


#back transform model summary with point estimates - Now Odds ratios (3:1 odds of bleaching)
model_parameters(model_beta_fish_functional_present_bayes_nest, exponentiate = T, ci = 0.95)



#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####

#Plotting isolation effects
p1 <- fish_functional %>% 
  modelr::data_grid(
    diff_isolation_sc = seq_range(diff_isolation_sc, n = 10),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_present_sst = median(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_functional_present_bayes_nest,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_isolation_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in isolation", y = "Beta diversity")

#
p2 <- fish_functional %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = seq_range(diff_reef_area_sc, n = 10),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_present_sst = median(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_functional_present_bayes_nest,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_reef_area_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in reef area", y = "Beta diversity")

# 
p3 <- fish_functional %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    diff_age_sc = seq_range(diff_age_sc, n = 10),
    dist_sc = mean(dist_sc),
    diff_present_sst = median(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_functional_present_bayes_nest,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_age_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in age", y = "Beta diversity")

# 
p4 <- fish_functional %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    dist_sc = seq_range(dist_sc, n = 10),
    diff_age_sc = mean(diff_age_sc),
    diff_present_sst = median(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_functional_present_bayes_nest,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = dist_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch distance", y = "Beta diversity")

#
p5 <- fish_functional %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    dist_sc = mean(dist_sc),
    diff_age_sc = mean(diff_age_sc),
    diff_present_sst = seq_range(diff_present_sst, n = 10)
  ) %>% 
  add_epred_draws(model_beta_fish_functional_present_bayes_nest,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_present_sst)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch SST", y = "Beta diversity")





#5. Terrestrial - PLANTS BETA TAXONOMIC ----


model_plants_taxonomic <- model_plants_taxonomic %>%
  filter(!Island1 %in% "Farquhar_Atoll",
         !Island2 %in% "Farquhar_Atoll")


# Scaling variables - difference between archipelagos islands
plants_taxonomic <- model_plants_taxonomic %>%
  mutate(
    beta_sne_adj = beta_sne - 0.000001,
    beta_sne_adj = ifelse(beta_sne_adj <= 0, 1e-6, beta_sne_adj),
    d_past_sc = scale(d_past.y), #emerse distance past
    diff_emerse_past_sc = scale(diff_emerse_past), #emerse island area present
    diff_emerse_pres_sc = scale(diff_emerse_pres), #emerse island area past
    diff_island_age_sc = scale(diff_island_age),
    distance_km_sc = scale(Distance_km), #emerse distance present
    diff_present_air = scale(diif_present_air),
    diff_past_air = scale(diff_past_air)
  )

# scaled_diff_emerse_pres <- plants_taxonomic$diff_emerse_pres - mean(plants_taxonomic$diff_emerse_pres)/sd(plants_taxonomic$diff_emerse_pres)
# 
# diff_emerse_pres <- scaled_diff_emerse_pres * sd(plants_taxonomic$diff_emerse_pres) + mean(plants_taxonomic$diff_emerse_pres)
# 
# diff_emerse_pres_sc <-  scale(plants_taxonomic$diff_island_age_sc)
# diff_emerse_pres_sc %>% glimpse


#------------------------- Models

# Null

model_beta_plants_taxonomic_present_nest_null_nest <- brm(
  beta_sne_adj ~ 1 +
    (1 | Archipelago),
  data = plants_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 5000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.98, max_treedepth = 14),
)


model_beta_plants_taxonomic_present_nest_null_nest <- add_criterion(model_beta_plants_taxonomic_present_nest_null_nest, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_plants_taxonomic_present_nest_null_nest,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_plants_taxonomic_present_nest_null_nest.rds"
)

# Assess model validity using LOO
loo(model_beta_plants_taxonomic_present_nest_null_nest)



# Present

model_beta_plants_taxonomic_present_nest <- brm(
  beta_sne_adj ~  diff_emerse_pres_sc + diff_island_age_sc + distance_km_sc + diff_present_air +
    (1 | Archipelago),
  data = plants_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_plants_taxonomic_present_nest <- add_criterion(model_beta_plants_taxonomic_present_nest, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_plants_taxonomic_present_nest,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_plants_taxonomic_present_nest.rds"
)

loo(model_beta_plants_taxonomic_present_nest)

#loo_compare(model_beta_plants_taxonomic_present_nest_null_nest, model_beta_plants_taxonomic_present_nest)

#Look at model summary its on 'logit' scale, because of Beta family link function
model_beta_plants_taxonomic_present_nest

## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_plants_taxonomic_present_nest, ndraws = 50) 
pp_check(model_beta_plants_taxonomic_present_nest, type = "ecdf_overlay", ndraws = 50) 

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_plants_taxonomic_present_nest$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_plants_taxonomic_present_nest, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_plants_taxonomic_present_nest, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_plants_taxonomic_present_nest, type = 'rhat') 


#back transform model summary with point estimates (Ratios)
model_parameters(model_beta_plants_taxonomic_present_nest, exponentiate = T, ci = 0.95)


#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####


#Plotting isolation effects
p1 <- plants_taxonomic %>% 
  modelr::data_grid(distance_km_sc = seq_range(distance_km_sc, n = 100) %>% as.numeric,
                    diff_emerse_pres_sc = median(diff_emerse_pres_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_present_air = median(diff_present_air, na.rm = T)) %>%
  add_epred_draws(model_beta_plants_taxonomic_present_nest, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = distance_km_sc)) +
  # geom_point(data = corals_taxonomic$beta_sne_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Distance between islands", y = "Beta diversity")

p2 <- plants_taxonomic %>% 
  modelr::data_grid(distance_km_sc = median(distance_km_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_pres_sc = seq_range(diff_emerse_pres_sc, n = 10),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_present_air = median(diff_present_air, na.rm = T)) %>%
  add_epred_draws(model_beta_plants_taxonomic_present_nest, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_emerse_pres_sc)) +
  # geom_point(data = corals_taxonomic$beta_sne_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Emersed Island Area", y = "Beta diversity")

p3 <- plants_taxonomic %>% 
  modelr::data_grid(distance_km_sc = median(distance_km_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_pres_sc = median(diff_emerse_pres_sc, na.rm = T),
                    diff_island_age_sc = seq_range(diff_island_age_sc, n = 10),
                    diff_present_air = median(diff_present_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_plants_taxonomic_present_nest, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_island_age_sc)) +
  # geom_point(data = corals_taxonomic$beta_sne_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Island Age", y = "Beta diversity")

p4 <- plants_taxonomic %>% 
  modelr::data_grid(distance_km_sc = median(distance_km_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_pres_sc = median(diff_emerse_pres_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_present_air = seq_range(diff_present_air, n = 10)) %>% 
  add_epred_draws(model_beta_plants_taxonomic_present_nest, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_present_air)) +
  # geom_point(data = corals_taxonomic$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Air Temperature", y = "Beta diversity")

library(cowplot)

p1 / p2 / p3 / p4






#6. Terrestrial - PLANTS BETA FUNCTIONAL ----


model_plants_functional <- model_plants_functional %>%
  filter(!Island1 %in% "Farquhar_Atoll",
         !Island2 %in% "Farquhar_Atoll")



# Scaling variables - difference between archipelagos islands
plants_functional <- model_plants_functional %>%
  mutate(
    beta_sne_adj = beta_nest - 0.000001, 
    beta_sne_adj = ifelse(beta_sne_adj <= 0, 1e-6, beta_sne_adj),
    d_past_sc = scale(d_past.y), #emerse distance past
    diff_emerse_past_sc = scale(diff_emerse_past), #emerse island area present
    diff_emerse_pres_sc = scale(diff_emerse_pres), #emerse island area past
    diff_island_age_sc = scale(diff_island_age),
    distance_km_sc = scale(Distance_km), #emerse distance present
    diff_present_air = scale(diif_present_air),
    diff_past_air = scale(diff_past_air)
  )


#------------------------- Models

# Null

model_beta_plants_functional_present_null_nest <- brm(
  beta_sne_adj ~ 1 +
    (1 | Archipelago),
  data = plants_functional,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 5000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.98, max_treedepth = 14),
)


model_beta_plants_functional_present_null_nest <- add_criterion(model_beta_plants_functional_present_null_nest, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_plants_functional_present_null_nest,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_plants_functional_present_null_nest.rds"
)

# Assess model validity using LOO
loo(model_beta_plants_functional_present_null_nest)



# Present

model_beta_plants_functional_present_nest <- brm(
  beta_sne_adj ~  diff_emerse_pres_sc + diff_island_age_sc + distance_km_sc + diff_present_air +
    (1 | Archipelago),
  data = plants_functional,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_plants_functional_present_nest <- add_criterion(model_beta_plants_functional_present_nest, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_plants_functional_present_nest,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_plants_functional_present_nest.rds"
)

loo(model_beta_plants_functional_present_nest)

#loo_compare(model_beta_plants_functional_present_null_nest, model_beta_plants_functional_present_nest)

#Look at model summary its on 'logit' scale, because of Beta family link function
model_beta_plants_functional_present_nest

## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_plants_functional_present_nest, ndraws = 50) 
pp_check(model_beta_plants_functional_present_nest, type = "ecdf_overlay", ndraws = 50) 

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_plants_functional_present_nest$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_plants_functional_present_nest, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_plants_functional_present_nest, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_plants_functional_present_nest, type = 'rhat') 


#back transform model summary with point estimates (Ratios)
model_parameters(model_beta_plants_functional_present_nest, exponentiate = T, ci = 0.95)


#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####


#Plotting isolation effects
p1 <- plants_functional %>% 
  modelr::data_grid(distance_km_sc = seq_range(distance_km_sc, n = 100) %>% as.numeric,
                    diff_emerse_pres_sc = median(diff_emerse_pres_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_present_air = median(diff_present_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_plants_functional_present_nest, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = distance_km_sc)) +
  # geom_point(data = corals_functional$beta_sne_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Distance between islands", y = "Beta diversity")

p2 <- plants_functional %>% 
  modelr::data_grid(distance_km_sc = median(distance_km_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_pres_sc = seq_range(diff_emerse_pres_sc, n = 10),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_present_air = median(diff_present_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_plants_functional_present_nest, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_emerse_pres_sc)) +
  # geom_point(data = corals_functional$beta_sne_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Emersed Island Area", y = "Beta diversity")

p3 <- plants_functional %>% 
  modelr::data_grid(distance_km_sc = median(distance_km_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_pres_sc = median(diff_emerse_pres_sc, na.rm = T),
                    diff_island_age_sc = seq_range(diff_island_age_sc, n = 10),
                    diff_present_air = median(diff_present_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_plants_functional_present_nest, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_island_age_sc)) +
  # geom_point(data = corals_functional$beta_sne_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Island Age", y = "Beta diversity")

p4 <- plants_functional %>% 
  modelr::data_grid(distance_km_sc = median(distance_km_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_pres_sc = median(diff_emerse_pres_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_present_air = seq_range(diff_present_air, n = 10)) %>% 
  add_epred_draws(model_beta_plants_functional_present_nest, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_present_air)) +
  # geom_point(data = corals_taxonomic$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Air Temperature", y = "Beta diversity")

library(cowplot)

p1 / p2 / p3 / p4





#7. Terrestrial - BIRDS BETA TAXONOMIC ----


model_birds_taxonomic <- model_birds_taxonomic %>%
  filter(!Island1 %in% "Farquhar_Atoll",
         !Island2 %in% "Farquhar_Atoll")



# Scaling variables - difference between archipelagos islands
birds_taxonomic <- model_birds_taxonomic %>%
  mutate(
    beta_sne_adj = beta_sne - 0.000001,
    beta_sne_adj = ifelse(beta_sne_adj <= 0, 1e-6, beta_sne_adj),
    d_past_sc = scale(d_past.y), #emerse distance past
    diff_emerse_past_sc = scale(diff_emerse_past), #emerse island area present
    diff_emerse_pres_sc = scale(diff_emerse_pres), #emerse island area past
    diff_island_age_sc = scale(diff_island_age),
    distance_km_sc = scale(Distance_km), #emerse distance present
    diff_present_air = scale(diif_present_air),
    diff_past_air = scale(diff_past_air)
  )




#------------------------- Models

# Null

model_beta_birds_taxonomic_present_null_nest <- brm(
  beta_sne_adj ~ 1 +
    (1 | Archipelago),
  data = birds_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 5000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.98, max_treedepth = 14),
)


model_beta_birds_taxonomic_present_null_nest <- add_criterion(model_beta_birds_taxonomic_present_null_nest, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_birds_taxonomic_present_null_nest,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_birds_taxonomic_present_null_nest.rds"
)

# Assess model validity using LOO
loo(model_beta_birds_taxonomic_present_null_nest)



# Present

model_beta_birds_taxonomic_present_nest <- brm(
  beta_sne_adj ~  diff_emerse_pres_sc + diff_island_age_sc + distance_km_sc + diff_present_air +
    (1 | Archipelago),
  data = birds_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_birds_taxonomic_present_nest <- add_criterion(model_beta_birds_taxonomic_present_nest, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_birds_taxonomic_present_nest,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_birds_taxonomic_present_nest.rds"
)

loo(model_beta_birds_taxonomic_present_nest)

#loo_compare(model_beta_birds_taxonomic_present_null_nest, model_beta_birds_taxonomic_present_nest)

 #Look at model summary its on 'logit' scale, because of Beta family link function
model_beta_birds_taxonomic_present_nest


## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_birds_taxonomic_present_nest, ndraws = 50) 
pp_check(model_beta_birds_taxonomic_present_nest, type = "ecdf_overlay", ndraws = 50) 

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_birds_taxonomic_present_nest$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_birds_taxonomic_present_nest, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_birds_taxonomic_present_nest, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_birds_taxonomic_present_nest, type = 'rhat') 


#back transform model summary with point estimates (Ratios)
model_parameters(model_beta_birds_taxonomic_present_nest, exponentiate = T, ci = 0.95)


#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####


#Plotting isolation effects
p1 <- birds_taxonomic %>% 
  modelr::data_grid(distance_km_sc = seq_range(distance_km_sc, n = 100) %>% as.numeric,
                    diff_emerse_pres_sc = median(diff_emerse_pres_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_present_air = median(diff_present_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_birds_taxonomic_present_nest, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = distance_km_sc)) +
  # geom_point(data = corals_taxonomic$beta_sne_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Distance between islands", y = "Beta diversity")

p2 <- birds_taxonomic %>% 
  modelr::data_grid(distance_km_sc = median(distance_km_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_pres_sc = seq_range(diff_emerse_pres_sc, n = 10),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_present_air = median(diff_present_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_birds_taxonomic_present_nest, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_emerse_pres_sc)) +
  # geom_point(data = corals_taxonomic$beta_sne_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Emersed Island Area", y = "Beta diversity")

p3 <- birds_taxonomic %>% 
  modelr::data_grid(distance_km_sc = median(distance_km_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_pres_sc = median(diff_emerse_pres_sc, na.rm = T),
                    diff_island_age_sc = seq_range(diff_island_age_sc, n = 10),
                    diff_present_air = median(diff_present_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_birds_taxonomic_present_nest, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_island_age_sc)) +
  # geom_point(data = corals_taxonomic$beta_sne_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Island Age", y = "Beta diversity")

p4 <- birds_taxonomic %>% 
  modelr::data_grid(distance_km_sc = median(distance_km_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_pres_sc = median(diff_emerse_pres_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_present_air = seq_range(diff_present_air, n = 10)) %>% 
  add_epred_draws(model_beta_birds_taxonomic_present_nest, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_present_air)) +
  # geom_point(data = corals_taxonomic$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Air Temperature", y = "Beta diversity")

library(cowplot)

p1 / p2 / p3 / p4




#8. Terrestrial - BIRDS BETA FUNCTIONAL ----


model_birds_functional <- model_birds_functional %>%
  filter(!Island1 %in% "Farquhar_Atoll",
         !Island2 %in% "Farquhar_Atoll")



# Scaling variables - difference between archipelagos islands
birds_functional <- model_birds_functional %>%
  mutate(
    beta_sne_adj = beta_nest - 0.000001,
    beta_sne_adj = ifelse(beta_sne_adj <= 0, 1e-6, beta_sne_adj),
    d_past_sc = scale(d_past.y), #emerse distance past
    diff_emerse_past_sc = scale(diff_emerse_past), #emerse island area present
    diff_emerse_pres_sc = scale(diff_emerse_pres), #emerse island area past
    diff_island_age_sc = scale(diff_island_age),
    distance_km_sc = scale(Distance_km), #emerse distance present
    diff_present_air = scale(diif_present_air),
    diff_past_air = scale(diff_past_air)
  )




#------------------------- Models

# Null

model_beta_birds_functional_present_null_nest <- brm(
  beta_sne_adj ~ 1 +
    (1 | Archipelago),
  data = birds_functional,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 5000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.98, max_treedepth = 14),
)


model_beta_birds_functional_present_null_nest <- add_criterion(model_beta_birds_functional_present_null_nest, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_birds_functional_present_null_nest,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_birds_functional_present_null_nest.rds"
)


# Assess model validity using LOO
loo(model_beta_birds_functional_present_null_nest)



# Present

model_beta_birds_functional_present_nest <- brm(
  beta_sne_adj ~  diff_emerse_pres_sc + diff_island_age_sc + distance_km_sc + diff_present_air +
    (1 | Archipelago),
  data = birds_functional,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_birds_functional_present_nest <- add_criterion(model_beta_birds_functional_present_nest, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_birds_functional_present_nest,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_birds_functional_present_nest.rds"
)


loo(model_beta_birds_functional_present_nest)

#loo_compare(model_beta_birds_functional_present_nest, model_beta_birds_functional_present_null_nest)

#Look at model summary its on 'logit' scale, because of Beta family link function
model_beta_birds_functional_present_nest

## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_birds_functional_present_nest, ndraws = 50) 
pp_check(model_beta_birds_functional_present_nest, type = "ecdf_overlay", ndraws = 50) 

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_birds_functional_present_nest$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_birds_functional_present_nest, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_birds_functional_present_nest, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_birds_functional_present_nest, type = 'rhat') 


#back transform model summary with point estimates (Ratios)
model_parameters(model_beta_birds_functional_present_nest, exponentiate = T, ci = 0.95)


#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####


#Plotting isolation effects
p1 <- birds_functional %>% 
  modelr::data_grid(distance_km_sc = seq_range(distance_km_sc, n = 100) %>% as.numeric,
                    diff_emerse_pres_sc = median(diff_emerse_pres_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_present_air = median(diff_present_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_birds_functional_present_nest, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = distance_km_sc)) +
  # geom_point(data = corals_functional$beta_sne_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Distance between islands", y = "Beta diversity")

p2 <- birds_functional %>% 
  modelr::data_grid(distance_km_sc = median(distance_km_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_pres_sc = seq_range(diff_emerse_pres_sc, n = 10),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_present_air = median(diff_present_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_birds_functional_present_nest, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_emerse_pres_sc)) +
  # geom_point(data = corals_functional$beta_sne_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Emersed Island Area", y = "Beta diversity")

p3 <- birds_functional %>% 
  modelr::data_grid(distance_km_sc = median(distance_km_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_pres_sc = median(diff_emerse_pres_sc, na.rm = T),
                    diff_island_age_sc = seq_range(diff_island_age_sc, n = 10),
                    diff_present_air = median(diff_present_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_birds_functional_present_nest, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_island_age_sc)) +
  # geom_point(data = corals_functional$beta_sne_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Island Age", y = "Beta diversity")

p4 <- birds_functional %>% 
  modelr::data_grid(distance_km_sc = median(distance_km_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_pres_sc = median(diff_emerse_pres_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_present_air = seq_range(diff_present_air, n = 10)) %>% 
  add_epred_draws(model_beta_birds_functional_present_nest, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_present_air)) +
  # geom_point(data = corals_taxonomic$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Air Temperature", y = "Beta diversity")

library(cowplot)

p1 / p2 / p3 / p4




# NESTEDNESS PAST ----

#-------------------------------------------------------------------------------


#1. Marine - CORALS BETA TAXONOMIC ----


## Scaling variables - difference between marine regions ----
corals_taxonomic <- model_corals_taxonomic %>%
  mutate(
    beta_sne_adj = beta_sne - 0.000001,
    beta_sne_adj = ifelse(beta_sne_adj <= 0, 1e-6, beta_sne_adj),
    diff_isolation_sc = scale(diff_isolation)[,1],
    diff_reef_area_sc = scale(diff_reef_area)[,1],
    diff_past_isolation_sc = scale(diff_past_isolation)[,1], #quaternary distance
    diff_past_reef_area_sc = scale(diff_past_reef_area)[,1],
    diff_age_sc = scale(diff_age)[,1],
    dist_sc = scale(Distance)[,1],
    diff_present_sst = scale(diff_sst)[,1],
    diff_past_sst = scale(diff_past_sst)[,1]
  )

corals_taxonomic %>% glimpse

#------------------------- Models

# past


# Bayes models
# Null

# model_beta_corals_taxonomic_past_null_nest_nest_gp3 <- brm(
#   beta_sne_adj ~ 1 + (1|Group1) + (1|Group2) + (1|key),
#   data = corals_taxonomic,
#   family = beta_family(),
#   chains = 4,
#   cores = 4,
#   iter = 4000,
#   warmup = 3000,
#   control = list(adapt_delta = 0.99, max_treedepth = 14),
#   save_pars = save_pars(all = T)
# )

# model_beta_corals_taxonomic_past_null_nest_nest_gp3 <- add_criterion(model_beta_corals_taxonomic_past_null_nest_nest_gp3, "loo", save_psis = TRUE, reloo = T)


# model_beta_corals_taxonomic_past_null_nest <- brm(
#   beta_sne_adj ~ 1,
#   data = corals_taxonomic,
#   family = beta_family(),
#   chains = 4,
#   cores = 4,
#   iter = 4000,
#   warmup = 3000,
#   control = list(adapt_delta = 0.99, max_treedepth = 14),
#   save_pars = save_pars(all = T)
# )

# model_beta_corals_taxonomic_past_null_nest <- add_criterion(model_beta_corals_taxonomic_past_null_nest, "loo", save_psis = TRUE, reloo = T)
# 
# loo_compare(model_beta_corals_taxonomic_past_null_nest, model_beta_corals_taxonomic_past_null_nest_nest_gp2, model_beta_corals_taxonomic_past_null_nest_nest_gp1, model_beta_corals_taxonomic_past_null_nest_nest_gp3)

model_beta_corals_taxonomic_past_null_nest_gp2 <- brm(
  beta_sne_adj ~ 1 + 
    (1+dist_sc|Group1) + (1+dist_sc|Group2),
  data = corals_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  control = list(adapt_delta = 0.99, max_treedepth = 14),
  save_pars = save_pars(all = T)
)

model_beta_corals_taxonomic_past_null_nest_gp2 <- add_criterion(model_beta_corals_taxonomic_past_null_nest_gp2, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_corals_taxonomic_past_null_nest_gp2,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_corals_taxonomic_past_null_nest_gp2.rds"
)


# Assess model validity using LOO
loo(model_beta_corals_taxonomic_past_null_nest_gp2)

# past
model_beta_corals_taxonomic_past_bayes_nest <- brm(
  beta_sne_adj ~ diff_past_isolation_sc + diff_past_reef_area_sc + diff_age_sc + dist_sc + diff_past_sst + 
    (1+dist_sc|Group1) + (1+dist_sc|Group2),
  data = corals_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_corals_taxonomic_past_bayes_nest <- add_criterion(model_beta_corals_taxonomic_past_bayes_nest, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_corals_taxonomic_past_bayes_nest,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_corals_taxonomic_past_bayes_nest.rds"
)

#Look at model summary its on 'logit' scale, because of bernouli family link function
model_beta_corals_taxonomic_past_bayes_nest

# Assess model validity using LOO
loo(model_beta_corals_taxonomic_past_bayes_nest)


#loo_compare(model_beta_corals_taxonomic_past_bayes_nest, model_beta_corals_taxonomic_past_null_nest_gp2)
## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_corals_taxonomic_past_bayes_nest , ndraws = 50) 
pp_check(model_beta_corals_taxonomic_past_bayes_nest, type = "ecdf_overlay", ndraws = 50) 

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_corals_taxonomic_past_bayes_nest$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_corals_taxonomic_past_bayes_nest, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_corals_taxonomic_past_bayes_nest, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_corals_taxonomic_past_bayes_nest, type = 'rhat') 


#back transform model summary with point estimates - Now Odds ratios (3:1 odds of bleaching)
model_parameters(model_beta_corals_taxonomic_past_bayes_nest, exponentiate = T, ci = 0.95)



#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####

#Plotting isolation effects
# 
p1 <- corals_taxonomic %>% 
  modelr::data_grid(
    diff_past_isolation_sc = seq_range(diff_past_isolation_sc, n = 10),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_past_sst = mean(diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_taxonomic_past_bayes_nest,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_past_isolation_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in isolation", y = "Beta diversity")

#
p2 <- corals_taxonomic %>% 
  modelr::data_grid(
    diff_past_isolation_sc = seq_range(diff_past_isolation_sc, n=10),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_past_sst = mean(diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_taxonomic_past_bayes_nest,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_past_isolation_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in area", y = "Beta diversity")

# 
p3 <- corals_taxonomic %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_past_isolation_sc),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    diff_age_sc = seq_range(diff_age_sc, n = 10),
    dist_sc = mean(dist_sc),
    diff_past_sst = mean(diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_taxonomic_past_bayes_nest,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_age_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in age", y = "Beta diversity")

# 
p4 <- corals_taxonomic %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_past_isolation_sc),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    dist_sc = seq_range(dist_sc, n = 10),
    diff_age_sc = mean(diff_age_sc),
    diff_past_sst = mean(diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_taxonomic_past_bayes_nest,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = dist_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch distance", y = "Beta diversity")

#
p5 <- corals_taxonomic %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_past_isolation_sc),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    dist_sc = mean(dist_sc),
    diff_age_sc = mean(diff_age_sc),
    diff_past_sst = seq_range(diff_past_sst, n = 10)
  ) %>% 
  add_epred_draws(model_beta_corals_taxonomic_past_bayes_nest,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_past_sst)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch SST", y = "Beta diversity")


p1 / p2 / p3 / p4 / p5



#2. Marine - CORALS BETA FUNCTIONAL ----


## Scaling variables - difference between marine regions ----
corals_functional <- model_corals_functional %>%
  mutate(
    beta_sne_adj = beta_sne - 0.000001,  
    diff_isolation_sc = scale(diff_isolation)[,1],
    diff_reef_area_sc = scale(diff_reef_area)[,1],
    diff_past_isolation_sc = scale(diff_past_isolation)[,1], #quaternary distance
    diff_past_reef_area_sc = scale(diff_past_reef_area)[,1],
    diff_age_sc = scale(diff_age)[,1],
    dist_sc = scale(Distance)[,1],
    diff_present_sst = scale(diff_sst)[,1],
    diff_past_sst = scale(diff_past_sst)[,1]
  )

corals_functional %>% glimpse

#------------------------- Models


# Bayes models
# Null

model_beta_corals_functional_past_null_nest_gp2 <- brm(
  beta_sne_adj ~ 1 + 
    (1+dist_sc|Group1) + (1+dist_sc|Group2),
  data = corals_functional,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  control = list(adapt_delta = 0.99, max_treedepth = 14),
  save_pars = save_pars(all = T)
)

model_beta_corals_functional_past_null_nest_gp2 <- add_criterion(model_beta_corals_functional_past_null_nest_gp2, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_corals_functional_past_null_nest_gp2,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_corals_functional_past_null_nest_gp2.rds"
)


# Assess model validity using LOO
loo(model_beta_corals_functional_past_null_nest_gp2)

# past
model_beta_corals_functional_past_bayes_nest <- brm(
  beta_sne_adj ~ diff_past_isolation_sc + diff_past_reef_area_sc + diff_age_sc + dist_sc + diff_past_sst +
  (1+dist_sc|Group1) + (1+dist_sc|Group2),
  data = corals_functional,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_corals_functional_past_bayes_nest <- add_criterion(model_beta_corals_functional_past_bayes_nest, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_corals_functional_past_bayes_nest,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_corals_functional_past_bayes_nest.rds"
)

# Assess model validity using LOO
loo(model_beta_corals_functional_past_bayes_nest)

#Look at model summary its on 'logit' scale, because of bernouli family link function
model_beta_corals_functional_past_bayes_nest


#loo_compare(model_beta_corals_functional_past_bayes_nest, model_beta_corals_functional_past_null_nest_gp2)

## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_corals_functional_past_null_nest_gp2 , ndraws = 50) 
pp_check(model_beta_corals_functional_past_bayes_nest, type = "ecdf_overlay", ndraws = 50) 

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_corals_functional_past_bayes_nest$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_corals_functional_past_bayes_nest, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_corals_functional_past_bayes_nest, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_corals_functional_past_bayes_nest, type = 'rhat') 


#back transform model summary with point estimates - Now Odds ratios (3:1 odds of bleaching)
model_parameters(model_beta_corals_functional_past_bayes_nest, exponentiate = T, ci = 0.95)


#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####

#Plotting isolation effects
p1 <- corals_functional %>% 
  modelr::data_grid(
    diff_past_isolation_sc = seq_range(diff_past_isolation_sc, n = 10),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_past_sst = mean (diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_functional_past_bayes_nest,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_past_isolation_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in isolation", y = "Beta diversity")

#
p2 <- corals_functional %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_past_isolation_sc),
    diff_past_reef_area_sc = seq_range(diff_past_reef_area_sc, n = 10),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_past_sst = mean (diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_functional_past_bayes_nest,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_past_reef_area_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in reef area", y = "Beta diversity")

# 
p3 <- corals_functional %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_past_isolation_sc),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    diff_age_sc = seq_range(diff_age_sc, n = 10),
    dist_sc = mean(dist_sc),
    diff_past_sst = mean (diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_functional_past_bayes_nest,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_age_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in age", y = "Beta diversity")

# 
p4 <- corals_functional %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_past_isolation_sc),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = seq_range(dist_sc, n = 10),
    diff_past_sst = mean (diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_functional_past_bayes_nest,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = dist_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch distance", y = "Beta diversity")

#  
p5 <- corals_taxonomic %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_past_isolation_sc),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    dist_sc = mean(dist_sc),
    diff_age_sc = mean(diff_age_sc),
    diff_past_sst = seq_range(diff_past_sst, n = 10)
  ) %>% 
  add_epred_draws(model_beta_corals_functional_past_bayes_nest,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_past_sst)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch SST", y = "Beta diversity")

p1 / p2 / p3 / p4 / p5


#3. Marine - FISH BETA TAXONOMIC ----


## Scaling variables - difference between marine regions ----
fish_taxonomic <- model_fish_taxonomic %>%
  mutate(
    beta_sne_adj = beta_sne - 0.000001,
    beta_sne_adj = ifelse(beta_sne_adj <= 0, 1e-6, beta_sne_adj),
    diff_isolation_sc = scale(diff_isolation)[,1],
    diff_reef_area_sc = scale(diff_reef_area)[,1],
    diff_past_isolation_sc = scale(diff_past_isolation)[,1], #quaternary distance
    diff_past_reef_area_sc = scale(diff_past_reef_area)[,1],
    diff_age_sc = scale(diff_age)[,1],
    dist_sc = scale(Distance)[,1],
    diff_present_sst = scale(diff_sst)[,1],
    diff_past_sst = scale(diff_past_sst)[,1]
  )

fish_taxonomic %>% glimpse



#------------------------- Models


# Bayes models
# Null


# model_beta_fish_taxonomic_past_null_nest_nest_gp3 <- add_criterion(model_beta_fish_taxonomic_past_null_nest_nest_gp3, "loo", save_psis = TRUE, reloo = T)
# model_beta_fish_taxonomic_past_null_nest_nest_gp2 <- brm(
#   beta_sne_adj ~ 1 + (1|Group1) + (1|Group2) ,
#   data = fish_taxonomic,
#   family = beta_family(),
#   chains = 4,
#   cores = 4,
#   iter = 4000,
#   warmup = 3000,
#   control = list(adapt_delta = 0.99, max_treedepth = 14),
#   save_pars = save_pars(all = T)
# )

model_beta_fish_taxonomic_past_null_nest_gp2 <- brm(
  beta_sne_adj ~ 1 + (1+dist_sc|Group1) + (1+dist_sc|Group2),
  data = fish_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  control = list(adapt_delta = 0.99, max_treedepth = 14),
  save_pars = save_pars(all = T)
)

model_beta_fish_taxonomic_past_null_nest_gp2 <- add_criterion(model_beta_fish_taxonomic_past_null_nest_gp2, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_fish_taxonomic_past_null_nest_gp2,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_fish_taxonomic_past_null_nest_gp2.rds"
)

# Assess model validity using LOO
loo(model_beta_fish_taxonomic_past_null_nest_gp2)

#Look at model summary its on 'logit' scale, because of bernouli family link function
model_beta_fish_taxonomic_past_null_nest_gp2


# past
model_beta_fish_taxonomic_past_bayes_nest <- brm(
  beta_sne_adj ~ diff_past_isolation_sc + diff_past_reef_area_sc + diff_age_sc + dist_sc + diff_past_sst +
  (1+dist_sc|Group1) + (1+dist_sc|Group2),
  data = fish_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_fish_taxonomic_past_bayes_nest <- add_criterion(model_beta_fish_taxonomic_past_bayes_nest, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_fish_taxonomic_past_bayes_nest,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_fish_taxonomic_past_bayes_nest.rds"
)

model_beta_fish_taxonomic_past_bayes_nest

# Assess model validity using LOO
loo(model_beta_fish_taxonomic_past_bayes_nest)


#loo_compare(model_beta_fish_taxonomic_past_bayes_nest, model_beta_fish_taxonomic_past_null_nest_gp2)
## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_fish_taxonomic_past_bayes_nest , ndraws = 50) 
pp_check(model_beta_fish_taxonomic_past_bayes_nest, type = "ecdf_overlay", ndraws = 50) 

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_fish_taxonomic_past_bayes_nest$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_fish_taxonomic_past_bayes_nest, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_fish_taxonomic_past_bayes_nest, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_fish_taxonomic_past_bayes_nest, type = 'rhat') 


#back transform model summary with point estimates - Now Odds ratios (3:1 odds of bleaching)
model_parameters(model_beta_fish_taxonomic_past_bayes_nest, exponentiate = T, ci = 0.95)


#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####

#Plotting isolation effects
p1 <- fish_taxonomic %>% 
  modelr::data_grid(
    diff_past_isolation_sc = seq_range(diff_past_isolation_sc, n = 10),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_past_sst = mean(diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_taxonomic_past_bayes_nest,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_past_isolation_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in isolation", y = "Beta diversity")

#
p2 <- fish_taxonomic %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_past_isolation_sc),
    diff_past_reef_area_sc = seq_range(diff_past_reef_area_sc, n = 10),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_past_sst = mean(diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_taxonomic_past_bayes_nest,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_past_reef_area_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in reef area", y = "Beta diversity")

# 
p3 <- fish_taxonomic %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_past_isolation_sc),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    diff_age_sc = seq_range(diff_age_sc, n = 10),
    dist_sc = mean(dist_sc),
    diff_past_sst = mean(diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_taxonomic_past_bayes_nest,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_age_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in age", y = "Beta diversity")

# 
p4 <- fish_taxonomic %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_past_isolation_sc),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    dist_sc = seq_range(dist_sc, n = 10),
    diff_age_sc = mean(diff_age_sc),
    diff_past_sst = mean(diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_taxonomic_past_bayes_nest,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = dist_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch distance", y = "Beta diversity")

#
p5 <- fish_taxonomic %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_past_isolation_sc),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    dist_sc = mean(dist_sc),
    diff_age_sc = mean(diff_age_sc),
    diff_past_sst = seq_range(diff_past_sst, n = 10)
  ) %>% 
  add_epred_draws(model_beta_fish_taxonomic_past_bayes_nest,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_past_sst)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch SST", y = "Beta diversity")

p1 / p2 / p3 / p4 / p5



#4. Marine - FISH BETA FUNCTIONAL ----


## Scaling variables - difference between marine regions ----
fish_functional <- model_fish_functional %>%
  mutate(
    beta_sne_adj = beta_sne - 0.000001,  
    beta_sne_adj = ifelse(beta_sne_adj <= 0, 1e-6, beta_sne_adj),
    diff_isolation_sc = scale(diff_isolation)[,1],
    diff_reef_area_sc = scale(diff_reef_area)[,1],
    diff_past_isolation_sc = scale(diff_past_isolation)[,1], #quaternary distance
    diff_past_reef_area_sc = scale(diff_past_reef_area)[,1],
    diff_age_sc = scale(diff_age)[,1],
    dist_sc = scale(Distance)[,1],
    diff_present_sst = scale(diff_sst)[,1],
    diff_past_sst = scale(diff_past_sst)[,1]
  ) 

fish_functional %>% glimpse

#------------------------- Models


# Bayes models
# Null
model_beta_fish_functional_past_null_nest_nest_gp2 <- brm(
  beta_sne_adj ~ 1 + (1+dist_sc|Group1) + (1+dist_sc|Group2) ,
  data = fish_functional,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  control = list(adapt_delta = 0.99, max_treedepth = 14),
  save_pars = save_pars(all = T)
)

model_beta_fish_functional_past_null_nest_nest_gp2 <- add_criterion(model_beta_fish_functional_past_null_nest_nest_gp2, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_fish_functional_past_null_nest_nest_gp2,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_fish_functional_past_null_nest_nest_gp2.rds"
)

# Assess model validity using LOO
loo(model_beta_fish_functional_past_null_nest_nest_gp2)


# past
model_beta_fish_functional_past_bayes_nest <- brm(
  beta_sne_adj ~  1 + diff_past_isolation_sc + diff_past_reef_area_sc + diff_age_sc + dist_sc + diff_past_sst + 
  (1+dist_sc|Group1) + (1+dist_sc|Group2),
  data = fish_functional,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_fish_functional_past_bayes_nest <- add_criterion(model_beta_fish_functional_past_bayes_nest, "loo", save_psis = TRUE, reloo = T, moment_match = TRUE)


saveRDS(
  model_beta_fish_functional_past_bayes_nest,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_fish_functional_past_bayes_nest.rds"
)


# Assess model validity using LOO
loo(model_beta_fish_functional_past_bayes_nest)


#loo_compare(model_beta_fish_functional_past_bayes_nest, model_beta_fish_functional_past_null_nest_nest_gp2)

## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_fish_functional_past_bayes_nest , ndraws = 50) 
pp_check(model_beta_fish_functional_past_bayes_nest, type = "ecdf_overlay", ndraws = 50) 
pp_check(model_beta_fish_functional_past_bayes_nest , ndraws = 50, type = "dens_overlay_grouped", group = "Group1")
pp_check(model_beta_fish_functional_past_bayes_nest , ndraws = 50, type = "scatter_avg_grouped", group = "Group1")
pp_check(model_beta_fish_functional_past_bayes_nest , ndraws = 50, type = "error_scatter_avg_grouped", group = "Group1")
pp_check(model_beta_fish_functional_past_bayes_nest, type='pit_ecdf_grouped', group = "Group1", prob = 0.95, plot_diff = F)
pp_check(model_beta_fish_functional_past_bayes_nest, type='pit_ecdf_grouped', group = "Group2", prob = 0.95, plot_diff = F)

pp_check(model_beta_fish_functional_past_bayes_nest , ndraws = 50, type = "scatter_avg_grouped", group = "Group2")

pairs(model_beta_fish_functional_past_bayes_nest, variable = variables(model_beta_fish_functional_past_bayes_nest)[1:5])
fish_functional %>% select(Group1, Group2) %>% unique %>% view

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_fish_functional_past_bayes_nest$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_fish_functional_past_bayes_nest, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_fish_functional_past_bayes_nest, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_fish_functional_past_bayes_nest, type = 'rhat') 


#back transform model summary with point estimates - Now Odds ratios (3:1 odds of bleaching)
model_parameters(model_beta_fish_functional_past_bayes_nest, exponentiate = T, ci = 0.95)



#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####

#Plotting isolation effects
p1 <- fish_functional %>% 
  modelr::data_grid(
    diff_past_isolation_sc = seq_range(diff_past_isolation_sc, n = 10),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_past_sst = median(diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_functional_past_bayes_nest,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_past_isolation_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in isolation", y = "Beta diversity")

#
p2 <- fish_functional %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_past_isolation_sc),
    diff_past_reef_area_sc = seq_range(diff_past_reef_area_sc, n = 10),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_past_sst = median(diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_functional_past_bayes_nest,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_past_reef_area_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in reef area", y = "Beta diversity")

# 
p3 <- fish_functional %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_past_isolation_sc),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    diff_age_sc = seq_range(diff_age_sc, n = 10),
    dist_sc = mean(dist_sc),
    diff_past_sst = median(diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_functional_past_bayes_nest,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_age_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in age", y = "Beta diversity")

# 
p4 <- fish_functional %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_past_isolation_sc),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    dist_sc = seq_range(dist_sc, n = 10),
    diff_age_sc = mean(diff_age_sc),
    diff_past_sst = median(diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_functional_past_bayes_nest,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = dist_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch distance", y = "Beta diversity")

#
p5 <- fish_functional %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_past_isolation_sc),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    dist_sc = mean(dist_sc),
    diff_age_sc = mean(diff_age_sc),
    diff_past_sst = seq_range(diff_past_sst, n = 10)
  ) %>% 
  add_epred_draws(model_beta_fish_functional_past_bayes_nest,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_past_sst)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch SST", y = "Beta diversity")

p1 / p2 / p3 / p4 / p5





#5. Terrestrial - PLANTS BETA TAXONOMIC ----


model_plants_taxonomic <- model_plants_taxonomic %>%
  filter(!Island1 %in% "Farquhar_Atoll",
         !Island2 %in% "Farquhar_Atoll")


# Scaling variables - difference between archipelagos islands
plants_taxonomic <- model_plants_taxonomic %>%
  mutate(
    beta_sne_adj = beta_sne - 0.000001,
    beta_sne_adj = ifelse(beta_sne_adj <= 0, 1e-6, beta_sne_adj),
    d_past_sc = scale(d_past.y), #emerse distance past
    diff_emerse_past_sc = scale(diff_emerse_past), #emerse island area past
    diff_emerse_pres_sc = scale(diff_emerse_pres), #emerse island area past
    diff_island_age_sc = scale(diff_island_age),
    distance_km_sc = scale(Distance_km), #emerse distance past
    diff_present_air = scale(diif_present_air),
    diff_past_air = scale(diff_past_air)
  )

# scaled_diff_emerse_pres <- plants_taxonomic$diff_emerse_pres - mean(plants_taxonomic$diff_emerse_pres)/sd(plants_taxonomic$diff_emerse_pres)
# 
# diff_emerse_pres <- scaled_diff_emerse_pres * sd(plants_taxonomic$diff_emerse_pres) + mean(plants_taxonomic$diff_emerse_pres)
# 
# diff_emerse_pres_sc <-  scale(plants_taxonomic$diff_island_age_sc)
# diff_emerse_pres_sc %>% glimpse


#------------------------- Models

# Null

model_beta_plants_taxonomic_past_null_nest <- brm(
  beta_sne_adj ~ 1 +
    (1 | Archipelago),
  data = plants_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 5000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.98, max_treedepth = 14),
)


model_beta_plants_taxonomic_past_null_nest <- add_criterion(model_beta_plants_taxonomic_past_null_nest, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_plants_taxonomic_past_null_nest,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_plants_taxonomic_past_null_nest.rds"
)

# Assess model validity using LOO
loo(model_beta_plants_taxonomic_past_null_nest)



# past

model_beta_plants_taxonomic_past_nest <- brm(
  beta_sne_adj ~  diff_emerse_past_sc + diff_island_age_sc + d_past_sc + diff_past_air +
    (1 | Archipelago),
  data = plants_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_plants_taxonomic_past_nest <- add_criterion(model_beta_plants_taxonomic_past_nest, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_plants_taxonomic_past_nest,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_plants_taxonomic_past_nest.rds"
)


loo(model_beta_plants_taxonomic_past_nest)

#loo_compare(model_beta_plants_taxonomic_past_nest, model_beta_plants_taxonomic_past_null_nest)

#Look at model summary its on 'logit' scale, because of Beta family link function
model_beta_plants_taxonomic_past_nest

## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_plants_taxonomic_past_nest, ndraws = 50) 
pp_check(model_beta_plants_taxonomic_past_nest, type = "ecdf_overlay", ndraws = 50) 

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_plants_taxonomic_past$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_plants_taxonomic_past_nest, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_plants_taxonomic_past_nest, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_plants_taxonomic_past_nest, type = 'rhat') 


#back transform model summary with point estimates (Ratios)
model_parameters(model_beta_plants_taxonomic_past_nest, exponentiate = T, ci = 0.95)


#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####


#Plotting isolation effects
p1 <- plants_taxonomic %>% 
  modelr::data_grid(d_past_sc = seq_range(d_past_sc, n = 100) %>% as.numeric,
                    diff_emerse_past_sc = median(diff_emerse_past_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_past_air = median(diff_past_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_plants_taxonomic_past_nest, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = d_past_sc)) +
  # geom_point(data = corals_taxonomic$beta_sne_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Distance between islands", y = "Beta diversity")

p2 <- plants_taxonomic %>% 
  modelr::data_grid(d_past_sc = median(d_past_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_past_sc = seq_range(diff_emerse_past_sc, n = 10),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_past_air = median(diff_past_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_plants_taxonomic_past_nest, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_emerse_past_sc)) +
  # geom_point(data = corals_taxonomic$beta_sne_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Emersed Island Area", y = "Beta diversity")

p3 <- plants_taxonomic %>% 
  modelr::data_grid(d_past_sc = median(d_past_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_past_sc = median(diff_emerse_past_sc, na.rm = T),
                    diff_island_age_sc = seq_range(diff_island_age_sc, n = 10),
                    diff_past_air = median(diff_past_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_plants_taxonomic_past_nest, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_island_age_sc)) +
  # geom_point(data = corals_taxonomic$beta_sne_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Island Age", y = "Beta diversity")

p4 <- plants_taxonomic %>% 
  modelr::data_grid(d_past_sc = median(d_past_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_past_sc = median(diff_emerse_past_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_past_air = seq_range(diff_past_air, n = 10)) %>% 
  add_epred_draws(model_beta_plants_taxonomic_past_nest, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_past_air)) +
  # geom_point(data = corals_taxonomic$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Air Temperature", y = "Beta diversity")

library(cowplot)

p1 / p2 / p3 / p4






#6. Terrestrial - PLANTS BETA FUNCTIONAL ----


model_plants_functional <- model_plants_functional %>%
  filter(!Island1 %in% "Farquhar_Atoll",
         !Island2 %in% "Farquhar_Atoll")



# Scaling variables - difference between archipelagos islands
plants_functional <- model_plants_functional %>%
  mutate(
    beta_sne_adj = beta_nest - 0.000001, 
    beta_sne_adj = ifelse(beta_sne_adj <= 0, 1e-6, beta_sne_adj),
    d_past_sc = scale(d_past.y), #emerse distance past
    diff_emerse_past_sc = scale(diff_emerse_past), #emerse island area past
    diff_emerse_pres_sc = scale(diff_emerse_pres), #emerse island area past
    diff_island_age_sc = scale(diff_island_age),
    distance_km_sc = scale(Distance_km), #emerse distance past
    diff_present_air = scale(diif_present_air),
    diff_past_air = scale(diff_past_air)
  )


#------------------------- Models

# Null

model_beta_plants_functional_past_null_nest <- brm(
  beta_sne_adj ~ 1 +
    (1 | Archipelago),
  data = plants_functional,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 5000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.98, max_treedepth = 14),
)


model_beta_plants_functional_past_null_nest <- add_criterion(model_beta_plants_functional_past_null_nest, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_plants_functional_past_null_nest,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_plants_functional_past_null_nest.rds"
)

# Assess model validity using LOO
loo(model_beta_plants_functional_past_null_nest)



# past

model_beta_plants_functional_past_nest <- brm(
  beta_sne_adj ~  diff_emerse_past_sc + diff_island_age_sc + d_past_sc + diff_past_air + 
    (1 | Archipelago),
  data = plants_functional,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_plants_functional_past_nest <- add_criterion(model_beta_plants_functional_past_nest, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_plants_functional_past_nest,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_plants_functional_past_nest.rds"
)

loo(model_beta_plants_functional_past_nest)

#loo_compare(model_beta_plants_functional_past_null_nest, model_beta_plants_functional_past)

#Look at model summary its on 'logit' scale, because of Beta family link function
model_beta_plants_functional_past_nest

## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_plants_functional_past_nest, ndraws = 50) 
pp_check(model_beta_plants_functional_past_nest, type = "ecdf_overlay", ndraws = 50) 

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_plants_functional_past$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_plants_functional_past_nest, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_plants_functional_past_nest, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_plants_functional_past_nest, type = 'rhat') 


#back transform model summary with point estimates (Ratios)
model_parameters(model_beta_plants_functional_past_nest, exponentiate = T, ci = 0.95)


#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####


#Plotting isolation effects
p1 <- plants_functional %>% 
  modelr::data_grid(d_past_sc = seq_range(d_past_sc, n = 100) %>% as.numeric,
                    diff_emerse_past_sc = median(diff_emerse_past_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_past_air = median(diff_past_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_plants_functional_past_nest, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = d_past_sc)) +
  # geom_point(data = corals_functional$beta_sne_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Distance between islands", y = "Beta diversity")

p2 <- plants_functional %>% 
  modelr::data_grid(d_past_sc = median(d_past_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_past_sc = seq_range(diff_emerse_past_sc, n = 10),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_past_air = median(diff_past_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_plants_functional_past_nest, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_emerse_past_sc)) +
  # geom_point(data = corals_functional$beta_sne_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Emersed Island Area", y = "Beta diversity")

p3 <- plants_functional %>% 
  modelr::data_grid(d_past_sc = median(d_past_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_past_sc = median(diff_emerse_past_sc, na.rm = T),
                    diff_island_age_sc = seq_range(diff_island_age_sc, n = 10),
                    diff_past_air = median(diff_past_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_plants_functional_past_nest, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_island_age_sc)) +
  # geom_point(data = corals_functional$beta_sne_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Island Age", y = "Beta diversity")

p4 <- plants_functional %>% 
  modelr::data_grid(d_past_sc = median(d_past_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_past_sc = median(diff_emerse_past_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_past_air = seq_range(diff_past_air, n = 10)) %>% 
  add_epred_draws(model_beta_plants_functional_past_nest, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_past_air)) +
  # geom_point(data = corals_taxonomic$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Air Temperature", y = "Beta diversity")

library(cowplot)

p1 / p2 / p3 / p4




#7. Terrestrial - BIRDS BETA TAXONOMIC ----


model_birds_taxonomic <- model_birds_taxonomic %>%
  filter(!Island1 %in% "Farquhar_Atoll",
         !Island2 %in% "Farquhar_Atoll")



# Scaling variables - difference between archipelagos islands
birds_taxonomic <- model_birds_taxonomic %>%
  mutate(
    beta_sne_adj = beta_sne - 0.000001,
    beta_sne_adj = ifelse(beta_sne_adj <= 0, 1e-6, beta_sne_adj),
    d_past_sc = scale(d_past.y), #emerse distance past
    diff_emerse_past_sc = scale(diff_emerse_past), #emerse island area past
    diff_emerse_pres_sc = scale(diff_emerse_pres), #emerse island area past
    diff_island_age_sc = scale(diff_island_age),
    distance_km_sc = scale(Distance_km), #emerse distance past
    diff_present_air = scale(diif_present_air),
    diff_past_air = scale(diff_past_air)
  )




#------------------------- Models

# Null

model_beta_birds_taxonomic_past_null_nest <- brm(
  beta_sne_adj ~ 1 +
    (1 | Archipelago),
  data = birds_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 5000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.98, max_treedepth = 14),
)


model_beta_birds_taxonomic_past_null_nest <- add_criterion(model_beta_birds_taxonomic_past_null_nest, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_birds_taxonomic_past_null_nest,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_birds_taxonomic_past_null_nest.rds"
)


# Assess model validity using LOO
loo(model_beta_birds_taxonomic_past_null_nest)


# past

model_beta_birds_taxonomic_past_nest <- brm(
  beta_sne_adj ~  diff_emerse_past_sc + diff_island_age_sc + d_past_sc + diff_past_air + 
    (1 | Archipelago),
  data = birds_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_birds_taxonomic_past_nest <- add_criterion(model_beta_birds_taxonomic_past_nest, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_birds_taxonomic_past_nest,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_birds_taxonomic_past_nest.rds"
)

loo(model_beta_birds_taxonomic_past_nest)

#loo_compare(model_beta_birds_taxonomic_past_null_nest, model_beta_birds_taxonomic_past_nest)

#Look at model summary its on 'logit' scale, because of Beta family link function
model_beta_birds_taxonomic_past_nest

## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_birds_taxonomic_past_nest, ndraws = 50) 

pp_check(model_beta_birds_taxonomic_past_nest, type = "ecdf_overlay", ndraws = 50) 

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_birds_taxonomic_past_nest$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_birds_taxonomic_past_nest, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_birds_taxonomic_past_nest, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_birds_taxonomic_past_nest, type = 'rhat') 


#back transform model summary with point estimates (Ratios)
model_parameters(model_beta_birds_taxonomic_past_nest, exponentiate = T, ci = 0.95)


#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####


#Plotting isolation effects
p1 <- birds_taxonomic %>% 
  modelr::data_grid(d_past_sc = seq_range(d_past_sc, n = 100) %>% as.numeric,
                    diff_emerse_past_sc = median(diff_emerse_past_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_past_air = median(diff_past_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_birds_taxonomic_past_nest, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = d_past_sc)) +
  # geom_point(data = corals_taxonomic$beta_sne_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Distance between islands", y = "Beta diversity")

p2 <- birds_taxonomic %>% 
  modelr::data_grid(d_past_sc = median(d_past_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_past_sc = seq_range(diff_emerse_past_sc, n = 10),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_past_air = median(diff_past_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_birds_taxonomic_past_nest, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_emerse_past_sc)) +
  # geom_point(data = corals_taxonomic$beta_sne_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Emersed Island Area", y = "Beta diversity")

p3 <- birds_taxonomic %>% 
  modelr::data_grid(d_past_sc = median(d_past_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_past_sc = median(diff_emerse_past_sc, na.rm = T),
                    diff_island_age_sc = seq_range(diff_island_age_sc, n = 10),
                    diff_past_air = median(diff_past_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_birds_taxonomic_past_nest, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_island_age_sc)) +
  # geom_point(data = corals_taxonomic$beta_sne_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Island Age", y = "Beta diversity")

p4 <- birds_taxonomic %>% 
  modelr::data_grid(d_past_sc = median(d_past_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_past_sc = median(diff_emerse_past_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_past_air = seq_range(diff_past_air, n = 10)) %>% 
  add_epred_draws(model_beta_birds_taxonomic_past_nest, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_past_air)) +
  # geom_point(data = corals_taxonomic$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Air Temperature", y = "Beta diversity")

library(cowplot)

p1 / p2 / p3 / p4




#8. Terrestrial - BIRDS BETA FUNCTIONAL ----


model_birds_functional <- model_birds_functional %>%
  filter(!Island1 %in% "Farquhar_Atoll",
         !Island2 %in% "Farquhar_Atoll")



# Scaling variables - difference between archipelagos islands
birds_functional <- model_birds_functional %>%
  mutate(
    beta_sne_adj = beta_nest - 0.000001,
    beta_sne_adj = ifelse(beta_sne_adj <= 0, 1e-6, beta_sne_adj),
    d_past_sc = scale(d_past.y), #emerse distance past
    diff_emerse_past_sc = scale(diff_emerse_past), #emerse island area past
    diff_emerse_pres_sc = scale(diff_emerse_pres), #emerse island area past
    diff_island_age_sc = scale(diff_island_age),
    distance_km_sc = scale(Distance_km), #emerse distance past
    diff_present_air = scale(diif_present_air),
    diff_past_air = scale(diff_past_air)
  )




#------------------------- Models

# Null

model_beta_birds_functional_past_null_nest <- brm(
  beta_sne_adj ~ 1 +
    (1 | Archipelago),
  data = birds_functional,
  family = beta_family(),
  chains = 4,
  cores = 4,
  iter = 5000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.98, max_treedepth = 14),
)


model_beta_birds_functional_past_null_nest <- add_criterion(model_beta_birds_functional_past_null_nest, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_birds_functional_past_null_nest,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_birds_functional_past_null_nest.rds"
)

# Assess model validity using LOO
loo(model_beta_birds_functional_past_null_nest)


# past

model_beta_birds_functional_past_nest <- brm(
  beta_sne_adj ~  diff_emerse_past_sc + diff_island_age_sc + d_past_sc + diff_past_air +
    (1 | Archipelago),
  data = birds_functional,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_birds_functional_past_nest <- add_criterion(model_beta_birds_functional_past_nest, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_birds_functional_past_nest,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_birds_functional_past_nest.rds"
)

loo(model_beta_birds_functional_past_nest)

#loo_compare(model_beta_birds_functional_past_null_nest, model_beta_birds_functional_past)

#Look at model summary its on 'logit' scale, because of Beta family link function
model_beta_birds_functional_past_nest

## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_birds_functional_past_nest, ndraws = 50) 
pp_check(model_beta_birds_functional_past_nest, type = "ecdf_overlay", ndraws = 50) 

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_birds_functional_past_nest$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_birds_functional_past_nest, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_birds_functional_past_nest, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_birds_functional_past_nest, type = 'rhat') 


#back transform model summary with point estimates (Ratios)
model_parameters(model_beta_birds_functional_past, exponentiate = T, ci = 0.95)


#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####


#Plotting isolation effects
p1 <- birds_functional %>% 
  modelr::data_grid(d_past_sc = seq_range(d_past_sc, n = 100) %>% as.numeric,
                    diff_emerse_past_sc = median(diff_emerse_past_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_past_air = median(diff_past_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_birds_functional_past_nest, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = d_past_sc)) +
  # geom_point(data = corals_functional$beta_sne_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Distance between islands", y = "Beta diversity")

p2 <- birds_functional %>% 
  modelr::data_grid(d_past_sc = median(d_past_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_past_sc = seq_range(diff_emerse_past_sc, n = 10),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_past_air = median(diff_past_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_birds_functional_past_nest, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_emerse_past_sc)) +
  # geom_point(data = corals_functional$beta_sne_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Emersed Island Area", y = "Beta diversity")

p3 <- birds_functional %>% 
  modelr::data_grid(d_past_sc = median(d_past_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_past_sc = median(diff_emerse_past_sc, na.rm = T),
                    diff_island_age_sc = seq_range(diff_island_age_sc, n = 10),
                    diff_past_air = median(diff_past_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_birds_functional_past_nest, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_island_age_sc)) +
  # geom_point(data = corals_functional$beta_sne_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Island Age", y = "Beta diversity")

p4 <- birds_functional %>% 
  modelr::data_grid(d_past_sc = median(d_past_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_past_sc = median(diff_emerse_past_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_past_air = seq_range(diff_past_air, n = 10)) %>% 
  add_epred_draws(model_beta_birds_functional_past_nest, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_past_air)) +
  # geom_point(data = corals_taxonomic$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Air Temperature", y = "Beta diversity")

library(cowplot)

p1 / p2 / p3 / p4





######################################################################################################################################################################


# TURNOVER ----

#-------------------------------------------------------------------------------


#1. Marine - CORALS BETA TAXONOMIC ----


## Scaling variables - difference between marine regions ----
corals_taxonomic <- model_corals_taxonomic %>%
  mutate(
    beta_sim_adj = beta_sim - 0.000001, 
    beta_sim_adj = ifelse(beta_sim_adj <= 0, 1e-6, beta_sim_adj),
    diff_isolation_sc = scale(diff_isolation)[,1],
    diff_reef_area_sc = scale(diff_reef_area)[,1],
    diff_past_isolation_sc = scale(diff_past_isolation)[,1], #quaternary distance
    diff_past_reef_area_sc = scale(diff_past_reef_area)[,1],
    diff_age_sc = scale(diff_age)[,1],
    dist_sc = scale(Distance)[,1],
    diff_present_sst = scale(diff_sst)[,1],
    diff_past_sst = scale(diff_past_sst)[,1]
  )

corals_taxonomic %>% glimpse

#------------------------- Models

# Present


# Bayes models
# Null

# model_beta_corals_taxonomic_present_null_turn_turngp3 <- brm(
#   beta_sim_adj ~ 1 + (1|Group1) + (1|Group2) + (1|key),
#   data = corals_taxonomic,
#   family = beta_family(),
#   chains = 4,
#   cores = 4,
#   iter = 4000,
#   warmup = 3000,
#   control = list(adapt_delta = 0.99, max_treedepth = 14),
#   save_pars = save_pars(all = T)
# )

# model_beta_corals_taxonomic_present_null_turn_turn_gp3 <- add_criterion(model_beta_corals_taxonomic_present_null_turn_turn_gp3, "loo", save_psis = TRUE, reloo = T)


# model_beta_corals_taxonomic_present_null <- brm(
#   beta_sim_adj ~ 1,
#   data = corals_taxonomic,
#   family = beta_family(),
#   chains = 4,
#   cores = 4,
#   iter = 4000,
#   warmup = 3000,
#   control = list(adapt_delta = 0.99, max_treedepth = 14),
#   save_pars = save_pars(all = T)
# )

# model_beta_corals_taxonomic_present_null <- add_criterion(model_beta_corals_taxonomic_present_null, "loo", save_psis = TRUE, reloo = T)
# 
# loo_compare(model_beta_corals_taxonomic_present_null, model_beta_corals_taxonomic_present_null_turn_turn_gp2, model_beta_corals_taxonomic_present_null_turn_turn_gp1, model_beta_corals_taxonomic_present_null_turn_turn_gp3)

model_beta_corals_taxonomic_present_null_turn <- brm(
  beta_sim_adj ~ 1 + 
    (1+dist_sc|Group1) + (1+dist_sc|Group2),
  data = corals_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  control = list(adapt_delta = 0.99, max_treedepth = 14),
  save_pars = save_pars(all = T)
)

model_beta_corals_taxonomic_present_null_turn <- add_criterion(model_beta_corals_taxonomic_present_null_turn, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_corals_taxonomic_present_null_turn,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_corals_taxonomic_present_null_turn.rds"
)

# Assess model validity using LOO
loo(model_beta_corals_taxonomic_present_null_turn)

# Present
model_beta_corals_taxonomic_present_bayes_turn <- brm(
  beta_sim_adj ~ diff_isolation_sc + diff_reef_area_sc + diff_age_sc + dist_sc + diff_present_sst +
    (1+dist_sc|Group1) + (1+dist_sc|Group2),
  data = corals_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_corals_taxonomic_present_bayes_turn <- add_criterion(model_beta_corals_taxonomic_present_bayes_turn, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_corals_taxonomic_present_bayes_turn,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_corals_taxonomic_present_bayes_turn.rds"
)

model_beta_corals_taxonomic_present_bayes_turn

# Assess model validity using LOO
loo(model_beta_corals_taxonomic_present_bayes_turn)

#Look at model summary its on 'logit' scale, because of bernouli family link function
model_beta_corals_taxonomic_present_bayes_turn

#loo_compare(model_beta_corals_taxonomic_present_bayes_turn, model_beta_corals_taxonomic_present_null_turn)
## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_corals_taxonomic_present_bayes_turn , ndraws = 50) 
pp_check(model_beta_corals_taxonomic_present_bayes_turn, type = "ecdf_overlay", ndraws = 50) 

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_corals_taxonomic_present_bayes_turn$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_corals_taxonomic_present_bayes_turn, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_corals_taxonomic_present_bayes_turn, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_corals_taxonomic_present_bayes_turn, type = 'rhat') 


#back transform model summary with point estimates - Now Odds ratios (3:1 odds of bleaching)
model_parameters(model_beta_corals_taxonomic_present_bayes_turn, exponentiate = T, ci = 0.95)



#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####

#Plotting isolation effects
# 
p1 <- corals_taxonomic %>% 
  modelr::data_grid(
    diff_isolation_sc = seq_range(diff_isolation_sc, n = 10),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_present_sst = median(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_taxonomic_present_bayes_turn,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_isolation_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in isolation", y = "Beta diversity")

#
p2 <- corals_taxonomic %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = seq_range(diff_reef_area_sc, n = 10),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_present_sst = median(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_taxonomic_present_bayes_turn,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_reef_area_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in reef area", y = "Beta diversity")

# 
p3 <- corals_taxonomic %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    diff_age_sc = seq_range(diff_age_sc, n = 10),
    dist_sc = mean(dist_sc),
    diff_present_sst = median(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_taxonomic_present_bayes_turn,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_age_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in age", y = "Beta diversity")

# 
p4 <- corals_taxonomic %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    dist_sc = seq_range(dist_sc, n = 10),
    diff_age_sc = mean(diff_age_sc),
    diff_present_sst = median(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_taxonomic_present_bayes_turn,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = dist_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch distance", y = "Beta diversity")

#
p5 <- corals_taxonomic %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    dist_sc = mean(dist_sc),
    diff_age_sc = mean(diff_age_sc),
    diff_present_sst = seq_range(diff_present_sst, n = 10)
  ) %>% 
  add_epred_draws(model_beta_corals_taxonomic_present_bayes_turn,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_present_sst)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch SST", y = "Beta diversity")


library(cowplot)

p1 / p2 / p3 / p4 / p5






#2. Marine - CORALS BETA FUNCTIONAL ----


## Scaling variables - difference between marine regions ----
corals_functional <- model_corals_functional %>%
  mutate(
    beta_sim_adj = beta_sim - 0.000001, 
    beta_sim_adj = ifelse(beta_sim_adj <= 0, 1e-6, beta_sim_adj),
    diff_isolation_sc = scale(diff_isolation)[,1],
    diff_reef_area_sc = scale(diff_reef_area)[,1],
    diff_past_isolation_sc = scale(diff_past_isolation)[,1], #quaternary distance
    diff_past_reef_area_sc = scale(diff_past_reef_area)[,1],
    diff_age_sc = scale(diff_age)[,1],
    dist_sc = scale(Distance)[,1],
    diff_present_sst = scale(diff_sst)[,1],
    diff_past_sst = scale(diff_past_sst)[,1]
  )

corals_functional %>% glimpse

#------------------------- Models


# Bayes models
# Null

model_beta_corals_functional_present_null_turn_gp2 <- brm(
  beta_sim_adj ~ 1 + 
    (1+dist_sc|Group1) + (1+dist_sc|Group2),
  data = corals_functional,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  control = list(adapt_delta = 0.99, max_treedepth = 14),
  save_pars = save_pars(all = T)
)

model_beta_corals_functional_present_null_turn_gp2 <- add_criterion(model_beta_corals_functional_present_null_turn_gp2, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_corals_functional_present_null_turn_gp2,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_corals_functional_present_null_turn_gp2.rds"
)

# Assess model validity using LOO
loo(model_beta_corals_functional_present_null_turn_gp2)

# Present
model_beta_corals_functional_present_bayes_turn <- brm(
  beta_sim_adj ~ diff_isolation_sc + diff_reef_area_sc + diff_age_sc + dist_sc + diff_present_sst + 
  (1+dist_sc|Group1) + (1+dist_sc|Group2),
  data = corals_functional,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_corals_functional_present_bayes_turn <- add_criterion(model_beta_corals_functional_present_bayes_turn, "loo", save_psis = TRUE, reloo = T)

# Assess model validity using LOO
loo(model_beta_corals_functional_present_bayes_turn)

#Look at model summary its on 'logit' scale, because of bernouli family link function
model_beta_corals_functional_present_bayes_turn

#loo_compare(model_beta_corals_functional_present_bayes_turn, model_beta_corals_functional_present_null_turn_gp2)

## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_corals_functional_present_null_turn_gp2 , ndraws = 50) 
pp_check(model_beta_corals_functional_present_bayes_turn, type = "ecdf_overlay", ndraws = 50) 

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_corals_functional_present_bayes_turn$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_corals_functional_present_bayes_turn, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_corals_functional_present_bayes_turn, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_corals_functional_present_bayes_turn, type = 'rhat') 


#back transform model summary with point estimates - Now Odds ratios (3:1 odds of bleaching)
model_parameters(model_beta_corals_functional_present_bayes_turn, exponentiate = T, ci = 0.95)


#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####

#Plotting isolation effects
p1 <- corals_functional %>% 
  modelr::data_grid(
    diff_isolation_sc = seq_range(diff_isolation_sc, n = 10),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_present_sst = mean(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_functional_present_bayes_turn,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_isolation_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in isolation", y = "Beta diversity")

#
p2 <- corals_functional %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = seq_range(diff_reef_area_sc, n = 10),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_present_sst = mean(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_functional_present_bayes_turn,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_reef_area_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in reef area", y = "Beta diversity")

# 
p3 <- corals_functional %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    diff_age_sc = seq_range(diff_age_sc, n = 10),
    dist_sc = mean(dist_sc),
    diff_present_sst = mean(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_functional_present_bayes_turn,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_age_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in age", y = "Beta diversity")

# 
p4 <- corals_functional %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    dist_sc = seq_range(dist_sc, n = 10),
    diff_age_sc = mean(diff_age_sc),
    diff_present_sst = mean(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_functional_present_bayes_turn,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = dist_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch distance", y = "Beta diversity")

#
p5 <- corals_functional %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    dist_sc = mean(dist_sc),
    diff_age_sc = mean(diff_age_sc),
    diff_present_sst = seq_range(diff_present_sst, n = 10)
  ) %>% 
  add_epred_draws(model_beta_corals_functional_present_bayes_turn,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_present_sst)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch SST", y = "Beta diversity")


library(cowplot)

p1 / p2 / p3 / p4 / p5




#3. Marine - FISH BETA TAXONOMIC ----


## Scaling variables - difference between marine regions ----
fish_taxonomic <- model_fish_taxonomic %>%
  mutate(
    beta_sim_adj = beta_sim - 0.000001,
    beta_sim_adj = ifelse(beta_sim_adj <= 0, 1e-6, beta_sim_adj),
    diff_isolation_sc = scale(diff_isolation)[,1],
    diff_reef_area_sc = scale(diff_reef_area)[,1],
    diff_past_isolation_sc = scale(diff_past_isolation)[,1], #quaternary distance
    diff_past_reef_area_sc = scale(diff_past_reef_area)[,1],
    diff_age_sc = scale(diff_age)[,1],
    dist_sc = scale(Distance)[,1],
    diff_present_sst = scale(diff_sst)[,1],
    diff_past_sst = scale(diff_past_sst)[,1]
     )

fish_taxonomic %>% glimpse



#------------------------- Models


# Bayes models
# Null


# model_beta_fish_taxonomic_present_null_turn_gp3 <- add_criterion(model_beta_fish_taxonomic_present_null_turn_gp3, "loo", save_psis = TRUE, reloo = T)
# model_beta_fish_taxonomic_present_null_turn_gp2 <- brm(
#   beta_sim_adj ~ 1 + (1|Group1) + (1|Group2) ,
#   data = fish_taxonomic,
#   family = beta_family(),
#   chains = 4,
#   cores = 4,
#   iter = 4000,
#   warmup = 3000,
#   control = list(adapt_delta = 0.99, max_treedepth = 14),
#   save_pars = save_pars(all = T)
# )

model_beta_fish_taxonomic_present_null_turn_gp2 <- brm(
  beta_sim_adj ~ 1 + (1+dist_sc|Group1) + (1+dist_sc|Group2),
  data = fish_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  control = list(adapt_delta = 0.99, max_treedepth = 14),
  save_pars = save_pars(all = T)
)

model_beta_fish_taxonomic_present_null_turn_gp2 <- add_criterion(model_beta_fish_taxonomic_present_null_turn_gp2, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_fish_taxonomic_present_null_turn_gp2,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_fish_taxonomic_present_null_turn_gp2.rds"
)

# Assess model validity using LOO
loo(model_beta_fish_taxonomic_present_null_turn_gp2)

#Look at model summary its on 'logit' scale, because of bernouli family link function
model_beta_fish_taxonomic_present_null_turn_gp2


# Present
model_beta_fish_taxonomic_present_bayes_turn <- brm(
  beta_sim_adj ~ diff_isolation_sc + diff_reef_area_sc + diff_age_sc + dist_sc + diff_present_sst + 
  (1+dist_sc|Group1) + (1+dist_sc|Group2),
  data = fish_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_fish_taxonomic_present_bayes_turn <- add_criterion(model_beta_fish_taxonomic_present_bayes_turn, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_fish_taxonomic_present_bayes_turn,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_fish_taxonomic_present_bayes_turn.rds"
)

model_beta_fish_taxonomic_present_bayes_turn

# Assess model validity using LOO
loo(model_beta_fish_taxonomic_present_bayes_turn)

#loo_compare(model_beta_fish_taxonomic_present_bayes_turn, model_beta_fish_taxonomic_present_null_turn_gp2)
## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_fish_taxonomic_present_bayes_turn , ndraws = 50) 
pp_check(model_beta_fish_taxonomic_present_bayes_turn, type = "ecdf_overlay", ndraws = 50) 

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_fish_taxonomic_present_bayes_turn$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_fish_taxonomic_present_bayes_turn, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_fish_taxonomic_present_bayes_turn, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_fish_taxonomic_present_bayes_turn, type = 'rhat') 


#back transform model summary with point estimates - Now Odds ratios (3:1 odds of bleaching)
model_parameters(model_beta_fish_taxonomic_present_bayes_turn, exponentiate = T, ci = 0.95)


#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####

#Plotting isolation effects
p1 <- fish_taxonomic %>% 
  modelr::data_grid(
    diff_isolation_sc = seq_range(diff_isolation_sc, n = 10),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_present_sst = mean(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_taxonomic_present_bayes_turn,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_isolation_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in isolation", y = "Beta diversity")

#
p2 <- fish_taxonomic %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = seq_range(diff_reef_area_sc, n = 10),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_present_sst = mean(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_taxonomic_present_bayes_turn,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_reef_area_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in reef area", y = "Beta diversity")

# 
p3 <- fish_taxonomic %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    diff_age_sc = seq_range(diff_age_sc, n = 10),
    dist_sc = mean(dist_sc),
    diff_present_sst = mean(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_taxonomic_present_bayes_turn,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_age_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in age", y = "Beta diversity")

# 
p4 <- fish_taxonomic %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    dist_sc = seq_range(dist_sc, n = 10),
    diff_age_sc = mean(diff_age_sc),
    diff_present_sst = mean(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_taxonomic_present_bayes_turn,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = dist_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch distance", y = "Beta diversity")

#
p5 <- fish_taxonomic %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    dist_sc = mean(dist_sc),
    diff_age_sc = mean(diff_age_sc),
    diff_present_sst = seq_range(diff_present_sst, n = 10)
  ) %>% 
  add_epred_draws(model_beta_fish_taxonomic_present_bayes_turn,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_present_sst)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch SST", y = "Beta diversity")


library(cowplot)

p1 / p2 / p3 / p4 / p5




#4. Marine - FISH BETA FUNCTIONAL ----


## Scaling variables - difference between marine regions ----
fish_functional <- model_fish_functional %>%
  mutate(
    beta_sim_adj = beta_sim - 0.000001,  
    beta_sim_adj = ifelse(beta_sim_adj <= 0, 1e-6, beta_sim_adj),
    diff_isolation_sc = scale(diff_isolation)[,1],
    diff_reef_area_sc = scale(diff_reef_area)[,1],
    diff_past_isolation_sc = scale(diff_past_isolation)[,1], #quaternary distance
    diff_past_reef_area_sc = scale(diff_past_reef_area)[,1],
    diff_age_sc = scale(diff_age)[,1],
    dist_sc = scale(Distance)[,1],
    diff_present_sst = scale(diff_sst)[,1],
    diff_past_sst = scale(diff_past_sst)[,1]
  ) 

fish_functional %>% glimpse

#------------------------- Models


# Bayes models
# Null
model_beta_fish_functional_present_null_turn_gp2 <- brm(
  beta_sim_adj ~ 1 + (1+dist_sc|Group1) + (1+dist_sc|Group2) ,
  data = fish_functional,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  control = list(adapt_delta = 0.99, max_treedepth = 14),
  save_pars = save_pars(all = T)
)

model_beta_fish_functional_present_null_turn_gp2 <- add_criterion(model_beta_fish_functional_present_null_turn_gp2, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_fish_functional_present_null_turn_gp2,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_fish_functional_present_null_turn_gp2.rds"
)

# Assess model validity using LOO
loo(model_beta_fish_functional_present_null_turn_gp2)


# Present
model_beta_fish_functional_present_bayes_turn <- brm(
  beta_sim_adj ~  1 + diff_isolation_sc + diff_reef_area_sc + diff_age_sc + dist_sc + diff_present_sst + 
  (1+dist_sc|Group1) + (1+dist_sc|Group2),
  data = fish_functional,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_fish_functional_present_bayes_turn <- add_criterion(model_beta_fish_functional_present_bayes_turn, "loo", save_psis = TRUE, reloo = T, moment_match = TRUE)

saveRDS(
  model_beta_fish_functional_present_bayes_turn,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_fish_functional_present_bayes_turn.rds"
)

# Assess model validity using LOO
loo(model_beta_fish_functional_present_bayes_turn)

#loo_compare(model_beta_fish_functional_present_bayes_turn, model_beta_fish_functional_present_null_turn_gp2)

## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_fish_functional_present_bayes_turn , ndraws = 50) 
pp_check(model_beta_fish_functional_present_bayes_turn, type = "ecdf_overlay", ndraws = 50) 
pp_check(model_beta_fish_functional_present_bayes_turn , ndraws = 50, type = "dens_overlay_grouped", group = "Group1")
# pp_check(model_beta_fish_functional_present_bayes_turn_g , ndraws = 50, type = "scatter_avg_grouped", group = "Group1")
# pp_check(model_beta_fish_functional_present_bayes_turn_g , ndraws = 50, type = "error_scatter_avg_grouped", group = "Group1")
# pp_check(model_beta_fish_functional_present_bayes_turn_g, type='pit_ecdf_grouped', group = "Group1", prob = 0.95, plot_diff = F)
# pp_check(model_beta_fish_functional_present_bayes_turn_g, type='pit_ecdf_grouped', group = "Group2", prob = 0.95, plot_diff = F)
# 
# pp_check(model_beta_fish_functional_present_bayes_turn_g , ndraws = 50, type = "scatter_avg_grouped", group = "Group2")

#pairs(model_beta_fish_functional_present_bayes_turn_g, variable = variables(model_beta_fish_functional_present_bayes_turn_g)[1:5])
#fish_functional %>% select(Group1, Group2) %>% unique %>% view

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_fish_functional_present_bayes_turn$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_fish_functional_present_bayes_turn, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_fish_functional_present_bayes_turn, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_fish_functional_present_bayes_turn, type = 'rhat') 


#back transform model summary with point estimates - Now Odds ratios (3:1 odds of bleaching)
model_parameters(model_beta_fish_functional_present_bayes_turn, exponentiate = T, ci = 0.95)



#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####

#Plotting isolation effects
p1 <- fish_functional %>% 
  modelr::data_grid(
    diff_isolation_sc = seq_range(diff_isolation_sc, n = 10),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_present_sst = mean(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_functional_present_bayes_turn,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_isolation_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in isolation", y = "Beta diversity")

#
p2 <- fish_functional %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = seq_range(diff_reef_area_sc, n = 10),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_present_sst = mean(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_functional_present_bayes_turn,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_reef_area_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in reef area", y = "Beta diversity")

# 
p3 <- fish_functional %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    diff_age_sc = seq_range(diff_age_sc, n = 10),
    dist_sc = mean(dist_sc),
    diff_present_sst = mean(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_functional_present_bayes_turn,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_age_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in age", y = "Beta diversity")

# 
p4 <- fish_functional %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    dist_sc = seq_range(dist_sc, n = 10),
    diff_age_sc = mean(diff_age_sc),
    diff_present_sst = mean(diff_present_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_functional_present_bayes_turn,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = dist_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch distance", y = "Beta diversity")

#
p5 <- fish_functional %>% 
  modelr::data_grid(
    diff_isolation_sc = mean(diff_isolation_sc),
    diff_reef_area_sc = mean(diff_reef_area_sc),
    dist_sc = mean(dist_sc),
    diff_age_sc = mean(diff_age_sc),
    diff_present_sst = seq_range(diff_present_sst, n = 10)
  ) %>% 
  add_epred_draws(model_beta_fish_functional_present_bayes_turn,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_present_sst)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch SST", y = "Beta diversity")


library(cowplot)

p1 / p2 / p3 / p4 / p5





#5. Terrestrial - PLANTS BETA TAXONOMIC ----


model_plants_taxonomic <- model_plants_taxonomic %>%
  filter(!Island1 %in% "Farquhar_Atoll",
         !Island2 %in% "Farquhar_Atoll")


# Scaling variables - difference between archipelagos islands
plants_taxonomic <- model_plants_taxonomic %>%
  mutate(
    beta_sim_adj = beta_sim - 0.000001,
    beta_sim_adj = ifelse(beta_sim_adj <= 0, 1e-6, beta_sim_adj),
    d_past_sc = scale(d_past.y), #emerse distance past
    diff_emerse_past_sc = scale(diff_emerse_past), #emerse island area present
    diff_emerse_pres_sc = scale(diff_emerse_pres), #emerse island area past
    diff_island_age_sc = scale(diff_island_age),
    distance_km_sc = scale(Distance_km), #emerse distance present
    diff_present_air = scale(diif_present_air),
    diff_past_air = scale(diff_past_air)
  )

# scaled_diff_emerse_pres <- plants_taxonomic$diff_emerse_pres - mean(plants_taxonomic$diff_emerse_pres)/sd(plants_taxonomic$diff_emerse_pres)
# 
# diff_emerse_pres <- scaled_diff_emerse_pres * sd(plants_taxonomic$diff_emerse_pres) + mean(plants_taxonomic$diff_emerse_pres)
# 
# diff_emerse_pres_sc <-  scale(plants_taxonomic$diff_island_age_sc)
# diff_emerse_pres_sc %>% glimpse


#------------------------- Models

# Null

model_beta_plants_taxonomic_present_turn_null_turn <- brm(
  beta_sim_adj ~ 1 +
    (1 | Archipelago),
  data = plants_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 5000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.98, max_treedepth = 14),
)


model_beta_plants_taxonomic_present_turn_null_turn <- add_criterion(model_beta_plants_taxonomic_present_turn_null_turn, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_plants_taxonomic_present_turn_null_turn,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_plants_taxonomic_present_turn_null_turn.rds"
)

# Assess model validity using LOO
loo(model_beta_plants_taxonomic_present_turn_null_turn)



# Present

model_beta_plants_taxonomic_present_turn <- brm(
  beta_sim_adj ~  diff_emerse_pres_sc + diff_island_age_sc + distance_km_sc + diff_present_air +
    (1 | Archipelago),
  data = plants_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_plants_taxonomic_present_turn <- add_criterion(model_beta_plants_taxonomic_present_turn, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_plants_taxonomic_present_turn,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_plants_taxonomic_present_turn.rds"
)

loo(model_beta_plants_taxonomic_present_turn)

#loo_compare(model_beta_plants_taxonomic_present_turn_null_turn, model_beta_plants_taxonomic_present_turn)

#Look at model summary its on 'logit' scale, because of Beta family link function
model_beta_plants_taxonomic_present_turn

## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_plants_taxonomic_present_turn, ndraws = 50) 
pp_check(model_beta_plants_taxonomic_present_turn, type = "ecdf_overlay", ndraws = 50) 

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_plants_taxonomic_present_turn$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_plants_taxonomic_present_turn, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_plants_taxonomic_present_turn, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_plants_taxonomic_present_turn, type = 'rhat') 


#back transform model summary with point estimates (Ratios)
model_parameters(model_beta_plants_taxonomic_present_turn, exponentiate = T, ci = 0.95)


#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####


#Plotting isolation effects
p1 <- plants_taxonomic %>% 
  modelr::data_grid(distance_km_sc = seq_range(distance_km_sc, n = 100) %>% as.numeric,
                    diff_emerse_pres_sc = median(diff_emerse_pres_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_present_air = median(diff_present_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_plants_taxonomic_present_turn, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = distance_km_sc)) +
  # geom_point(data = corals_taxonomic$beta_sim_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Distance between islands", y = "Beta diversity")

p2 <- plants_taxonomic %>% 
  modelr::data_grid(distance_km_sc = median(distance_km_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_pres_sc = seq_range(diff_emerse_pres_sc, n = 10),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_present_air = median(diff_present_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_plants_taxonomic_present_turn, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_emerse_pres_sc)) +
  # geom_point(data = corals_taxonomic$beta_sim_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Emersed Island Area", y = "Beta diversity")

p3 <- plants_taxonomic %>% 
  modelr::data_grid(distance_km_sc = median(distance_km_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_pres_sc = median(diff_emerse_pres_sc, na.rm = T),
                    diff_island_age_sc = seq_range(diff_island_age_sc, n = 10),
                    diff_present_air = median(diff_present_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_plants_taxonomic_present_turn, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_island_age_sc)) +
  # geom_point(data = corals_taxonomic$beta_sim_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Island Age", y = "Beta diversity")

p4 <- plants_taxonomic %>% 
  modelr::data_grid(distance_km_sc = median(distance_km_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_pres_sc = median(diff_emerse_pres_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_present_air = seq_range(diff_present_air, n = 10)) %>% 
  add_epred_draws(model_beta_plants_taxonomic_present_turn, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_present_air)) +
  # geom_point(data = corals_taxonomic$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Air Temperature", y = "Beta diversity")

library(cowplot)

p1 / p2 / p3 / p4






#6. Terrestrial - PLANTS BETA FUNCTIONAL ----


model_plants_functional <- model_plants_functional %>%
  filter(!Island1 %in% "Farquhar_Atoll",
         !Island2 %in% "Farquhar_Atoll")



# Scaling variables - difference between archipelagos islands
plants_functional <- model_plants_functional %>%
  mutate(
    beta_sim_adj = beta_turn - 0.000001, 
    beta_sim_adj = ifelse(beta_sim_adj <= 0, 1e-6, beta_sim_adj),
    d_past_sc = scale(d_past.y), #emerse distance past
    diff_emerse_past_sc = scale(diff_emerse_past), #emerse island area present
    diff_emerse_pres_sc = scale(diff_emerse_pres), #emerse island area past
    diff_island_age_sc = scale(diff_island_age),
    distance_km_sc = scale(Distance_km), #emerse distance present
    diff_present_air = scale(diif_present_air),
    diff_past_air = scale(diff_past_air)
  )


#------------------------- Models

# Null

model_beta_plants_functional_present_null_turn <- brm(
  beta_sim_adj ~ 1 +
    (1 | Archipelago),
  data = plants_functional,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 5000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.98, max_treedepth = 14),
)


model_beta_plants_functional_present_null_turn <- add_criterion(model_beta_plants_functional_present_null_turn, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_plants_functional_present_null_turn,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_plants_functional_present_null_turn.rds"
)

# Assess model validity using LOO
loo(model_beta_plants_functional_present_null_turn)



# Present

model_beta_plants_functional_present_turn <- brm(
  beta_sim_adj ~  diff_emerse_pres_sc + diff_island_age_sc + distance_km_sc + diff_present_air +
    (1 | Archipelago),
  data = plants_functional,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_plants_functional_present_turn <- add_criterion(model_beta_plants_functional_present_turn, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_plants_functional_present_turn,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_plants_functional_present_turn.rds"
)

loo(model_beta_plants_functional_present_turn)

#loo_compare(model_beta_plants_functional_present_turn, model_beta_plants_functional_present)

#Look at model summary its on 'logit' scale, because of Beta family link function
model_beta_plants_functional_present_turn

## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_plants_functional_present_turn, ndraws = 50) 
pp_check(model_beta_plants_functional_present_turn, type = "ecdf_overlay", ndraws = 50) 

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_plants_functional_present_turn$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_plants_functional_present_turn, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_plants_functional_present_turn, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_plants_functional_present_turn, type = 'rhat') 


#back transform model summary with point estimates (Ratios)
model_parameters(model_beta_plants_functional_present_turn, exponentiate = T, ci = 0.95)


#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####


#Plotting isolation effects
p1 <- plants_functional %>% 
  modelr::data_grid(distance_km_sc = seq_range(distance_km_sc, n = 100) %>% as.numeric,
                    diff_emerse_pres_sc = median(diff_emerse_pres_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_present_air = median(diff_present_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_plants_functional_present_turn, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = distance_km_sc)) +
  # geom_point(data = corals_functional$beta_sim_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Distance between islands", y = "Beta diversity")

p2 <- plants_functional %>% 
  modelr::data_grid(distance_km_sc = median(distance_km_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_pres_sc = seq_range(diff_emerse_pres_sc, n = 10),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_present_air = median(diff_present_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_plants_functional_present_turn, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_emerse_pres_sc)) +
  # geom_point(data = corals_functional$beta_sim_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Emersed Island Area", y = "Beta diversity")

p3 <- plants_functional %>% 
  modelr::data_grid(distance_km_sc = median(distance_km_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_pres_sc = median(diff_emerse_pres_sc, na.rm = T),
                    diff_island_age_sc = seq_range(diff_island_age_sc, n = 10),
                    diff_present_air = median(diff_present_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_plants_functional_present_turn, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_island_age_sc)) +
  # geom_point(data = corals_functional$beta_sim_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Island Age", y = "Beta diversity")

p4 <- plants_functional %>% 
  modelr::data_grid(distance_km_sc = median(distance_km_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_pres_sc = median(diff_emerse_pres_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_present_air = seq_range(diff_present_air, n = 10)) %>% 
  add_epred_draws(model_beta_plants_functional_present_turn, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_present_air)) +
  # geom_point(data = corals_taxonomic$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Air Temperature", y = "Beta diversity")

p1 / p2 / p3 / p4





#7. Terrestrial - BIRDS BETA TAXONOMIC ----


model_birds_taxonomic <- model_birds_taxonomic %>%
  filter(!Island1 %in% "Farquhar_Atoll",
         !Island2 %in% "Farquhar_Atoll")



# Scaling variables - difference between archipelagos islands
birds_taxonomic <- model_birds_taxonomic %>%
  mutate(
    beta_sim_adj = beta_sim - 0.000001,
    beta_sim_adj = ifelse(beta_sim_adj <= 0, 1e-6, beta_sim_adj),
    d_past_sc = scale(d_past.y), #emerse distance past
    diff_emerse_past_sc = scale(diff_emerse_past), #emerse island area present
    diff_emerse_pres_sc = scale(diff_emerse_pres), #emerse island area past
    diff_island_age_sc = scale(diff_island_age),
    distance_km_sc = scale(Distance_km), #emerse distance present
    diff_present_air = scale(diif_present_air),
    diff_past_air = scale(diff_past_air)
  )




#------------------------- Models

# Null

model_beta_birds_taxonomic_present_null_turn <- brm(
  beta_sim_adj ~ 1 +
    (1 | Archipelago),
  data = birds_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 5000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.98, max_treedepth = 14),
)


model_beta_birds_taxonomic_present_null_turn <- add_criterion(model_beta_birds_taxonomic_present_null_turn, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_birds_taxonomic_present_null_turn,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_birds_taxonomic_present_null_turn.rds"
)

# Assess model validity using LOO
loo(model_beta_birds_taxonomic_present_null_turn)



# Present

model_beta_birds_taxonomic_present_turn <- brm(
  beta_sim_adj ~  diff_emerse_pres_sc + diff_island_age_sc + distance_km_sc + diff_present_air +
    (1 | Archipelago),
  data = birds_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_birds_taxonomic_present_turn <- add_criterion(model_beta_birds_taxonomic_present_turn, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_birds_taxonomic_present_turn,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_birds_taxonomic_present_turn.rds"
)

loo(model_beta_birds_taxonomic_present_turn)

#loo_compare(model_beta_birds_taxonomic_present_null_turn, model_beta_birds_taxonomic_present_turn)

#Look at model summary its on 'logit' scale, because of Beta family link function
model_beta_birds_taxonomic_present_turn


## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_birds_taxonomic_present_turn, ndraws = 50) 
pp_check(model_beta_birds_taxonomic_present_turn, type = "ecdf_overlay", ndraws = 50) 

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_birds_taxonomic_present_turn$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_birds_taxonomic_present_turn, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_birds_taxonomic_present_turn, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_birds_taxonomic_present_turn, type = 'rhat') 


#back transform model summary with point estimates (Ratios)
model_parameters(model_beta_birds_taxonomic_present_turn, exponentiate = T, ci = 0.95)


#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####


#Plotting isolation effects
p1 <- birds_taxonomic %>% 
  modelr::data_grid(distance_km_sc = seq_range(distance_km_sc, n = 100) %>% as.numeric,
                    diff_emerse_pres_sc = median(diff_emerse_pres_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_present_air = median(diff_present_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_birds_taxonomic_present_turn, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = distance_km_sc)) +
  # geom_point(data = corals_taxonomic$beta_sim_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Distance between islands", y = "Beta diversity")

p2 <- birds_taxonomic %>% 
  modelr::data_grid(distance_km_sc = median(distance_km_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_pres_sc = seq_range(diff_emerse_pres_sc, n = 10),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_present_air = median(diff_present_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_birds_taxonomic_present_turn, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_emerse_pres_sc)) +
  # geom_point(data = corals_taxonomic$beta_sim_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Emersed Island Area", y = "Beta diversity")

p3 <- birds_taxonomic %>% 
  modelr::data_grid(distance_km_sc = median(distance_km_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_pres_sc = median(diff_emerse_pres_sc, na.rm = T),
                    diff_island_age_sc = seq_range(diff_island_age_sc, n = 10),
                    diff_present_air = median(diff_present_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_birds_taxonomic_present_turn, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_island_age_sc)) +
  # geom_point(data = corals_taxonomic$beta_sim_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Island Age", y = "Beta diversity")

p4 <- birds_taxonomic %>% 
  modelr::data_grid(distance_km_sc = median(distance_km_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_pres_sc = median(diff_emerse_pres_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_present_air = seq_range(diff_present_air, n = 10)) %>% 
  add_epred_draws(model_beta_birds_taxonomic_present_turn, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_present_air)) +
  # geom_point(data = corals_taxonomic$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Air Temperature", y = "Beta diversity")

p1 / p2 / p3 / p4




#8. Terrestrial - BIRDS BETA FUNCTIONAL ----


model_birds_functional <- model_birds_functional %>%
  filter(!Island1 %in% "Farquhar_Atoll",
         !Island2 %in% "Farquhar_Atoll")



# Scaling variables - difference between archipelagos islands
birds_functional <- model_birds_functional %>%
  mutate(
    beta_sim_adj = beta_turn - 0.000001,
    beta_sim_adj = ifelse(beta_sim_adj <= 0, 1e-6, beta_sim_adj),
    d_past_sc = scale(d_past.y), #emerse distance past
    diff_emerse_past_sc = scale(diff_emerse_past), #emerse island area present
    diff_emerse_pres_sc = scale(diff_emerse_pres), #emerse island area past
    diff_island_age_sc = scale(diff_island_age),
    distance_km_sc = scale(Distance_km), #emerse distance present
    diff_present_air = scale(diif_present_air),
    diff_past_air = scale(diff_past_air)
  )




#------------------------- Models

# Null

model_beta_birds_functional_present_null_turn <- brm(
  beta_sim_adj ~ 1 +
    (1 | Archipelago),
  data = birds_functional,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 5000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.98, max_treedepth = 14),
)


model_beta_birds_functional_present_null_turn <- add_criterion(model_beta_birds_functional_present_null_turn, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_birds_functional_present_null_turn,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_birds_functional_present_null_turn.rds"
)

# Assess model validity using LOO
loo(model_beta_birds_functional_present_null_turn)



# Present

model_beta_birds_functional_present_turn <- brm(
  beta_sim_adj ~  diff_emerse_pres_sc + diff_island_age_sc + distance_km_sc + diff_present_air +
    (1 | Archipelago),
  data = birds_functional,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_birds_functional_present_turn <- add_criterion(model_beta_birds_functional_present_turn, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_birds_functional_present_turn,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_birds_functional_present_turn.rds"
)

loo(model_beta_birds_functional_present_turn)

#loo_compare(model_beta_birds_functional_present_turn, model_beta_birds_functional_present_null_turn)

#Look at model summary its on 'logit' scale, because of Beta family link function
model_beta_birds_functional_present_turn

## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_birds_functional_present_turn, ndraws = 50) 
pp_check(model_beta_birds_functional_present_turn, type = "ecdf_overlay", ndraws = 50) 

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_birds_functional_present_turn$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_birds_functional_present_turn, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_birds_functional_present_turn, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_birds_functional_present_turn, type = 'rhat') 


#back transform model summary with point estimates (Ratios)
model_parameters(model_beta_birds_functional_present_turn, exponentiate = T, ci = 0.95)


#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####


#Plotting isolation effects
p1 <- birds_functional %>% 
  modelr::data_grid(distance_km_sc = seq_range(distance_km_sc, n = 100) %>% as.numeric,
                    diff_emerse_pres_sc = median(diff_emerse_pres_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_present_air = median(diff_present_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_birds_functional_present_turn, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = distance_km_sc)) +
  # geom_point(data = corals_functional$beta_sim_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Distance between islands", y = "Beta diversity")

p2 <- birds_functional %>% 
  modelr::data_grid(distance_km_sc = median(distance_km_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_pres_sc = seq_range(diff_emerse_pres_sc, n = 10),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_present_air = median(diff_present_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_birds_functional_present_turn, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_emerse_pres_sc)) +
  # geom_point(data = corals_functional$beta_sim_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Emersed Island Area", y = "Beta diversity")

p4 <- birds_functional %>% 
  modelr::data_grid(distance_km_sc = median(distance_km_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_pres_sc = median(diff_emerse_pres_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_present_air = seq_range(diff_present_air, n = 10)) %>% 
  add_epred_draws(model_beta_birds_functional_present_turn, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_present_air)) +
  # geom_point(data = corals_taxonomic$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Air Temperature", y = "Beta diversity")

p1 / p2 / p3 / p4




# TURNOVER PAST ----

#-------------------------------------------------------------------------------


#1. Marine - CORALS BETA TAXONOMIC ----


## Scaling variables - difference between marine regions ----
corals_taxonomic <- model_corals_taxonomic %>%
  mutate(
    beta_sim_adj = beta_sim - 0.000001,
    beta_sim_adj = ifelse(beta_sim_adj <= 0, 1e-6, beta_sim_adj),
    diff_isolation_sc = scale(diff_isolation)[,1],
    diff_reef_area_sc = scale(diff_reef_area)[,1],
    diff_past_isolation_sc = scale(diff_past_isolation)[,1], #quaternary distance
    diff_past_reef_area_sc = scale(diff_past_reef_area)[,1],
    diff_age_sc = scale(diff_age)[,1],
    dist_sc = scale(Distance)[,1],
    diff_present_sst = scale(diff_sst)[,1],
    diff_past_sst = scale(diff_past_sst)[,1]
  )

corals_taxonomic %>% glimpse

#------------------------- Models

# past


# Bayes models
# Null

# model_beta_corals_taxonomic_past_null_turn_turn_gp3 <- brm(
#   beta_sim_adj ~ 1 + (1|Group1) + (1|Group2) + (1|key),
#   data = corals_taxonomic,
#   family = beta_family(),
#   chains = 4,
#   cores = 4,
#   iter = 4000,
#   warmup = 3000,
#   control = list(adapt_delta = 0.99, max_treedepth = 14),
#   save_pars = save_pars(all = T)
# )

# model_beta_corals_taxonomic_past_null_turn_turn_gp3 <- add_criterion(model_beta_corals_taxonomic_past_null_turn_turn_gp3, "loo", save_psis = TRUE, reloo = T)


# model_beta_corals_taxonomic_past_null_turn <- brm(
#   beta_sim_adj ~ 1,
#   data = corals_taxonomic,
#   family = beta_family(),
#   chains = 4,
#   cores = 4,
#   iter = 4000,
#   warmup = 3000,
#   control = list(adapt_delta = 0.99, max_treedepth = 14),
#   save_pars = save_pars(all = T)
# )

# model_beta_corals_taxonomic_past_null_turn <- add_criterion(model_beta_corals_taxonomic_past_null_turn, "loo", save_psis = TRUE, reloo = T)
# 
# loo_compare(model_beta_corals_taxonomic_past_null_turn, model_beta_corals_taxonomic_past_null_turn_turn_gp2, model_beta_corals_taxonomic_past_null_turn_turn_gp1, model_beta_corals_taxonomic_past_null_turn_turn_gp3)

model_beta_corals_taxonomic_past_null_turn_gp2 <- brm(
  beta_sim_adj ~ 1 + 
    (1+dist_sc|Group1) + (1+dist_sc|Group2),
  data = corals_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  control = list(adapt_delta = 0.99, max_treedepth = 14),
  save_pars = save_pars(all = T)
)

model_beta_corals_taxonomic_past_null_turn_gp2 <- add_criterion(model_beta_corals_taxonomic_past_null_turn_gp2, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_corals_taxonomic_past_null_turn_gp2,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_corals_taxonomic_past_null_turn_gp2.rds"
)

# Assess model validity using LOO
loo(model_beta_corals_taxonomic_past_null_turn_gp2)

# past
model_beta_corals_taxonomic_past_bayes_turn <- brm(
  beta_sim_adj ~ diff_past_isolation_sc + diff_past_reef_area_sc + diff_age_sc + dist_sc + diff_past_sst +
    (1+dist_sc|Group1) + (1+dist_sc|Group2),
  data = corals_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_corals_taxonomic_past_bayes_turn <- add_criterion(model_beta_corals_taxonomic_past_bayes_turn, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_corals_taxonomic_past_bayes_turn,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_corals_taxonomic_past_bayes_turn.rds"
)

#Look at model summary its on 'logit' scale, because of bernouli family link function
model_beta_corals_taxonomic_past_bayes_turn

# Assess model validity using LOO
loo(model_beta_corals_taxonomic_past_bayes_turn)


#loo_compare(model_beta_corals_taxonomic_past_bayes_turn, model_beta_corals_taxonomic_past_null_turn_gp2)
## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_corals_taxonomic_past_bayes_turn , ndraws = 50) 
pp_check(model_beta_corals_taxonomic_past_bayes_turn, type = "ecdf_overlay", ndraws = 50) 

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_corals_taxonomic_past_bayes_turn$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_corals_taxonomic_past_bayes_turn, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_corals_taxonomic_past_bayes_turn, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_corals_taxonomic_past_bayes_turn, type = 'rhat') 


#back transform model summary with point estimates - Now Odds ratios (3:1 odds of bleaching)
model_parameters(model_beta_corals_taxonomic_past_bayes_turn, exponentiate = T, ci = 0.95)



#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####

#Plotting isolation effects
# 
p1 <- corals_taxonomic %>% 
  modelr::data_grid(
    diff_past_isolation_sc = seq_range(diff_isolation_sc, n = 10),
    diff_past_reef_area_sc = mean(diff_reef_area_sc),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_past_sst = mean(diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_taxonomic_past_bayes_turn,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_past_isolation_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in isolation", y = "Beta diversity")

#
p2 <- corals_taxonomic %>% 
  modelr::data_grid(
    diff_past_isolation_sc = seq_range(diff_isolation_sc, n=10),
    diff_past_reef_area_sc = mean(diff_reef_area_sc),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_past_sst = mean(diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_taxonomic_past_bayes_turn,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_past_isolation_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in area", y = "Beta diversity")

# 
p3 <- corals_taxonomic %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_isolation_sc),
    diff_past_reef_area_sc = mean(diff_reef_area_sc),
    diff_age_sc = seq_range(diff_age_sc, n = 10),
    dist_sc = mean(dist_sc),
    diff_past_sst = mean(diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_taxonomic_past_bayes_turn,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_age_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in age", y = "Beta diversity")

# 
p4 <- corals_taxonomic %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_past_isolation_sc),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    dist_sc = seq_range(dist_sc, n = 10),
    diff_age_sc = mean(diff_age_sc),
    diff_past_sst = mean(diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_taxonomic_past_bayes_turn,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = dist_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch distance", y = "Beta diversity")

p5 <- corals_taxonomic %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_past_isolation_sc),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    dist_sc = mean(dist_sc),
    diff_age_sc = mean(diff_age_sc),
    diff_past_sst = seq_range(diff_past_sst, n = 10)
  ) %>% 
  add_epred_draws(model_beta_corals_taxonomic_past_bayes_turn,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_past_sst)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch SST", y = "Beta diversity")


p1 / p2 / p3 / p4 / p5






#2. Marine - CORALS BETA FUNCTIONAL ----



## Scaling variables - difference between marine regions ----
corals_functional <- model_corals_functional %>%
  mutate(
    beta_sim_adj = beta_sim - 0.000001,  
    beta_sim_adj = ifelse(beta_sim_adj <= 0, 1e-6, beta_sim_adj),
    diff_isolation_sc = scale(diff_isolation)[,1],
    diff_reef_area_sc = scale(diff_reef_area)[,1],
    diff_past_isolation_sc = scale(diff_past_isolation)[,1], #quaternary distance
    diff_past_reef_area_sc = scale(diff_past_reef_area)[,1],
    diff_age_sc = scale(diff_age)[,1],
    dist_sc = scale(Distance)[,1],
    diff_present_sst = scale(diff_sst)[,1],
    diff_past_sst = scale(diff_past_sst)[,1]
  )

corals_functional %>% glimpse

#------------------------- Models


# Bayes models
# Null

model_beta_corals_functional_past_null_turn_gp2 <- brm(
  beta_sim_adj ~ 1 + 
    (1+dist_sc|Group1) + (1+dist_sc|Group2),
  data = corals_functional,
  family = beta_family(),
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 3000,
  control = list(adapt_delta = 0.99, max_treedepth = 14),
  save_pars = save_pars(all = T)
)

model_beta_corals_functional_past_null_turn_gp2 <- add_criterion(model_beta_corals_functional_past_null_turn_gp2, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_corals_functional_past_null_turn_gp2,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_corals_functional_past_null_turn_gp2.rds"
)

# Assess model validity using LOO
loo(model_beta_corals_functional_past_null_turn_gp2)

# past
model_beta_corals_functional_past_bayes_turn <- brm(
  beta_sim_adj ~ diff_past_isolation_sc + diff_past_reef_area_sc + diff_age_sc + dist_sc + diff_past_sst +
  (1+dist_sc|Group1) + (1+dist_sc|Group2),
  data = corals_functional,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_corals_functional_past_bayes_turn <- add_criterion(model_beta_corals_functional_past_bayes_turn, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_corals_functional_past_bayes_turn,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_corals_functional_past_bayes_turn.rds"
)

# Assess model validity using LOO
loo(model_beta_corals_functional_past_bayes_turn)

#Look at model summary its on 'logit' scale, because of bernouli family link function
model_beta_corals_functional_past_bayes_turn

#loo_compare(model_beta_corals_functional_past_bayes_turn, model_beta_corals_functional_past_null_turn_gp2)

## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_corals_functional_past_null_turn_gp2 , ndraws = 50) 
pp_check(model_beta_corals_functional_past_bayes_turn, type = "ecdf_overlay", ndraws = 50) 

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_corals_functional_past_bayes_turn$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_corals_functional_past_bayes_turn, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_corals_functional_past_bayes_turn, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_corals_functional_past_bayes_turn, type = 'rhat') 


#back transform model summary with point estimates - Now Odds ratios (3:1 odds of bleaching)
model_parameters(model_beta_corals_functional_past_bayes_turn, exponentiate = T, ci = 0.95)


#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####

#Plotting isolation effects
p1 <- corals_functional %>% 
  modelr::data_grid(
    diff_past_isolation_sc = seq_range(diff_past_isolation_sc, n = 10),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_past_sst = median(diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_functional_past_bayes_turn,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_past_isolation_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in isolation", y = "Beta diversity")

#
p2 <- corals_functional %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_past_isolation_sc),
    diff_past_reef_area_sc = seq_range(diff_past_reef_area_sc, n = 10),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_past_sst = median(diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_functional_past_bayes_turn,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_past_reef_area_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in reef area", y = "Beta diversity")

# 
p3 <- corals_functional %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_past_isolation_sc),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    diff_age_sc = seq_range(diff_age_sc, n = 10),
    dist_sc = mean(dist_sc),
    diff_past_sst = median(diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_functional_past_bayes_turn,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_age_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in age", y = "Beta diversity")

# 
p4 <- corals_functional %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_past_isolation_sc),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = seq_range(dist_sc, n = 10),
    diff_past_sst = median(diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_corals_functional_past_bayes_turn,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = dist_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch distance", y = "Beta diversity")

p5 <- corals_functional %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_past_isolation_sc),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    dist_sc = mean(dist_sc),
    diff_age_sc = mean(diff_age_sc),
    diff_past_sst = seq_range(diff_past_sst, n = 10)
  ) %>% 
  add_epred_draws(model_beta_corals_functional_past_bayes_turn,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_past_sst)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch SST", y = "Beta diversity")


p1 / p2 / p3 / p4 / p5





#3. Marine - FISH BETA TAXONOMIC ----


## Scaling variables - difference between marine regions ----
fish_taxonomic <- model_fish_taxonomic %>%
  mutate(
    beta_sim_adj = beta_sim - 0.000001,
    beta_sim_adj = ifelse(beta_sim_adj <= 0, 1e-6, beta_sim_adj),
    diff_isolation_sc = scale(diff_isolation)[,1],
    diff_reef_area_sc = scale(diff_reef_area)[,1],
    diff_past_isolation_sc = scale(diff_past_isolation)[,1], #quaternary distance
    diff_past_reef_area_sc = scale(diff_past_reef_area)[,1],
    diff_age_sc = scale(diff_age)[,1],
    dist_sc = scale(Distance)[,1],
    diff_present_sst = scale(diff_sst)[,1],
    diff_past_sst = scale(diff_past_sst)[,1]
  )

fish_taxonomic %>% glimpse



#------------------------- Models


# Bayes models
# Null


# model_beta_fish_taxonomic_past_null_turn_turn_gp3 <- add_criterion(model_beta_fish_taxonomic_past_null_turn_turn_gp3, "loo", save_psis = TRUE, reloo = T)
# model_beta_fish_taxonomic_past_null_turn_turn_gp2 <- brm(
#   beta_sim_adj ~ 1 + (1|Group1) + (1|Group2) ,
#   data = fish_taxonomic,
#   family = beta_family(),
#   chains = 4,
#   cores = 4,
#   iter = 4000,
#   warmup = 3000,
#   control = list(adapt_delta = 0.99, max_treedepth = 14),
#   save_pars = save_pars(all = T)
# )

model_beta_fish_taxonomic_past_null_turn_gp2 <- brm(
  beta_sim_adj ~ 1 + (1+dist_sc|Group1) + (1+dist_sc|Group2),
  data = fish_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  control = list(adapt_delta = 0.99, max_treedepth = 14),
  save_pars = save_pars(all = T)
)

model_beta_fish_taxonomic_past_null_turn_gp2 <- add_criterion(model_beta_fish_taxonomic_past_null_turn_gp2, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_fish_taxonomic_past_null_turn_gp2,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_fish_taxonomic_past_null_turn_gp2.rds"
)


# Assess model validity using LOO
loo(model_beta_fish_taxonomic_past_null_turn_gp2)

#Look at model summary its on 'logit' scale, because of bernouli family link function
model_beta_fish_taxonomic_past_null_turn_gp2


# past
model_beta_fish_taxonomic_past_bayes_turn <- brm(
  beta_sim_adj ~ diff_past_isolation_sc + diff_past_reef_area_sc + diff_age_sc + dist_sc + diff_past_sst +
  (1+dist_sc|Group1) + (1+dist_sc|Group2),
  data = fish_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_fish_taxonomic_past_bayes_turn <- add_criterion(model_beta_fish_taxonomic_past_bayes_turn, "loo", save_psis = TRUE, reloo = T)

model_beta_fish_taxonomic_past_bayes_turn

saveRDS(
  model_beta_fish_taxonomic_past_bayes_turn,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_fish_taxonomic_past_bayes_turn.rds"
)

# Assess model validity using LOO
loo(model_beta_fish_taxonomic_past_bayes_turn)

#loo_compare(model_beta_fish_taxonomic_past_bayes_turn, model_beta_fish_taxonomic_past_null_turn_gp2)
## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_fish_taxonomic_past_bayes_turn , ndraws = 50) 
pp_check(model_beta_fish_taxonomic_past_bayes_turn, type = "ecdf_overlay", ndraws = 50) 

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_fish_taxonomic_past_bayes_turn$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_fish_taxonomic_past_bayes_turn, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_fish_taxonomic_past_bayes_turn, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_fish_taxonomic_past_bayes_turn, type = 'rhat') 


#back transform model summary with point estimates - Now Odds ratios (3:1 odds of bleaching)
model_parameters(model_beta_fish_taxonomic_past_bayes_turn, exponentiate = T, ci = 0.95)


#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####

#Plotting isolation effects
p1 <- fish_taxonomic %>% 
  modelr::data_grid(
    diff_past_isolation_sc = seq_range(diff_past_isolation_sc, n = 10),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_past_sst = median(diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_taxonomic_past_bayes_turn,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_past_isolation_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in isolation", y = "Beta diversity")

#
p2 <- fish_taxonomic %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_past_isolation_sc),
    diff_past_reef_area_sc = seq_range(diff_past_reef_area_sc, n = 10),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_past_sst = median(diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_taxonomic_past_bayes_turn,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_past_reef_area_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in reef area", y = "Beta diversity")

# 
p3 <- fish_taxonomic %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_past_isolation_sc),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    diff_age_sc = seq_range(diff_age_sc, n = 10),
    dist_sc = mean(dist_sc),
    diff_past_sst = median(diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_taxonomic_past_bayes_turn,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_age_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in age", y = "Beta diversity")

# 
p4 <- fish_taxonomic %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_past_isolation_sc),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    dist_sc = seq_range(dist_sc, n = 10),
    diff_age_sc = mean(diff_age_sc),
    diff_past_sst = median(diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_taxonomic_past_bayes_turn,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = dist_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch distance", y = "Beta diversity")

p5 <- fish_taxonomic %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_past_isolation_sc),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    dist_sc = mean(dist_sc),
    diff_age_sc = mean(diff_age_sc),
    diff_past_sst = seq_range(diff_past_sst, n = 10)
  ) %>% 
  add_epred_draws(model_beta_fish_taxonomic_past_bayes_turn,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_past_sst)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch SST", y = "Beta diversity")


p1 / p2 / p3 / p4 / p5




#4. Marine - FISH BETA FUNCTIONAL ----


## Scaling variables - difference between marine regions ----
fish_functional <- model_fish_functional %>%
  mutate(
    beta_sim_adj = beta_sim - 0.000001,  
    beta_sim_adj = ifelse(beta_sim_adj <= 0, 1e-6, beta_sim_adj),
    diff_isolation_sc = scale(diff_isolation)[,1],
    diff_reef_area_sc = scale(diff_reef_area)[,1],
    diff_past_isolation_sc = scale(diff_past_isolation)[,1], #quaternary distance
    diff_past_reef_area_sc = scale(diff_past_reef_area)[,1],
    diff_age_sc = scale(diff_age)[,1],
    dist_sc = scale(Distance)[,1],
    diff_present_sst = scale(diff_sst)[,1],
    diff_past_sst = scale(diff_past_sst)[,1]
  ) 

fish_functional %>% glimpse

#------------------------- Models


# Bayes models
# Null
model_beta_fish_functional_past_null_turn_gp2 <- brm(
  beta_sim_adj ~ 1 + (1+dist_sc|Group1) + (1+dist_sc|Group2) ,
  data = fish_functional,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  control = list(adapt_delta = 0.99, max_treedepth = 14),
  save_pars = save_pars(all = T)
)

model_beta_fish_functional_past_null_turn_gp2 <- add_criterion(model_beta_fish_functional_past_null_turn_gp2, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_fish_functional_past_null_turn_gp2,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_fish_functional_past_null_turn_gp2.rds"
)

# Assess model validity using LOO
loo(model_beta_fish_functional_past_null_turn_gp2)


# past
model_beta_fish_functional_past_bayes_turn <- brm(
  beta_sim_adj ~  1 + diff_past_isolation_sc + diff_past_reef_area_sc + diff_age_sc + dist_sc + diff_past_sst +
  (1+dist_sc|Group1) + (1+dist_sc|Group2),
  data = fish_functional,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)



# Add Leave One Out cross validation criterion
model_beta_fish_functional_past_bayes_turn <- add_criterion(model_beta_fish_functional_past_bayes_turn, "loo", save_psis = TRUE, reloo = T, moment_match = TRUE)

saveRDS(
  model_beta_fish_functional_past_bayes_turn,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_fish_functional_past_bayes_turn.rds"
)

# Assess model validity using LOO
loo(model_beta_fish_functional_past_bayes_turn)

#loo_compare(model_beta_fish_functional_past_bayes_turn, model_beta_fish_functional_past_null_turn_gp2)

## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_fish_functional_past_bayes_turn , ndraws = 50) 
pp_check(model_beta_fish_functional_past_bayes_turn, type = "ecdf_overlay", ndraws = 50) 
pp_check(model_beta_fish_functional_past_bayes_turn , ndraws = 50, type = "dens_overlay_grouped", group = "Group1")
pp_check(model_beta_fish_functional_past_bayes_turn , ndraws = 50, type = "scatter_avg_grouped", group = "Group1")
pp_check(model_beta_fish_functional_past_bayes_turn , ndraws = 50, type = "error_scatter_avg_grouped", group = "Group1")
pp_check(model_beta_fish_functional_past_bayes_turn, type='pit_ecdf_grouped', group = "Group1", prob = 0.95, plot_diff = F)
pp_check(model_beta_fish_functional_past_bayes_turn, type='pit_ecdf_grouped', group = "Group2", prob = 0.95, plot_diff = F)

pp_check(model_beta_fish_functional_past_bayes_turn , ndraws = 50, type = "scatter_avg_grouped", group = "Group2")

pairs(model_beta_fish_functional_past_bayes_turn, variable = variables(model_beta_fish_functional_past_bayes_turn)[1:5])
fish_functional %>% select(Group1, Group2) %>% unique %>% view

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_fish_functional_past_bayes_turn$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_fish_functional_past_bayes_turn, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_fish_functional_past_bayes_turn, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_fish_functional_past_bayes_turn, type = 'rhat') 


#back transform model summary with point estimates - Now Odds ratios (3:1 odds of bleaching)
model_parameters(model_beta_fish_functional_past_bayes_turn, exponentiate = T, ci = 0.95)



#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####

#Plotting isolation effects
p1 <- fish_functional %>% 
  modelr::data_grid(
    diff_past_isolation_sc = seq_range(diff_past_isolation_sc, n = 10),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_past_sst = median(diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_functional_past_bayes_turn,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_past_isolation_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in isolation", y = "Beta diversity")

#
p2 <- fish_functional %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_past_isolation_sc),
    diff_past_reef_area_sc = seq_range(diff_past_reef_area_sc, n = 10),
    diff_age_sc = mean(diff_age_sc),
    dist_sc = mean(dist_sc),
    diff_past_sst = median(diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_functional_past_bayes_turn,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_past_reef_area_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in reef area", y = "Beta diversity")

# 
p3 <- fish_functional %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_past_isolation_sc),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    diff_age_sc = seq_range(diff_age_sc, n = 10),
    dist_sc = mean(dist_sc),
    diff_past_sst = median(diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_functional_past_bayes_turn,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_age_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in age", y = "Beta diversity")

# 
p4 <- fish_functional %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_past_isolation_sc),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    dist_sc = seq_range(dist_sc, n = 10),
    diff_age_sc = mean(diff_age_sc),
    diff_past_sst = median(diff_past_sst)
  ) %>% 
  add_epred_draws(model_beta_fish_functional_past_bayes_turn,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = dist_sc)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch distance", y = "Beta diversity")

p5 <- fish_functional %>% 
  modelr::data_grid(
    diff_past_isolation_sc = mean(diff_past_isolation_sc),
    diff_past_reef_area_sc = mean(diff_past_reef_area_sc),
    dist_sc = mean(dist_sc),
    diff_age_sc = mean(diff_age_sc),
    diff_past_sst = seq_range(diff_past_sst, n = 10)
  ) %>% 
  add_epred_draws(model_beta_fish_functional_past_bayes_turn,
                  re_formula = NA, ndraws = 500) %>% 
  unique %>% 
  ggplot(aes(y = .epred, x = diff_past_sst)) +
  geom_line(aes(group = .draw), alpha = .08) +
  theme_tidybayes() +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Difference in Arch SST", y = "Beta diversity")


p1 / p2 / p3 / p4 / p5





#5. Terrestrial - PLANTS BETA TAXONOMIC ----


model_plants_taxonomic <- model_plants_taxonomic %>%
  filter(!Island1 %in% "Farquhar_Atoll",
         !Island2 %in% "Farquhar_Atoll")


# Scaling variables - difference between archipelagos islands
plants_taxonomic <- model_plants_taxonomic %>%
  mutate(
    beta_sim_adj = beta_sim - 0.000001,
    beta_sim_adj = ifelse(beta_sim_adj <= 0, 1e-6, beta_sim_adj),
    d_past_sc = scale(d_past.y), #emerse distance past
    diff_emerse_past_sc = scale(diff_emerse_past), #emerse island area past
    diff_emerse_pres_sc = scale(diff_emerse_pres), #emerse island area past
    diff_island_age_sc = scale(diff_island_age),
    distance_km_sc = scale(Distance_km), #emerse distance past
    diff_present_air = scale(diif_present_air),
    diff_past_air = scale(diff_past_air)
  )

# scaled_diff_emerse_pres <- plants_taxonomic$diff_emerse_pres - mean(plants_taxonomic$diff_emerse_pres)/sd(plants_taxonomic$diff_emerse_pres)
# 
# diff_emerse_pres <- scaled_diff_emerse_pres * sd(plants_taxonomic$diff_emerse_pres) + mean(plants_taxonomic$diff_emerse_pres)
# 
# diff_emerse_pres_sc <-  scale(plants_taxonomic$diff_island_age_sc)
# diff_emerse_pres_sc %>% glimpse


#------------------------- Models

# Null

model_beta_plants_taxonomic_past_null_turn <- brm(
  beta_sim_adj ~ 1 +
    (1 | Archipelago),
  data = plants_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 5000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.98, max_treedepth = 14),
)


model_beta_plants_taxonomic_past_null_turn <- add_criterion(model_beta_plants_taxonomic_past_null_turn, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_plants_taxonomic_past_null_turn,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_plants_taxonomic_past_null_turn.rds"
)

# Assess model validity using LOO
loo(model_beta_plants_taxonomic_past_null_turn)


# past

model_beta_plants_taxonomic_past_turn <- brm(
  beta_sim_adj ~  diff_emerse_past_sc + diff_island_age_sc + d_past_sc + diff_past_air + 
    (1 | Archipelago),
  data = plants_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_plants_taxonomic_past_turn <- add_criterion(model_beta_plants_taxonomic_past_turn, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_plants_taxonomic_past_turn,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_plants_taxonomic_past_turn.rds"
)


loo(model_beta_plants_taxonomic_past_turn)

#loo_compare(model_beta_plants_taxonomic_past_turn, model_beta_plants_taxonomic_past_null_turn)

#Look at model summary its on 'logit' scale, because of Beta family link function
model_beta_plants_taxonomic_past_turn

## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_plants_taxonomic_past_turn, ndraws = 50) 
pp_check(model_beta_plants_taxonomic_past_turn, type = "ecdf_overlay", ndraws = 50) 

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_plants_taxonomic_past$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_plants_taxonomic_past_turn, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_plants_taxonomic_past_turn, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_plants_taxonomic_past_turn, type = 'rhat') 


#back transform model summary with point estimates (Ratios)
model_parameters(model_beta_plants_taxonomic_past_turn, exponentiate = T, ci = 0.95)


#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####


#Plotting isolation effects
p1 <- plants_taxonomic %>% 
  modelr::data_grid(d_past_sc = seq_range(d_past_sc, n = 100) %>% as.numeric,
                    diff_emerse_past_sc = median(diff_emerse_past_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_past_air = median(diff_past_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_plants_taxonomic_past_turn, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = d_past_sc)) +
  # geom_point(data = corals_taxonomic$beta_sim_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Distance between islands", y = "Beta diversity")

p2 <- plants_taxonomic %>% 
  modelr::data_grid(d_past_sc = median(d_past_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_past_sc = seq_range(diff_emerse_past_sc, n = 10),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_past_air = median(diff_past_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_plants_taxonomic_past_turn, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_emerse_past_sc)) +
  # geom_point(data = corals_taxonomic$beta_sim_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Emersed Island Area", y = "Beta diversity")

p3 <- plants_taxonomic %>% 
  modelr::data_grid(d_past_sc = median(d_past_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_past_sc = median(diff_emerse_past_sc, na.rm = T),
                    diff_island_age_sc = seq_range(diff_island_age_sc, n = 10),
                    diff_past_air = median(diff_past_air, na.rm = T)) %>% 
  add_epred_draws(model_beta_plants_taxonomic_past_turn, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_island_age_sc)) +
  # geom_point(data = corals_taxonomic$beta_sim_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Island Age", y = "Beta diversity")

p4 <- plants_taxonomic %>% 
  modelr::data_grid(d_past_sc = median(d_past_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_past_sc = median(diff_emerse_past_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_past_air = seq_range(diff_past_air, n = 10)) %>% 
  add_epred_draws(model_beta_plants_taxonomic_past_turn, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_past_air)) +
  # geom_point(data = corals_taxonomic$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Air Temperature", y = "Beta diversity")

library(cowplot)

p1 / p2 / p3 / p4






#6. Terrestrial - PLANTS BETA FUNCTIONAL ----


model_plants_functional <- model_plants_functional %>%
  filter(!Island1 %in% "Farquhar_Atoll",
         !Island2 %in% "Farquhar_Atoll")



# Scaling variables - difference between archipelagos islands
plants_functional <- model_plants_functional %>%
  mutate(
    beta_sim_adj = beta_turn - 0.000001, 
    beta_sim_adj = ifelse(beta_sim_adj <= 0, 1e-6, beta_sim_adj),
    d_past_sc = scale(d_past.y), #emerse distance past
    diff_emerse_past_sc = scale(diff_emerse_past), #emerse island area past
    diff_emerse_pres_sc = scale(diff_emerse_pres), #emerse island area past
    diff_island_age_sc = scale(diff_island_age),
    distance_km_sc = scale(Distance_km), #emerse distance past
    diff_present_air = scale(diif_present_air),
    diff_past_air = scale(diff_past_air)
  )


#------------------------- Models

# Null

model_beta_plants_functional_past_null_turn <- brm(
  beta_sim_adj ~ 1 +
    (1 | Archipelago),
  data = plants_functional,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 5000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.98, max_treedepth = 14),
)


model_beta_plants_functional_past_null_turn <- add_criterion(model_beta_plants_functional_past_null_turn, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_plants_functional_past_null_turn,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_plants_functional_past_null_turn.rds"
)

# Assess model validity using LOO
loo(model_beta_plants_functional_past_null_turn)



# past

model_beta_plants_functional_past_turn <- brm(
  beta_sim_adj ~  diff_emerse_past_sc + diff_island_age_sc + d_past_sc + diff_past_air +
    (1 | Archipelago),
  data = plants_functional,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_plants_functional_past_turn <- add_criterion(model_beta_plants_functional_past_turn, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_plants_functional_past_turn,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_plants_functional_past_turn.rds"
)

loo(model_beta_plants_functional_past_turn)

#loo_compare(model_beta_plants_functional_past_null_turn, model_beta_plants_functional_past)

#Look at model summary its on 'logit' scale, because of Beta family link function
model_beta_plants_functional_past_turn

## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_plants_functional_past_turn, ndraws = 50) 
pp_check(model_beta_plants_functional_past_turn, type = "ecdf_overlay", ndraws = 50) 

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_plants_functional_past_turn$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_plants_functional_past_turn, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_plants_functional_past_turn, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_plants_functional_past_turn, type = 'rhat') 


#back transform model summary with point estimates (Ratios)
model_parameters(model_beta_plants_functional_past_turn, exponentiate = T, ci = 0.95)


#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####


#Plotting isolation effects
p1 <- plants_functional %>% 
  modelr::data_grid(d_past_sc = seq_range(d_past_sc, n = 100) %>% as.numeric,
                    diff_emerse_past_sc = median(diff_emerse_past_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_past_air = median(diff_past_air, na.rm =T)) %>% 
  add_epred_draws(model_beta_plants_functional_past_turn, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = d_past_sc)) +
  # geom_point(data = corals_functional$beta_sim_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Distance between islands", y = "Beta diversity")

p2 <- plants_functional %>% 
  modelr::data_grid(d_past_sc = median(d_past_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_past_sc = seq_range(diff_emerse_past_sc, n = 10),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_past_air = median(diff_past_air, na.rm =T)) %>% 
  add_epred_draws(model_beta_plants_functional_past_turn, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_emerse_past_sc)) +
  # geom_point(data = corals_functional$beta_sim_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Emersed Island Area", y = "Beta diversity")

p3 <- plants_functional %>% 
  modelr::data_grid(d_past_sc = median(d_past_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_past_sc = median(diff_emerse_past_sc, na.rm = T),
                    diff_island_age_sc = seq_range(diff_island_age_sc, n = 10),
                    diff_past_air = median(diff_past_air, na.rm =T)) %>% 
  add_epred_draws(model_beta_plants_functional_past_turn, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_island_age_sc)) +
  # geom_point(data = corals_functional$beta_sim_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Island Age", y = "Beta diversity")

p4 <- plants_functional %>% 
  modelr::data_grid(d_past_sc = median(d_past_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_past_sc = median(diff_emerse_past_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_past_air = seq_range(diff_past_air, n = 10)) %>% 
  add_epred_draws(model_beta_plants_functional_past_turn, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_past_air)) +
  # geom_point(data = corals_taxonomic$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Air Temperature", y = "Beta diversity")

library(cowplot)

p1 / p2 / p3 / p4





#7. Terrestrial - BIRDS BETA TAXONOMIC ----


model_birds_taxonomic <- model_birds_taxonomic %>%
  filter(!Island1 %in% "Farquhar_Atoll",
         !Island2 %in% "Farquhar_Atoll")



# Scaling variables - difference between archipelagos islands
birds_taxonomic <- model_birds_taxonomic %>%
  mutate(
    beta_sim_adj = beta_sim - 0.000001,
    beta_sim_adj = ifelse(beta_sim_adj <= 0, 1e-6, beta_sim_adj),
    d_past_sc = scale(d_past.y), #emerse distance past
    diff_emerse_past_sc = scale(diff_emerse_past), #emerse island area past
    diff_emerse_pres_sc = scale(diff_emerse_pres), #emerse island area past
    diff_island_age_sc = scale(diff_island_age),
    distance_km_sc = scale(Distance_km), #emerse distance past
    diff_present_air = scale(diif_present_air),
    diff_past_air = scale(diff_past_air)
  )




#------------------------- Models

# Null

model_beta_birds_taxonomic_past_null_turn <- brm(
  beta_sim_adj ~ 1 +
    (1 | Archipelago),
  data = birds_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 5000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.98, max_treedepth = 14),
)


model_beta_birds_taxonomic_past_null_turn <- add_criterion(model_beta_birds_taxonomic_past_null_turn, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_birds_taxonomic_past_null_turn,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_birds_taxonomic_past_null_turn.rds"
)


# Assess model validity using LOO
loo(model_beta_birds_taxonomic_past_null_turn)



# past

model_beta_birds_taxonomic_past_turn <- brm(
  beta_sim_adj ~  diff_emerse_past_sc + diff_island_age_sc + d_past_sc + diff_past_air + 
    (1 | Archipelago),
  data = birds_taxonomic,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_birds_taxonomic_past_turn <- add_criterion(model_beta_birds_taxonomic_past_turn, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_birds_taxonomic_past_turn,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_birds_taxonomic_past_turn.rds"
)


loo(model_beta_birds_taxonomic_past_turn)

#loo_compare(model_beta_birds_taxonomic_past_null_turn, model_beta_birds_taxonomic_past_turn)

#Look at model summary its on 'logit' scale, because of Beta family link function
model_beta_birds_taxonomic_past_turn

## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_birds_taxonomic_past_turn, ndraws = 50) 

pp_check(model_beta_birds_taxonomic_past_turn, type = "ecdf_overlay", ndraws = 50) 

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_birds_taxonomic_past_turn$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_birds_taxonomic_past_turn, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_birds_taxonomic_past_turn, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_birds_taxonomic_past_turn, type = 'rhat') 


#back transform model summary with point estimates (Ratios)
model_parameters(model_beta_birds_taxonomic_past_turn, exponentiate = T, ci = 0.95)


#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####


#Plotting isolation effects
p1 <- birds_taxonomic %>% 
  modelr::data_grid(d_past_sc = seq_range(d_past_sc, n = 10) %>% as.numeric,
                    diff_emerse_past_sc = median(diff_emerse_past_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_past_air = median(diff_past_air, na.rm =T)) %>% 
  add_epred_draws(model_beta_birds_taxonomic_past_turn, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = d_past_sc)) +
  # geom_point(data = corals_taxonomic$beta_sim_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Distance between islands", y = "Beta diversity")

p2 <- birds_taxonomic %>% 
  modelr::data_grid(d_past_sc = median(d_past_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_past_sc = seq_range(diff_emerse_past_sc, n = 10),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_past_air = median(diff_past_air, na.rm =T)) %>% 
  add_epred_draws(model_beta_birds_taxonomic_past_turn, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_emerse_past_sc)) +
  # geom_point(data = corals_taxonomic$beta_sim_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Emersed Island Area", y = "Beta diversity")

p3 <- birds_taxonomic %>% 
  modelr::data_grid(d_past_sc = median(d_past_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_past_sc = median(diff_emerse_past_sc, na.rm = T),
                    diff_island_age_sc = seq_range(diff_island_age_sc, n = 10),
                    diff_past_air = median(diff_past_air, na.rm =T)) %>% 
  add_epred_draws(model_beta_birds_taxonomic_past_turn, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_island_age_sc)) +
  # geom_point(data = corals_taxonomic$beta_sim_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Island Age", y = "Beta diversity")

p4 <- birds_taxonomic %>% 
  modelr::data_grid(d_past_sc = median(d_past_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_past_sc = median(diff_emerse_past_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_past_air = seq_range(diff_past_air, n = 10)) %>% 
  add_epred_draws(model_beta_birds_taxonomic_past_turn, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_past_air)) +
  # geom_point(data = corals_taxonomic$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Air Temperature", y = "Beta diversity")

library(cowplot)

p1 / p2 / p3 / p4





#8. Terrestrial - BIRDS BETA FUNCTIONAL ----


model_birds_functional <- model_birds_functional %>%
  filter(!Island1 %in% "Farquhar_Atoll",
         !Island2 %in% "Farquhar_Atoll")



# Scaling variables - difference between archipelagos islands
birds_functional <- model_birds_functional %>%
  mutate(
    beta_sim_adj = beta_turn - 0.000001,
    beta_sim_adj = ifelse(beta_sim_adj <= 0, 1e-6, beta_sim_adj),
    d_past_sc = scale(d_past.y), #emerse distance past
    diff_emerse_past_sc = scale(diff_emerse_past), #emerse island area past
    diff_emerse_pres_sc = scale(diff_emerse_pres), #emerse island area past
    diff_island_age_sc = scale(diff_island_age),
    distance_km_sc = scale(Distance_km), #emerse distance past
    diff_present_air = scale(diif_present_air),
    diff_past_air = scale(diff_past_air)
  )
    





#------------------------- Models

# Null

model_beta_birds_functional_past_null_turn <- brm(
  beta_sim_adj ~ 1 +
    (1 | Archipelago),
  data = birds_functional,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 5000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.98, max_treedepth = 14),
)


model_beta_birds_functional_past_null_turn <- add_criterion(model_beta_birds_functional_past_null_turn, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_birds_functional_past_null_turn,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_birds_functional_past_null_turn.rds"
)

# Assess model validity using LOO
loo(model_beta_birds_functional_past_null_turn)



# past

model_beta_birds_functional_past_turn <- brm(
  beta_sim_adj ~  diff_emerse_past_sc + diff_island_age_sc + d_past_sc + diff_past_air +
    (1 | Archipelago),
  data = birds_functional,
  family = beta_family(),
  chains = chains,
  cores = chains,           
  backend = "cmdstanr",         
  threads = threading(threads_per_chain),  
  iter = 4000,
  warmup = 3000,
  save_pars = save_pars(all = T),
  control = list(adapt_delta = 0.99, max_treedepth = 14),
)


# Add Leave One Out cross validation criterion
model_beta_birds_functional_past_turn <- add_criterion(model_beta_birds_functional_past_turn, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  model_beta_birds_functional_past_turn,
  "/home/student/Documents/Luiza/Islands_Biogeography/Inter_Data/model_beta_birds_functional_past_turn.rds"
)


loo(model_beta_birds_functional_past_turn)

#loo_compare(model_beta_birds_functional_past_null_turn, model_beta_birds_functional_past_turn)

#Look at model summary its on 'logit' scale, because of Beta family link function
model_beta_birds_functional_past_turn

## Model checks ####
#posterior predictive check, general agreement in distribution
#i.e, could the model accurately produce these data?
pp_check(model_beta_birds_functional_past_turn, ndraws = 50) 
pp_check(model_beta_birds_functional_past_turn, type = "ecdf_overlay", ndraws = 50) 

#any pathological model behaviour?
rstan::check_hmc_diagnostics(model_beta_birds_functional_past_turn$fit) 

# HAIRY CATERPILLARS  - well mixed trace plot means model explored all off the data space equally
mcmc_plot(model_beta_birds_functional_past_turn, type = "trace") 

# No correlation between iterations 
mcmc_plot(model_beta_birds_functional_past_turn, type = "acf")  

#Rhat values <1.05 means chains converged on estimate of a given parameter
mcmc_plot(model_beta_birds_functional_past_turn, type = 'rhat') 


#back transform model summary with point estimates (Ratios)
model_parameters(model_beta_birds_functional_past_turn, exponentiate = T, ci = 0.95)


#--------------------------
# Plot

##Plot distributions of bleaching probability for each `Taxa`####


#Plotting isolation effects
p1 <- birds_functional %>% 
  modelr::data_grid(d_past_sc = seq_range(d_past_sc, n = 100) %>% as.numeric,
                    diff_emerse_past_sc = median(diff_emerse_past_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_past_air = median(diff_past_air, na.rm=T)) %>% 
  add_epred_draws(model_beta_birds_functional_past_turn, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = d_past_sc)) +
  # geom_point(data = corals_functional$beta_sim_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Distance between islands", y = "Beta diversity")

p2 <- birds_functional %>% 
  modelr::data_grid(d_past_sc = median(d_past_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_past_sc = seq_range(diff_emerse_past_sc, n = 10),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_past_air = median(diff_past_air, na.rm=T)) %>% 
  add_epred_draws(model_beta_birds_functional_past_turn, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_emerse_past_sc)) +
  # geom_point(data = corals_functional$beta_sim_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Emersed Island Area", y = "Beta diversity")

p3 <- birds_functional %>% 
  modelr::data_grid(d_past_sc = median(d_past_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_past_sc = median(diff_emerse_past_sc, na.rm = T),
                    diff_island_age_sc = seq_range(diff_island_age_sc, n = 10),
                    diff_past_air = median(diff_past_air, na.rm=T)) %>% 
  add_epred_draws(model_beta_birds_functional_past_turn, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_island_age_sc)) +
  # geom_point(data = corals_functional$beta_sim_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Island Age", y = "Beta diversity")

p4 <- birds_functional %>% 
  modelr::data_grid(d_past_sc = median(d_past_sc, na.rm = T) %>% as.numeric,
                    diff_emerse_past_sc = median(diff_emerse_past_sc, na.rm = T),
                    diff_island_age_sc = median(diff_island_age_sc, na.rm = T),
                    diff_past_air = seq_range(diff_past_air, n = 10)) %>% 
  add_epred_draws(model_beta_birds_functional_past_turn, re_formula = NA, ndraws = 1000) %>% unique %>% 
  ggplot(aes(y = .epred, x = diff_past_air)) +
  # geom_point(data = corals_taxonomic$beta_sor_adj, alpha = 0.2)+
  geom_line(aes(y = .epred, group = .draw), alpha = .03)+
  # geom_rug(linewidth = .01)+
  theme_tidybayes()+
  # coord_cartesian(ylim = c(0,1))+
  labs(x = "Difference Air Temperature", y = "Beta diversity")

library(cowplot)

p1 / p2 / p3 / p4



