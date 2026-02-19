
### Packages

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


#Beta diversity and variables data ----


#### TAXA * VARIABLES

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
model_plants_taxonomic <- read.csv("model_data_plants_taxonomic_arch.csv", header = TRUE) %>% 
  mutate_if(is.character, as.factor) %>% glimpse

# Plants functional -----
model_plants_functional <- read.csv("model_data_plants_functional_arch.csv", header = TRUE) %>% 
  mutate_if(is.character, as.factor) %>% glimpse

# Birds taxonomic -----
model_birds_taxonomic <- read.csv("model_data_birds_taxonomic_arch.csv", header = TRUE) %>% 
  mutate_if(is.character, as.factor) %>% glimpse

# Birds functional -----
model_birds_functional <- read.csv("model_data_bird_functional_arch.csv", header = TRUE) %>% 
  mutate_if(is.character, as.factor) %>% glimpse


###

# Parallel 

total_cores <- 150L
chains <- 4L
threads_per_chain <- floor(total_cores / chains)  # = 25
message("using chains=", chains, " threads_per_chain=", threads_per_chain,
        " total threads=", chains * threads_per_chain)



#### Beta taxonomic NEST PRESENT ----------------------------------------------------

# Corais taxonomic
coral_tax <- model_corals_taxonomic %>%
  mutate(
    Archipelago1 = Archipelago1.model,
    Archipelago2 = Archipelago2.model,
    taxa = "Coral",
    dist = Distance,      
    area = diff_reef_area,
    temp = diff_sst,
    age = diff_age,
    group2 = Group2,
    group1 = Group1
  )

# Fish taxonomic
fish_tax <- model_fish_taxonomic %>%
  mutate(
    Archipelago1 = Archipelago1,
    Archipelago2 = Archipelago2,
    taxa = "Fish",
    dist = Distance,      
    area = diff_reef_area,
    temp = diff_sst,
    age = diff_age,
    group2 = Group2,
    group1 = Group1
  )

# Plants taxonomic
plant_tax <- model_plants_taxonomic %>%
  mutate(
    Archipelago1 = Archipelago1,
    Archipelago2 = Archipelago2,
    taxa = "Plant",
    dist = Distance,      
    area = diff_present_area,
    temp = diff_present_air,
    age = diff_age,
    group2 = Group2,
    group1 = Group1
  )

# Birds taxonomic
bird_tax <- model_birds_taxonomic %>%
  mutate(
    Archipelago1 = Archipelago1,
    Archipelago2 = Archipelago2,
    taxa = "Bird",
    dist = Distance,     
    area = diff_present_area,
    temp = diff_present_air,
    age = diff_age,
    group2 = Group2,
    group1 = Group1
  )

# selecting columns
coral_tax <- coral_tax %>% select(beta_sne, Archipelago1, Archipelago2, taxa, area, dist, age, temp, group2, group1)
fish_tax  <- fish_tax  %>% select(beta_sne, Archipelago1, Archipelago2, taxa, area, dist, age, temp, group2, group1)
plant_tax <- plant_tax %>% select(beta_sne, Archipelago1, Archipelago2, taxa, area, dist, age, temp, group2, group1)
bird_tax  <- bird_tax  %>% select(beta_sne, Archipelago1, Archipelago2, taxa, area, dist, age, temp, group2, group1)

# combining 
model_all_tax <- bind_rows(coral_tax, fish_tax, plant_tax, bird_tax)

model_all_tax$taxa <- factor(model_all_tax$taxa, levels = c("Coral", "Fish", "Plant", "Bird"))

###

# Scaling up variables per ecosystem type

# 0 and 1
rescale_01 <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}


model_all_tax <- model_all_tax %>%
  mutate(
    ecosystem = case_when(
      taxa %in% c("Coral", "Fish")  ~ "Marine",
      taxa %in% c("Plant", "Bird")  ~ "Terrestrial"
    )
  )

model_all_tax$ecosystem <- factor(
  model_all_tax$ecosystem,
  levels = c("Marine", "Terrestrial")
)


##

model_all_tax_scaled <- model_all_tax %>%
  
  # global variables
  mutate(
    age  = rescale_01(age),
    dist = rescale_01(dist)
  ) %>%
  
  # ecosystem variables
  group_by(ecosystem) %>%
  mutate(
    area = rescale_01(area),
    temp = rescale_01(temp)
  ) %>%
  ungroup()


### Scaling again

epsilon <- 1e-4  # um valor pequeno

model_all_tax_scaled <- model_all_tax_scaled %>%
  mutate(
    beta_sne_rescaled = beta_sne * (1 - 2 * epsilon) + epsilon
  )

##

mod_tax_bayes_model2 <- brm(
  beta_sne_rescaled ~ area * taxa + dist * taxa + age * taxa + temp * taxa + (1 | group2),
  data = model_all_tax_scaled,
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


##
mod_tax_bayes_model2 <- add_criterion(mod_tax_bayes_model2, "loo", save_psis = TRUE, reloo = T)


saveRDS(
  mod_tax_bayes_model2,
  "/home/student/Islands_Biogeography/Inter_Data/model_beta_all_taxonomic_nest_present.rds"
)


# Assess model validity using LOO
loo(mod_tax_bayes_model2)


### Plot

post_env <- mod_tax_bayes_model2 %>%
  gather_draws(
    b_area,
    `b_area:taxaFish`, `b_area:taxaPlant`, `b_area:taxaBird`,
    b_dist,
    `b_taxaFish:dist`, `b_taxaPlant:dist`, `b_taxaBird:dist`,
    b_age,
    `b_taxaFish:age`,  `b_taxaPlant:age`,  `b_taxaBird:age`,
    b_temp,
    `b_taxaFish:temp`, `b_taxaPlant:temp`, `b_taxaBird:temp`
  )


#

post_env_eff <- post_env %>%
  mutate(
    Variable = case_when(
      grepl("area", .variable) ~ "Area",
      grepl("dist", .variable) ~ "Distance",
      grepl("age",  .variable) ~ "Age",
      grepl("temp", .variable) ~ "Temperature"
    )
  ) %>%
  group_by(.draw, Variable) %>%
  summarise(
    Coral = sum(.value[.variable %in% c("b_area", "b_dist", "b_age", "b_temp")]),
    Fish  = Coral + sum(.value[grepl("taxaFish",  .variable)]),
    Plant = Coral + sum(.value[grepl("taxaPlant", .variable)]),
    Bird  = Coral + sum(.value[grepl("taxaBird",  .variable)]),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(Coral, Fish, Plant, Bird),
    names_to = "Taxon",
    values_to = "Effect"
  )


#

tax_nest_pres <-  ggplot(
  post_env_eff,
  aes(x = Effect, fill = Taxon)
) +
  stat_halfeye(
    adjust = 0.6,
    alpha = 0.6,
    .width = c(0.5, 0.8, 0.95),
    slab_color = NA
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_grid(Variable ~ .) +
  coord_cartesian(xlim = c(-4, 10)) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    strip.text.y = element_text(angle = 0, face = "bold"),
    legend.position = "top"
  ) +
  labs(
    x = "Posterior effect size on β-diversity nest",
    y = NULL,
    fill = "Taxonomic group"
  )




#### Beta taxonomic TURN PRESENT ----------------------------------------------------

# coral

coral_turn_pres <- model_corals_taxonomic %>%
  mutate(
    Archipelago1 = Archipelago1.model,
    Archipelago2 = Archipelago2.model,
    taxa = "Coral",
    dist = Distance,
    area = diff_reef_area,
    temp = diff_sst,
    age  = diff_age,
    group2 = Group2,
    group1 = Group1
  )


# fish

fish_turn_pres <- model_fish_taxonomic %>%
  mutate(
    Archipelago1 = Archipelago1,
    Archipelago2 = Archipelago2,
    taxa = "Fish",
    dist = Distance,
    area = diff_reef_area,
    temp = diff_sst,
    age  = diff_age,
    group2 = Group2,
    group1 = Group1
  )


# plant

plant_turn_pres <- model_plants_taxonomic %>%
  mutate(
    Archipelago1 = Archipelago1,
    Archipelago2 = Archipelago2,
    taxa = "Plant",
    dist = Distance,
    area = diff_present_area,
    temp = diff_present_air,
    age  = diff_age,
    group2 = Group2,
    group1 = Group1
  )


# bird

bird_turn_pres <- model_birds_taxonomic %>%
  mutate(
    Archipelago1 = Archipelago1,
    Archipelago2 = Archipelago2,
    taxa = "Bird",
    dist = Distance,
    area = diff_present_area,
    temp = diff_present_air,
    age  = diff_age,
    group2 = Group2,
    group1 = Group1
  )


###

coral_turn_pres <- coral_turn_pres %>% select(beta_sim, Archipelago1, Archipelago2, taxa, area, dist, age, temp, group2, group1)
fish_turn_pres  <- fish_turn_pres  %>% select(beta_sim, Archipelago1, Archipelago2, taxa, area, dist, age, temp, group2, group1)
plant_turn_pres <- plant_turn_pres %>% select(beta_sim, Archipelago1, Archipelago2, taxa, area, dist, age, temp, group2, group1)
bird_turn_pres  <- bird_turn_pres  %>% select(beta_sim, Archipelago1, Archipelago2, taxa, area, dist, age, temp, group2, group1)


# merging data

model_all_tax_turn_pres <- bind_rows(coral_turn_pres, fish_turn_pres, plant_turn_pres, bird_turn_pres)

model_all_tax_turn_pres$taxa <- factor(
  model_all_tax_turn_pres$taxa,
  levels = c("Coral", "Fish", "Plant", "Bird")
)


###

# Scaling up variables per ecosystem type

# 0 and 1
rescale_01 <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}


model_all_tax <- model_all_tax_turn_pres %>%
  mutate(
    ecosystem = case_when(
      taxa %in% c("Coral", "Fish")  ~ "Marine",
      taxa %in% c("Plant", "Bird")  ~ "Terrestrial"
    )
  )

model_all_tax$ecosystem <- factor(
  model_all_tax$ecosystem,
  levels = c("Marine", "Terrestrial")
)


##

model_all_tax_scaled <- model_all_tax %>%
  
  # global variables
  mutate(
    age  = rescale_01(age),
    dist = rescale_01(dist)
  ) %>%
  
  # ecosystem variables
  group_by(ecosystem) %>%
  mutate(
    area = rescale_01(area),
    temp = rescale_01(temp)
  ) %>%
  ungroup()


### Scaling again

epsilon <- 1e-4  # um valor pequeno

model_all_tax_turn_pres_scaled <- model_all_tax_scaled %>%
  mutate(
    beta_sim_rescaled = beta_sim * (1 - 2 * epsilon) + epsilon
  )


# Model

mod_turn_pres_bayes_model2 <- brm(
  beta_sim_rescaled ~ area * taxa + dist * taxa + age * taxa + temp * taxa + (1 | group2),
  data = model_all_tax_turn_pres_scaled,
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


##
mod_turn_pres_bayes_model2 <- add_criterion(mod_turn_pres_bayes_model2, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  mod_turn_pres_bayes_model2,
  "/home/student/Islands_Biogeography/Inter_Data/model_beta_all_taxonomic_turn_present.rds"
)


# Assess model validity using LOO
loo(mod_turn_pres_bayes_model2)


# Plot

post_env <- mod_turn_pres_bayes_model2 %>%
  gather_draws(
    b_area,
    `b_area:taxaFish`, `b_area:taxaPlant`, `b_area:taxaBird`,
    b_dist,
    `b_taxaFish:dist`, `b_taxaPlant:dist`, `b_taxaBird:dist`,
    b_age,
    `b_taxaFish:age`,  `b_taxaPlant:age`,  `b_taxaBird:age`,
    b_temp,
    `b_taxaFish:temp`, `b_taxaPlant:temp`, `b_taxaBird:temp`
  )


#

post_env_eff <- post_env %>%
  mutate(
    Variable = case_when(
      grepl("area", .variable) ~ "Area",
      grepl("dist", .variable) ~ "Distance",
      grepl("age",  .variable) ~ "Age",
      grepl("temp", .variable) ~ "Temperature"
    )
  ) %>%
  group_by(.draw, Variable) %>%
  summarise(
    Coral = sum(.value[.variable %in% c("b_area", "b_dist", "b_age", "b_temp")]),
    Fish  = Coral + sum(.value[grepl("taxaFish",  .variable)]),
    Plant = Coral + sum(.value[grepl("taxaPlant", .variable)]),
    Bird  = Coral + sum(.value[grepl("taxaBird",  .variable)]),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(Coral, Fish, Plant, Bird),
    names_to = "Taxon",
    values_to = "Effect"
  )


#
tax_turn_pres <- ggplot(
  post_env_eff,
  aes(x = Effect, fill = Taxon)
) +
  stat_halfeye(
    adjust = 0.6,
    alpha = 0.6,
    .width = c(0.5, 0.8, 0.95),
    slab_color = NA
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_grid(Variable ~ .) +
  coord_cartesian(xlim = c(-20, 3)) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    strip.text.y = element_text(angle = 0, face = "bold"),
    legend.position = "top"
  ) +
  labs(
    x = "Posterior effect size on β-diversity turn",
    y = NULL,
    fill = "Taxonomic group"
  )




#### Beta taxonomic NEST PAST ----------------------------------------------------


# coral

coral_tax_past <- model_corals_taxonomic %>%
  mutate(
    Archipelago1 = Archipelago1.model,
    Archipelago2 = Archipelago2.model,
    taxa = "Coral",
    dist = Distance,      
    area = diff_past_reef_area,
    temp = diff_past_sst,
    age  = diff_age,
    group2 = Group2,
    group1 = Group1
  )


# fish

fish_tax_past <- model_fish_taxonomic %>%
  mutate(
    Archipelago1 = Archipelago1,
    Archipelago2 = Archipelago2,
    taxa = "Fish",
    dist = Distance,      
    area = diff_past_reef_area,
    temp = diff_past_sst,
    age  = diff_age,
    group2 = Group2,
    group1 = Group1
  )


# plant

plant_tax_past <- model_plants_taxonomic %>%
  mutate(
    Archipelago1 = Archipelago1,
    Archipelago2 = Archipelago2,
    taxa = "Plant",
    dist = Distance,      
    area = diff_past_area,
    temp = diff_past_air,
    age  = diff_age,
    group2 = Group2,
    group1 = Group1
  )


# bird

bird_tax_past <- model_birds_taxonomic %>%
  mutate(
    Archipelago1 = Archipelago1,
    Archipelago2 = Archipelago2,
    taxa = "Bird",
    dist = Distance,     
    area = diff_past_area,
    temp = diff_past_air,
    age  = diff_age,
    group2 = Group2,
    group1 = Group1
  )


##

coral_tax_past <- coral_tax_past %>%
  select(beta_sne, Archipelago1, Archipelago2, taxa, area, dist, age, temp, group2, group1)

fish_tax_past <- fish_tax_past %>%
  select(beta_sne, Archipelago1, Archipelago2, taxa, area, dist, age, temp, group2, group1)

plant_tax_past <- plant_tax_past %>%
  select(beta_sne, Archipelago1, Archipelago2, taxa, area, dist, age, temp, group2, group1)

bird_tax_past <- bird_tax_past %>%
  select(beta_sne, Archipelago1, Archipelago2, taxa, area, dist, age, temp, group2, group1)


## Merging

model_all_tax_past <- bind_rows(
  coral_tax_past,
  fish_tax_past,
  plant_tax_past,
  bird_tax_past
)

model_all_tax_past$taxa <- factor(
  model_all_tax_past$taxa,
  levels = c("Coral", "Fish", "Plant", "Bird")
)

# Scaling up variables per ecosystem type

# 0 and 1
rescale_01 <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}


model_all_tax <- model_all_tax_past %>%
  mutate(
    ecosystem = case_when(
      taxa %in% c("Coral", "Fish")  ~ "Marine",
      taxa %in% c("Plant", "Bird")  ~ "Terrestrial"
    )
  )

model_all_tax$ecosystem <- factor(
  model_all_tax$ecosystem,
  levels = c("Marine", "Terrestrial")
)


##

model_all_tax_scaled <- model_all_tax %>%
  
  # global variables
  mutate(
    age  = rescale_01(age),
    dist = rescale_01(dist)
  ) %>%
  
  # ecosystem variables
  group_by(ecosystem) %>%
  mutate(
    area = rescale_01(area),
    temp = rescale_01(temp)
  ) %>%
  ungroup()


### Scaling again

epsilon <- 1e-4  # um valor pequeno

model_all_tax_past_scaled <- model_all_tax_scaled %>%
  mutate(
    beta_sne_rescaled = beta_sne * (1 - 2 * epsilon) + epsilon
  )


## Model

mod_tax_bayes_past_model2 <- brm(
  beta_sne_rescaled ~ area * taxa + dist * taxa + age * taxa + temp * taxa + (1 | group2),
  data = model_all_tax_past_scaled,
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

 ##
mod_tax_bayes_past_model2 <- add_criterion(mod_tax_bayes_past_model2, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  mod_tax_bayes_past_model2,
  "/home/student/Islands_Biogeography/Inter_Data/model_beta_all_taxonomic_nest_past.rds"
)

# Assess model validity using LOO
loo(mod_tax_bayes_past_model2)


## Plot

post_env <- mod_tax_bayes_past_model2 %>%
  gather_draws(
    b_area,
    `b_area:taxaFish`, `b_area:taxaPlant`, `b_area:taxaBird`,
    b_dist,
    `b_taxaFish:dist`, `b_taxaPlant:dist`, `b_taxaBird:dist`,
    b_age,
    `b_taxaFish:age`,  `b_taxaPlant:age`,  `b_taxaBird:age`,
    b_temp,
    `b_taxaFish:temp`, `b_taxaPlant:temp`, `b_taxaBird:temp`
  )


#

post_env_eff <- post_env %>%
  mutate(
    Variable = case_when(
      grepl("area", .variable) ~ "Area",
      grepl("dist", .variable) ~ "Distance",
      grepl("age",  .variable) ~ "Age",
      grepl("temp", .variable) ~ "Temperature"
    )
  ) %>%
  group_by(.draw, Variable) %>%
  summarise(
    Coral = sum(.value[.variable %in% c("b_area", "b_dist", "b_age", "b_temp")]),
    Fish  = Coral + sum(.value[grepl("taxaFish",  .variable)]),
    Plant = Coral + sum(.value[grepl("taxaPlant", .variable)]),
    Bird  = Coral + sum(.value[grepl("taxaBird",  .variable)]),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(Coral, Fish, Plant, Bird),
    names_to = "Taxon",
    values_to = "Effect"
  )


#
tax_nest_past <- ggplot(
  post_env_eff,
  aes(x = Effect, fill = Taxon)
) +
  stat_halfeye(
    adjust = 0.6,
    alpha = 0.6,
    .width = c(0.5, 0.8, 0.95),
    slab_color = NA
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_grid(Variable ~ .) +
  coord_cartesian(xlim = c(-20, 3)) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    strip.text.y = element_text(angle = 0, face = "bold"),
    legend.position = "top"
  ) +
  labs(
    x = "Posterior effect size on β-diversity nest",
    y = NULL,
    fill = "Taxonomic group"
  )



  
#### Beta taxonomic TURN PAST ----------------------------------------------------


# coral

coral_tax_past <- model_corals_taxonomic %>%
  mutate(
    Archipelago1 = Archipelago1.model,
    Archipelago2 = Archipelago2.model,
    taxa = "Coral",
    dist = Distance,      
    area = diff_past_reef_area,
    temp = diff_past_sst,
    age  = diff_age,
    group2 = Group2,
    group1 = Group1
  )


# fish

fish_tax_past <- model_fish_taxonomic %>%
  mutate(
    Archipelago1 = Archipelago1,
    Archipelago2 = Archipelago2,
    taxa = "Fish",
    dist = Distance,      
    area = diff_past_reef_area,
    temp = diff_past_sst,
    age  = diff_age,
    group2 = Group2,
    group1 = Group1
  )


# plant

plant_tax_past <- model_plants_taxonomic %>%
  mutate(
    Archipelago1 = Archipelago1,
    Archipelago2 = Archipelago2,
    taxa = "Plant",
    dist = Distance,      
    area = diff_past_area,
    temp = diff_past_air,
    age  = diff_age,
    group2 = Group2,
    group1 = Group1
  )


# bird

bird_tax_past <- model_birds_taxonomic %>%
  mutate(
    Archipelago1 = Archipelago1,
    Archipelago2 = Archipelago2,
    taxa = "Bird",
    dist = Distance,     
    area = diff_past_area,
    temp = diff_past_air,
    age  = diff_age,
    group2 = Group2,
    group1 = Group1
  )



coral_tax_past <- coral_tax_past %>%
  select(beta_sim, Archipelago1, Archipelago2, taxa, area, dist, age, temp, group2, group1)

fish_tax_past <- fish_tax_past %>%
  select(beta_sim, Archipelago1, Archipelago2, taxa, area, dist, age, temp, group2, group1)

plant_tax_past <- plant_tax_past %>%
  select(beta_sim, Archipelago1, Archipelago2, taxa, area, dist, age, temp, group2, group1)

bird_tax_past <- bird_tax_past %>%
  select(beta_sim, Archipelago1, Archipelago2, taxa, area, dist, age, temp, group2, group1)

# Merging

model_all_tax_past_sim <- bind_rows(
  coral_tax_past,
  fish_tax_past,
  plant_tax_past,
  bird_tax_past
)

model_all_tax_past_sim$taxa <- factor(
  model_all_tax_past_sim$taxa,
  levels = c("Coral", "Fish", "Plant", "Bird")
)


# Scaling up variables per ecosystem type

# 0 and 1
rescale_01 <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}


model_all_tax <- model_all_tax_past_sim %>%
  mutate(
    ecosystem = case_when(
      taxa %in% c("Coral", "Fish")  ~ "Marine",
      taxa %in% c("Plant", "Bird")  ~ "Terrestrial"
    )
  )

model_all_tax$ecosystem <- factor(
  model_all_tax$ecosystem,
  levels = c("Marine", "Terrestrial")
)


##

model_all_tax_scaled <- model_all_tax %>%
  
  # global variables
  mutate(
    age  = rescale_01(age),
    dist = rescale_01(dist)
  ) %>%
  
  # ecosystem variables
  group_by(ecosystem) %>%
  mutate(
    area = rescale_01(area),
    temp = rescale_01(temp)
  ) %>%
  ungroup()


### Scaling again

epsilon <- 1e-4  # um valor pequeno

model_all_tax_past_sim_scaled <- model_all_tax_scaled %>%
  mutate(
    beta_sim_rescaled = beta_sim * (1 - 2 * epsilon) + epsilon
  )


# Model

mod_tax_bayes_past_sim_model2 <- brm(
  beta_sim_rescaled ~ area * taxa + dist * taxa + age * taxa + temp * taxa + (1 | group2),
  data = model_all_tax_past_sim_scaled,
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

##
mod_tax_bayes_past_sim_model2 <- add_criterion(mod_tax_bayes_past_sim_model2, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  mod_tax_bayes_past_sim_model2,
  "/home/student/Islands_Biogeography/Inter_Data/model_beta_all_taxonomic_turn_past.rds"
)

# Assess model validity using LOO
loo(mod_tax_bayes_past_sim_model2)


#Plot

post_env <- mod_tax_bayes_past_sim_model2 %>%
  gather_draws(
    b_area,
    `b_area:taxaFish`, `b_area:taxaPlant`, `b_area:taxaBird`,
    b_dist,
    `b_taxaFish:dist`, `b_taxaPlant:dist`, `b_taxaBird:dist`,
    b_age,
    `b_taxaFish:age`,  `b_taxaPlant:age`,  `b_taxaBird:age`,
    b_temp,
    `b_taxaFish:temp`, `b_taxaPlant:temp`, `b_taxaBird:temp`
  )


#

post_env_eff <- post_env %>%
  mutate(
    Variable = case_when(
      grepl("area", .variable) ~ "Area",
      grepl("dist", .variable) ~ "Distance",
      grepl("age",  .variable) ~ "Age",
      grepl("temp", .variable) ~ "Temperature"
    )
  ) %>%
  group_by(.draw, Variable) %>%
  summarise(
    Coral = sum(.value[.variable %in% c("b_area", "b_dist", "b_age", "b_temp")]),
    Fish  = Coral + sum(.value[grepl("taxaFish",  .variable)]),
    Plant = Coral + sum(.value[grepl("taxaPlant", .variable)]),
    Bird  = Coral + sum(.value[grepl("taxaBird",  .variable)]),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(Coral, Fish, Plant, Bird),
    names_to = "Taxon",
    values_to = "Effect"
  )


#
tax_turn_past <- ggplot(
  post_env_eff,
  aes(x = Effect, fill = Taxon)
) +
  stat_halfeye(
    adjust = 0.6,
    alpha = 0.6,
    .width = c(0.5, 0.8, 0.95),
    slab_color = NA
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_grid(Variable ~ .) +
  coord_cartesian(xlim = c(-20, 3)) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    strip.text.y = element_text(angle = 0, face = "bold"),
    legend.position = "top"
  ) +
  labs(
    x = "Posterior effect size on β-diversity turn",
    y = NULL,
    fill = "Taxonomic group"
  )




################################################################################
################################################################################

## Ploting all past and present taxonomic indeces

# Nest

plot_tax_nest <- tax_nest_pres +
  labs(title = "Present") +
  tax_nest_past +
  labs(title = "Past") +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "Taxonomic beta diversity – Nestedness (βSNE)",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  )

# Turn

plot_tax_turn <- tax_turn_pres +
  labs(title = "Present") +
  tax_turn_past +
  labs(title = "Past") +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "Taxonomic beta diversity – Turnover (βSIM)",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  )


## Taxa overlapping

# plot_beta_with_taxa <- function(model, xlab = "Posterior effect size") {
#   
#   ## Taxonomic effects (sobrepostos)
#   post_tax <- model %>%
#     gather_draws(
#       b_Intercept,
#       b_taxaBird,
#       b_taxaFish,
#       b_taxaPlant
#     ) %>%
#     mutate(
#       parameter = case_when(
#         .variable == "b_Intercept" ~ "Intercept",
#         .variable == "b_taxaBird"  ~ "Birds",
#         .variable == "b_taxaFish"  ~ "Fish",
#         .variable == "b_taxaPlant" ~ "Plants"
#       ),
#       y_plot = "Taxonomic groups"
#     )
#   
#   ## Environmental effects
#   post_env <- model %>%
#     gather_draws(
#       b_area,
#       b_dist,
#       b_age,
#       b_temp
#     ) %>%
#     mutate(
#       y_plot = case_when(
#         .variable == "b_area" ~ "Area",
#         .variable == "b_dist" ~ "Distance",
#         .variable == "b_age"  ~ "Age",
#         .variable == "b_temp" ~ "Temperature"
#       )
#     )
#   
#   ## Combine
#   post_all <- bind_rows(post_tax, post_env) %>%
#     mutate(
#       y_plot = factor(
#         y_plot,
#         levels = c(
#           "Taxonomic groups",
#           "Area",
#           "Distance",
#           "Age",
#           "Temperature"
#         )
#       )
#     )
#   
#   ## Plot
#   ggplot(post_all, aes(x = .value, y = y_plot)) +
#     
#     ## Taxa (sobrepostos)
#     stat_halfeye(
#       data = subset(post_all, y_plot == "Taxonomic groups"),
#       aes(fill = parameter),
#       adjust = 0.5,
#       .width = c(0.5, 0.8, 0.95),
#       alpha = 0.7
#     ) +
#     
#     ## Ambiente
#     stat_halfeye(
#       data = subset(post_all, y_plot != "Taxonomic groups"),
#       fill = "grey70",
#       adjust = 0.5,
#       .width = c(0.5, 0.8, 0.95),
#       alpha = 0.7
#     ) +
#     
#     geom_vline(xintercept = 0, linetype = "dashed") +
#     scale_y_discrete(limits = rev) +
#     theme_minimal() +
#     theme(
#       legend.position = "top",
#       panel.grid.major.y = element_blank()
#     ) +
#     labs(
#       x = xlab,
#       y = NULL,
#       fill = "Taxa"
#     )
# }
# 
# ##
# tax_nest_pres <- plot_beta_with_taxa(
#   mod_tax_bayes,
#   xlab = "Posterior effect size (taxonomic βSNE)"
# )
# 
# ##
# tax_nest_past <- plot_beta_with_taxa(
#   mod_tax_bayes_past,
#   xlab = "Posterior effect size (taxonomic βSNE)"
# )
# 
# ##
# tax_turn_pres <- plot_beta_with_taxa(
#   mod_turn_pres_bayes,
#   xlab = "Posterior effect size (taxonomic βSIM)"
# )
# 
# ##
# tax_turn_past <- plot_beta_with_taxa(
#   mod_tax_bayes_past_sim,
#   xlab = "Posterior effect size (taxonomic βSIM)"
# )
# 
# 
# ## final plot
# 
# # Nest
# plot_tax_nest <- tax_nest_past +
#   labs(title = "Past") +
#   tax_nest_pres +
#   labs(title = "Present") +
#   plot_layout(ncol = 2) +
#   plot_annotation(
#     title = "Taxonomic beta diversity – Nestedness (βSNE)",
#     theme = theme(
#       plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
#     )
#   )
# 
# # Turn
# plot_tax_nest <- tax_turn_past +
#   labs(title = "Past") +
#   tax_turn_pres +
#   labs(title = "Present") +
#   plot_layout(ncol = 2) +
#   plot_annotation(
#     title = "Taxonomic beta diversity – Turnover (βSIM)",
#     theme = theme(
#       plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
#     )
#   )
# 
# 
# 






#### Beta functional NEST PRESENT ----------------------------------------------------

# Coral
coral_fun_nest <- model_corals_functional %>%
  mutate(
    Archipelago1 = Archipelago1.model,
    Archipelago2 = Archipelago2.model,
    taxa = "Coral",
    dist = Distance,
    area = diff_reef_area,
    temp = diff_sst,
    age  = diff_age,
    group2 = Group2,
    group1 = Group1
  ) %>%
  select(beta_sne, Archipelago1, Archipelago2, taxa, area, dist, age, temp, group2, group1)

# Fish
fish_fun_nest <- model_fish_functional %>%
  mutate(
    Archipelago1 = Archipelago1,
    Archipelago2 = Archipelago2,
    taxa = "Fish",
    dist = Distance,
    area = diff_reef_area,
    temp = diff_sst,
    age  = diff_age,
    group2 = Group2,
    group1 = Group1
  ) %>%
  select(beta_sne, Archipelago1, Archipelago2, taxa, area, dist, age, temp, group2, group1)

# Plant
plant_fun_nest <- model_plants_functional %>%
  mutate(
    Archipelago1 = Archipelago1,
    Archipelago2 = Archipelago2,
    taxa = "Plant",
    dist = Distance,
    area = diff_present_area,
    temp = diff_present_air,
    age  = diff_age,
    group2 = Group2,
    group1 = Group1
  ) %>%
  select(beta_sne, Archipelago1, Archipelago2, taxa, area, dist, age, temp, group2, group1)

# Bird
bird_fun_nest <- model_birds_functional %>%
  mutate(
    Archipelago1 = Archipelago1,
    Archipelago2 = Archipelago2,
    taxa = "Bird",
    dist = Distance,
    area = diff_present_area,
    temp = diff_present_air,
    age  = diff_age,
    group2 = Group2,
    group1 = Group1
  ) %>%
  select(beta_sne, Archipelago1, Archipelago2, taxa, area, dist, age, temp, group2, group1)

# Merge all
model_all_fun_nest <- bind_rows(coral_fun_nest, fish_fun_nest, plant_fun_nest, bird_fun_nest)
model_all_fun_nest$taxa <- factor(model_all_fun_nest$taxa, levels = c("Coral", "Fish", "Plant", "Bird"))

# Scaling up variables per ecosystem type

# 0 and 1
rescale_01 <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}


model_all_tax <- model_all_fun_nest %>%
  mutate(
    ecosystem = case_when(
      taxa %in% c("Coral", "Fish")  ~ "Marine",
      taxa %in% c("Plant", "Bird")  ~ "Terrestrial"
    )
  )

model_all_tax$ecosystem <- factor(
  model_all_tax$ecosystem,
  levels = c("Marine", "Terrestrial")
)


##

model_all_tax_scaled <- model_all_tax %>%
  
  # global variables
  mutate(
    age  = rescale_01(age),
    dist = rescale_01(dist)
  ) %>%
  
  # ecosystem variables
  group_by(ecosystem) %>%
  mutate(
    area = rescale_01(area),
    temp = rescale_01(temp)
  ) %>%
  ungroup()


### Scaling again

epsilon <- 1e-4  # um valor pequeno

model_all_fun_nest_scaled <- model_all_tax_scaled %>%
  mutate(
    beta_sne_rescaled = beta_sne * (1 - 2 * epsilon) + epsilon
  )


# Model
mod_fun_bayes_nest_model2 <- brm(
  beta_sne_rescaled ~ area * taxa + dist * taxa + age * taxa + temp * taxa + (1 | group2),
  data = model_all_fun_nest_scaled,
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


##
mod_fun_bayes_nest_model2 <- add_criterion(mod_fun_bayes_nest_model2, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  mod_fun_bayes_nest_model2,
  "/home/student/Islands_Biogeography/Inter_Data/model_beta_all_functional_nest_present.rds"
)


# Assess model validity using LOO
loo(mod_fun_bayes_nest_model2)


# Plot
post_env_fun_nest <- mod_fun_bayes_nest_model2 %>%
  gather_draws(b_area, b_dist, b_age, b_temp) %>%
  mutate(parameter = case_when(
    .variable == "b_area" ~ "Area",
    .variable == "b_dist" ~ "Distance",
    .variable == "b_age"  ~ "Age",
    .variable == "b_temp" ~ "Temperature"
  ))

fun_nest_pres <- ggplot(post_env_fun_nest, aes(x = .value, y = parameter, fill = parameter)) +
  stat_halfeye(adjust = 0.5, .width = c(0.5, 0.8, 0.95), alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  theme(legend.position = "none", panel.grid.major.y = element_blank()) +
  labs(x = "Posterior effect size", y = NULL)



#### Beta functional TURN PRESENT ----------------------------------------------------

coral_fun_turn <- model_corals_functional %>%
  mutate(taxa = "Coral", dist = Distance, area = diff_reef_area, temp = diff_sst,
         age = diff_age, group2 = Group2, group1 = Group1) %>%
  select(beta_sim, Archipelago1.model, Archipelago2.model, taxa, area, dist, age, temp, group2, group1)

fish_fun_turn <- model_fish_functional %>%
  mutate(taxa = "Fish", dist = Distance, area = diff_reef_area, temp = diff_sst,
         age = diff_age, group2 = Group2, group1 = Group1) %>%
  select(beta_sim, Archipelago1, Archipelago2, taxa, area, dist, age, temp, group2, group1)

plant_fun_turn <- model_plants_functional %>%
  mutate(taxa = "Plant", dist = Distance, area = diff_present_area, temp = diff_present_air,
         age = diff_age, group2 = Group2, group1 = Group1) %>%
  select(beta_sim, Archipelago1, Archipelago2, taxa, area, dist, age, temp, group2, group1)

bird_fun_turn <- model_birds_functional %>%
  mutate(taxa = "Bird", dist = Distance, area = diff_present_area, temp = diff_present_air,
         age = diff_age, group2 = Group2, group1 = Group1) %>%
  select(beta_sim, Archipelago1, Archipelago2, taxa, area, dist, age, temp, group2, group1)

##
model_all_fun_turn <- bind_rows(coral_fun_turn, fish_fun_turn, plant_fun_turn, bird_fun_turn)
model_all_fun_turn$taxa <- factor(model_all_fun_turn$taxa, levels = c("Coral", "Fish", "Plant", "Bird"))

# Scaling up variables per ecosystem type

# 0 and 1
rescale_01 <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}


model_all_tax <- model_all_fun_turn %>%
  mutate(
    ecosystem = case_when(
      taxa %in% c("Coral", "Fish")  ~ "Marine",
      taxa %in% c("Plant", "Bird")  ~ "Terrestrial"
    )
  )

model_all_tax$ecosystem <- factor(
  model_all_tax$ecosystem,
  levels = c("Marine", "Terrestrial")
)


##

model_all_tax_scaled <- model_all_tax %>%
  
  # global variables
  mutate(
    age  = rescale_01(age),
    dist = rescale_01(dist)
  ) %>%
  
  # ecosystem variables
  group_by(ecosystem) %>%
  mutate(
    area = rescale_01(area),
    temp = rescale_01(temp)
  ) %>%
  ungroup()


### Scaling again

epsilon <- 1e-4  # um valor pequeno

model_all_fun_turn_scaled <- model_all_tax_scaled %>%
  mutate(
    beta_sim_rescaled = beta_sim * (1 - 2 * epsilon) + epsilon
  )

# Model
mod_fun_bayes_turn_model2 <- brm(
  beta_sim_rescaled ~ area * taxa + dist * taxa + age * taxa + temp * taxa + (1 | group2),
  data = model_all_fun_turn_scaled,
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

##
mod_fun_bayes_turn_model2 <- add_criterion(mod_fun_bayes_turn_model2, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  mod_fun_bayes_turn_model2,
  "/home/student/Islands_Biogeography/Inter_Data/model_beta_all_functional_turn_present.rds"
)


# Assess model validity using LOO
loo(mod_fun_bayes_turn_model2)


# Plot
post_env_fun_turn <- mod_fun_bayes_turn_model2 %>%
  gather_draws(b_area, b_dist, b_age, b_temp) %>%
  mutate(parameter = case_when(
    .variable == "b_area" ~ "Area",
    .variable == "b_dist" ~ "Distance",
    .variable == "b_age"  ~ "Age",
    .variable == "b_temp" ~ "Temperature"
  ))

#
fun_turn_pres <- ggplot(post_env_fun_turn, aes(x = .value, y = parameter, fill = parameter)) +
  stat_halfeye(adjust = 0.5, .width = c(0.5, 0.8, 0.95), alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  theme(legend.position = "none", panel.grid.major.y = element_blank()) +
  labs(x = "Posterior effect size (functional βSIM)", y = NULL)




#### Beta functional NEST PAST ----------------------------------------------------

# Coral
coral_fun_nest_past <- model_corals_functional %>%
  mutate(
    Archipelago1 = Archipelago1.model,
    Archipelago2 = Archipelago2.model,
    taxa = "Coral",
    dist = Distance,
    area = diff_past_reef_area,
    temp = diff_past_sst,
    age  = diff_age,
    group2 = Group2,
    group1 = Group1
  ) %>%
  select(beta_sne, Archipelago1, Archipelago2, taxa, area, dist, age, temp, group2, group1)

# Fish
fish_fun_nest_past <- model_fish_functional %>%
  mutate(
    Archipelago1 = Archipelago1,
    Archipelago2 = Archipelago2,
    taxa = "Fish",
    dist = Distance,
    area = diff_past_reef_area,
    temp = diff_past_sst,
    age  = diff_age,
    group2 = Group2,
    group1 = Group1
  ) %>%
  select(beta_sne, Archipelago1, Archipelago2, taxa, area, dist, age, temp, group2, group1)

# Plants
plant_fun_nest_past <- model_plants_functional %>%
  mutate(
    Archipelago1 = Archipelago1,
    Archipelago2 = Archipelago2,
    taxa = "Plant",
    dist = Distance,
    area = diff_past_area,
    temp = diff_past_air,
    age  = diff_age,
    group2 = Group2,
    group1 = Group1
  ) %>%
  select(beta_sne, Archipelago1, Archipelago2, taxa, area, dist, age, temp, group2, group1)

# Birds
bird_fun_nest_past <- model_birds_functional %>%
  mutate(
    Archipelago1 = Archipelago1,
    Archipelago2 = Archipelago2,
    taxa = "Bird",
    dist = Distance,
    area = diff_past_area,
    temp = diff_past_air,
    age  = diff_age,
    group2 = Group2,
    group1 = Group1
  ) %>%
  select(beta_sne, Archipelago1, Archipelago2, taxa, area, dist, age, temp, group2, group1)

# Merge all
model_all_fun_nest_past <- bind_rows(coral_fun_nest_past, fish_fun_nest_past, plant_fun_nest_past, bird_fun_nest_past)
model_all_fun_nest_past$taxa <- factor(model_all_fun_nest_past$taxa, levels = c("Coral", "Fish", "Plant", "Bird"))

# Scaling up variables per ecosystem type

# 0 and 1
rescale_01 <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}


model_all_tax <- model_all_fun_nest_past %>%
  mutate(
    ecosystem = case_when(
      taxa %in% c("Coral", "Fish")  ~ "Marine",
      taxa %in% c("Plant", "Bird")  ~ "Terrestrial"
    )
  )

model_all_tax$ecosystem <- factor(
  model_all_tax$ecosystem,
  levels = c("Marine", "Terrestrial")
)


##

model_all_tax_scaled <- model_all_tax %>%
  
  # global variables
  mutate(
    age  = rescale_01(age),
    dist = rescale_01(dist)
  ) %>%
  
  # ecosystem variables
  group_by(ecosystem) %>%
  mutate(
    area = rescale_01(area),
    temp = rescale_01(temp)
  ) %>%
  ungroup()


### Scaling again

epsilon <- 1e-4  

model_all_fun_nest_scaled_past <- model_all_tax_scaled %>%
  mutate(
    beta_sne_rescaled_past = beta_sne * (1 - 2 * epsilon) + epsilon
  )


# Model
mod_fun_bayes_nest_past_model2 <- brm(
  beta_sne_rescaled_past ~ area * taxa + dist * taxa + age * taxa + temp * taxa + (1 | group2),
  data = model_all_fun_nest_scaled_past,
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


##
mod_fun_bayes_nest_past_model2 <- add_criterion(mod_fun_bayes_nest_past_model2, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  mod_fun_bayes_nest_past_model2,
  "/home/student/Islands_Biogeography/Inter_Data/model_beta_all_functional_nest_past.rds"
)


# Assess model validity using LOO
loo(mod_fun_bayes_nest_past_model2)


# Plot
post_env_fun_nest_past <- mod_fun_bayes_nest_past_model2 %>%
  gather_draws(b_area, b_dist, b_age, b_temp) %>%
  mutate(parameter = case_when(
    .variable == "b_area" ~ "Area",
    .variable == "b_dist" ~ "Distance",
    .variable == "b_age"  ~ "Age",
    .variable == "b_temp" ~ "Temperature"
  ))

fun_nest_past <- ggplot(post_env_fun_nest_past, aes(x = .value, y = parameter, fill = parameter)) +
  stat_halfeye(adjust = 0.5, .width = c(0.5, 0.8, 0.95), alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  theme(legend.position = "none", panel.grid.major.y = element_blank()) +
  labs(x = "Posterior effect size", y = NULL)



#### Beta functional TURN PAST ----------------------------------------------------

coral_fun_turn_past <- model_corals_functional %>%
  mutate(taxa = "Coral", dist = Distance, area = diff_past_reef_area, temp = diff_past_sst,
         age = diff_age, group2 = Group2, group1 = Group1) %>%
  select(beta_sim, Archipelago1.model, Archipelago2.model, taxa, area, dist, age, temp, group2, group1)

fish_fun_turn_past <- model_fish_functional %>%
  mutate(taxa = "Fish", dist = Distance, area = diff_past_reef_area, temp = diff_past_sst,
         age = diff_age, group2 = Group2, group1 = Group1) %>%
  select(beta_sim, Archipelago1, Archipelago2, taxa, area, dist, age, temp, group2, group1)

plant_fun_turn_past <- model_plants_functional %>%
  mutate(taxa = "Plant", dist = Distance, area = diff_past_area, temp = diff_past_air,
         age = diff_age, group2 = Group2, group1 = Group1) %>%
  select(beta_sim, Archipelago1, Archipelago2, taxa, area, dist, age, temp, group2, group1)

bird_fun_turn_past <- model_birds_functional %>%
  mutate(taxa = "Bird", dist = Distance, area = diff_past_area, temp = diff_past_air,
         age = diff_age, group2 = Group2, group1 = Group1) %>%
  select(beta_sim, Archipelago1, Archipelago2, taxa, area, dist, age, temp, group2, group1)

##
model_all_fun_turn_past <- bind_rows(coral_fun_turn_past, fish_fun_turn_past, plant_fun_turn_past, bird_fun_turn_past)
model_all_fun_turn_past$taxa <- factor(model_all_fun_turn_past$taxa, levels = c("Coral", "Fish", "Plant", "Bird"))

# Scaling up variables per ecosystem type

# 0 and 1
rescale_01 <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}


model_all_tax <- model_all_fun_turn_past %>%
  mutate(
    ecosystem = case_when(
      taxa %in% c("Coral", "Fish")  ~ "Marine",
      taxa %in% c("Plant", "Bird")  ~ "Terrestrial"
    )
  )

model_all_tax$ecosystem <- factor(
  model_all_tax$ecosystem,
  levels = c("Marine", "Terrestrial")
)


##

model_all_tax_scaled <- model_all_tax %>%
  
  # global variables
  mutate(
    age  = rescale_01(age),
    dist = rescale_01(dist)
  ) %>%
  
  # ecosystem variables
  group_by(ecosystem) %>%
  mutate(
    area = rescale_01(area),
    temp = rescale_01(temp)
  ) %>%
  ungroup()


### Scaling again

epsilon <- 1e-4  # um valor pequeno

model_all_fun_turn_past_scaled <- model_all_tax_scaled %>%
  mutate(
    beta_sim_rescaled_past = beta_sim * (1 - 2 * epsilon) + epsilon
  )


# Model
mod_fun_bayes_turn_past_model2 <- brm(
  beta_sim_rescaled_past ~ area * taxa + dist * taxa + age * taxa + temp * taxa + (1 | group2),
  data = model_all_fun_turn_past_scaled,
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


##
mod_fun_bayes_turn_past_model2 <- add_criterion(mod_fun_bayes_turn_past_model2, "loo", save_psis = TRUE, reloo = T)

saveRDS(
  mod_fun_bayes_turn_past_model2,
  "/home/student/Islands_Biogeography/Inter_Data/model_beta_all_functional_turn_past.rds"
)


# Assess model validity using LOO
loo(mod_fun_bayes_turn_past_model2)


# Plot
post_env_fun_turn_past <- mod_fun_bayes_turn_past_model2 %>%
  gather_draws(b_area, b_dist, b_age, b_temp) %>%
  mutate(parameter = case_when(
    .variable == "b_area" ~ "Area",
    .variable == "b_dist" ~ "Distance",
    .variable == "b_age"  ~ "Age",
    .variable == "b_temp" ~ "Temperature"
  ))

#
fun_turn_past <- ggplot(post_env_fun_turn_past, aes(x = .value, y = parameter, fill = parameter)) +
  stat_halfeye(adjust = 0.5, .width = c(0.5, 0.8, 0.95), alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  theme(legend.position = "none", panel.grid.major.y = element_blank()) +
  labs(x = "Posterior effect size (functional βSIM)", y = NULL)


################################################################################
################################################################################

## Ploting all past and present functional indeces

# Nest

plot_fun_nest_past <- fun_nest_pres +
  labs(title = "Present") +
  tax_nest_past +
  labs(title = "Past") +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "Functional beta diversity – Nestedness (βSNE)",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  )

# Turn

plot_fun_turn <- fun_turn_pres +
  labs(title = "Present") +
  fun_turn_past +
  labs(title = "Past") +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "Taxonomic beta diversity – Turnover (βSIM)",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  )


## Taxa overlapping
# 
# plot_beta_functional <- function(model, xlab = "Posterior effect size") {
#   
#   ## Taxonomic effects (sobrepostos)
#   post_tax <- model %>%
#     gather_draws(
#       b_Intercept,
#       b_taxaBird,
#       b_taxaFish,
#       b_taxaPlant
#     ) %>%
#     mutate(
#       parameter = case_when(
#         .variable == "b_Intercept" ~ "Intercept",
#         .variable == "b_taxaBird"  ~ "Birds",
#         .variable == "b_taxaFish"  ~ "Fish",
#         .variable == "b_taxaPlant" ~ "Plants"
#       ),
#       y_plot = "Taxonomic groups"
#     )
#   
#   ## Environmental effects
#   post_env <- model %>%
#     gather_draws(
#       b_area,
#       b_dist,
#       b_age,
#       b_temp
#     ) %>%
#     mutate(
#       y_plot = case_when(
#         .variable == "b_area" ~ "Area",
#         .variable == "b_dist" ~ "Distance",
#         .variable == "b_age"  ~ "Age",
#         .variable == "b_temp" ~ "Temperature"
#       )
#     )
#   
#   ## Combine
#   post_all <- bind_rows(post_tax, post_env) %>%
#     mutate(
#       y_plot = factor(
#         y_plot,
#         levels = c(
#           "Taxonomic groups",
#           "Area",
#           "Distance",
#           "Age",
#           "Temperature"
#         )
#       )
#     )
#   
#   ## Plot
#   ggplot(post_all, aes(x = .value, y = y_plot)) +
#     
#     ## Taxa (sobrepostos)
#     stat_halfeye(
#       data = subset(post_all, y_plot == "Taxonomic groups"),
#       aes(fill = parameter),
#       adjust = 0.5,
#       .width = c(0.5, 0.8, 0.95),
#       alpha = 0.7
#     ) +
#     
#     ## Ambiente
#     stat_halfeye(
#       data = subset(post_all, y_plot != "Taxonomic groups"),
#       fill = "grey70",
#       adjust = 0.5,
#       .width = c(0.5, 0.8, 0.95),
#       alpha = 0.7
#     ) +
#     
#     geom_vline(xintercept = 0, linetype = "dashed") +
#     scale_y_discrete(limits = rev) +
#     theme_minimal() +
#     theme(
#       legend.position = "top",
#       panel.grid.major.y = element_blank()
#     ) +
#     labs(
#       x = xlab,
#       y = NULL,
#       fill = "Taxa"
#     )
# }
# 
# ##
# fun_nest_pres <- plot_beta_functional(
#   mod_fun_bayes_nest,
#   xlab = "Posterior effect size (functional βSNE)"
# )
# 
# ##
# fun_nest_past <- plot_beta_functional(
#   mod_fun_bayes_nest_past,
#   xlab = "Posterior effect size (functional βSNE)"
# )
# 
# ##
# fun_turn_pres <- plot_beta_functional(
#   mod_fun_bayes_turn,
#   xlab = "Posterior effect size (functional βSIM)"
# )
# 
# 
# ##
# fun_turn_past <- plot_beta_functional(
#   mod_fun_bayes_turn_past,
#   xlab = "Posterior effect size (functional βSIM)"
# )
# 
# 
# # Plot
# 
# #Nest
# plot_fun_nest <- fun_nest_past +
#   labs(title = "Past") +
#   fun_nest_pres +
#   labs(title = "Present") +
#   plot_layout(ncol = 2) +
#   plot_annotation(
#     title = "Functional beta diversity – Nestedness (βSNE)",
#     theme = theme(
#       plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
#     )
#   )
# 
# #Turn
# plot_fun_nest <- fun_turn_past +
#   labs(title = "Past") +
#   fun_turn_pres +
#   labs(title = "Present") +
#   plot_layout(ncol = 2) +
#   plot_annotation(
#     title = "Functional beta diversity – Turnover (βSIM)",
#     theme = theme(
#       plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
#     )
#   )


##############################################################################
##############################################################################

# Final plot - Taxonomic + Functional = Past + Present

extract_effects <- function(model, response_label, time_label, comp_label) {
  
  # Extrair draws como dataframe
  draws <- as_draws_df(model)
  
  # Manter apenas coeficientes fixos (b_)
  draws_long <- draws %>%
    select(starts_with("b_")) %>%
    mutate(.draw = row_number()) %>%
    pivot_longer(
      cols = - .draw,
      names_to = ".variable",
      values_to = ".value"
    )
  
  # Identificar variável ambiental
  draws_long <- draws_long %>%
    mutate(
      Variable = case_when(
        str_detect(.variable, "area") ~ "Area",
        str_detect(.variable, "dist") ~ "Distance",
        str_detect(.variable, "age")  ~ "Age",
        str_detect(.variable, "temp") ~ "Temperature",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(Variable))
  
  # efeitos principais (Coral = baseline)
  base_effects <- draws_long %>%
    filter(!str_detect(.variable, "taxa")) %>%
    select(.draw, Variable, base = .value)
  
  # interações
  interaction_effects <- draws_long %>%
    filter(str_detect(.variable, "taxa")) %>%
    mutate(
      Taxon = case_when(
        str_detect(.variable, "Fish")  ~ "Fish",
        str_detect(.variable, "Plant") ~ "Plant",
        str_detect(.variable, "Bird")  ~ "Bird"
      )
    ) %>%
    select(.draw, Variable, Taxon, interaction = .value)
  
  # Coral = apenas efeito base
  coral <- base_effects %>%
    mutate(
      Taxon = "Coral",
      Effect = base
    ) %>%
    select(.draw, Variable, Taxon, Effect)
  
  # Outros taxa = base + interação
  others <- interaction_effects %>%
    left_join(base_effects, by = c(".draw", "Variable")) %>%
    mutate(
      Effect = base + interaction
    ) %>%
    select(.draw, Variable, Taxon, Effect)
  
  bind_rows(coral, others) %>%
    mutate(
      Response = response_label,
      Time = time_label,
      Component = comp_label
    )
}


###

tax_nest_present  <- extract_effects(mod_tax_bayes_model2, "Taxonomic", "Present", "Nest")
tax_turn_present  <- extract_effects(mod_turn_pres_bayes_model2, "Taxonomic", "Present", "Turn")
tax_nest_past     <- extract_effects(mod_tax_bayes_past_model2, "Taxonomic", "Past", "Nest")
tax_turn_past     <- extract_effects(mod_tax_bayes_past_sim_model2, "Taxonomic", "Past", "Turn")

all_taxonomic <- bind_rows(
  tax_nest_present,
  tax_turn_present,
  tax_nest_past,
  tax_turn_past
)


## Plot

ggplot(all_taxonomic,
       aes(x = Effect, fill = Taxon)) +
  stat_halfeye(
    adjust = 0.6,
    alpha = 0.6,
    .width = c(0.5, 0.8, 0.95),
    slab_color = NA
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  coord_cartesian(xlim = c(-2.5, 3.6)) +
  facet_grid(Variable ~ Time + Component, scales = "free_x") +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  ) +
  labs(
    x = "Posterior effect size",
    y = NULL,
    fill = "Taxonomic group"
  )

####