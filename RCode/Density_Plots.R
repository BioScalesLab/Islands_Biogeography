

######
####
##
# ---- Script for plotting Bayesian models results as density plots

## Packages

library(brms)
library(bayesplot)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidybayes)


### EX. CORALS TAXONOMIC PRESENT - TOTAL BETA 

# extract all posterior draws
posterior <- as_draws_df(model_beta_corals_taxonomic_present_bayes_all)

# check parameter names
colnames(posterior)


##

# select only parameters of interest (B_ fixed)
posterior_params <- posterior %>%
  select(starts_with("b_")) %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value")


##

ggplot(posterior_params, aes(x = value, y = ..density..)) +
  geom_histogram(aes(y=..density..), bins = 50, fill = "gray", color = "black", alpha = 0.7) +
  facet_wrap(~parameter, scales = "free") +
  geom_vline(xintercept = 0, linetype="dashed", color="red") +
  theme_minimal() +
  labs(x = "Parameter value", y = "Posterior density")


##

# Add posterior probabilities or quantiles


posterior %>%
  gather_draws(b_diff_isolation_sc, b_diff_reef_area_sc, b_diff_age_sc, b_dist_sc, b_diff_present_sst) %>%
  ggplot(aes(x = .value)) +
  stat_halfeye(adjust = 0.5, fill = "gray", .width = c(0.5, 0.8, 0.95)) +
  facet_wrap(~.variable, scales = "free") +
  geom_vline(xintercept = 0, linetype="dashed", color="red") +
  theme_minimal()





#######################################

### EX. CORALS TAXONOMIC PAST - TOTAL BETA 


# The coefficient names are different:
#   
# Present: b_diff_isolation_sc
# 
# Past: b_diff_past_isolation_sc
# 
# To compare them, we need to standardize the names.
# 
# Extract and rename – PRESENT model


post_present <- model_beta_corals_taxonomic_present_bayes_all %>%
  gather_draws(
    b_diff_isolation_sc,
    b_diff_reef_area_sc,
    b_diff_age_sc,
    b_dist_sc
  ) %>%
  mutate(
    time = "Present",
    parameter = case_when(
      .variable == "b_diff_isolation_sc" ~ "Isolation",
      .variable == "b_diff_reef_area_sc" ~ "Reef area",
      .variable == "b_diff_age_sc"       ~ "Age",
      .variable == "b_dist_sc"           ~ "Distance"
    )
  )

# # Extract and rename – PAST model
post_past <- model_beta_corals_taxonomic_past_bayes_all %>%
  gather_draws(
    b_diff_past_isolation_sc,
    b_diff_past_reef_area_sc,
    b_diff_age_sc,
    b_dist_sc
  ) %>%
  mutate(
    time = "Past",
    parameter = case_when(
      .variable == "b_diff_past_isolation_sc" ~ "Isolation",
      .variable == "b_diff_past_reef_area_sc" ~ "Reef area",
      .variable == "b_diff_age_sc"            ~ "Age",
      .variable == "b_dist_sc"                ~ "Distance"
    )
  )

##

##

# This object contains:
#   
#   .value → posterior draws
# 
# parameter → Isolation, Reef area, Age, Distance
# 
# time → Present / Past


ggplot(posterior_compare,
       aes(x = .value, y = time, fill = time)) +
  stat_halfeye(
    adjust = 0.5,
    .width = c(0.5, 0.8, 0.95),
    alpha = 0.7
  ) +
  facet_wrap(~parameter, scales = "free_x") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Posterior effect size",
    y = NULL,
    fill = "Time period"
  )


## Other layout

ggplot(posterior_compare,
       aes(x = .value, y = time, fill = time)) +
  stat_halfeye(
    adjust = 0.5,
    .width = c(0.5, 0.8, 0.95),
    alpha = 0.7
  ) +
  facet_grid(parameter ~ ., scales = "free_x") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Posterior effect size",
    y = NULL,
    fill = "Time period"
  )


## X range between -0.25 and 2.5 (because of distance)


ggplot(posterior_compare,
       aes(x = .value, y = time, fill = time)) +
  stat_halfeye(
    adjust = 0.5,
    .width = c(0.5, 0.8, 0.95),
    alpha = 0.7
  ) +
  facet_grid(parameter ~ ., scales = "free_x") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(limits = c(-0.25, 2.5)) +  
  theme_minimal() +
  labs(
    x = "Posterior effect size",
    y = NULL,
    fill = "Time period"
  )




##

# Overlapping Past and Present distributions
ggplot(posterior_compare,
       aes(x = .value, color = time, fill = time)) +
  stat_density(alpha = 0.3, adjust = 1) +
  facet_wrap(~parameter, scales = "free") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Posterior effect size",
    y = "Density",
    color = "Time period",
    fill = "Time period"
  )

