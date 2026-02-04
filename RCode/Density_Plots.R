
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
library(stringr)
library(patchwork)



#1. Corals Taxonomic - Past and Pres ----


extract_betas <- function(model, facet, time){
  as_draws_df(model) %>%
    select(starts_with("b_")) %>%
    pivot_longer(
      cols = everything(),
      names_to = "parameter",
      values_to = "value"
    ) %>%
    mutate(
      parameter = str_remove(parameter, "^b_"),
      facet = facet,
      time = time
    )
}


## Posterior PRESENT

post_present <- bind_rows(
  extract_betas(model_beta_corals_taxonomic_present_bayes_all,
                facet = "Total beta", time = "Present"),
  
  extract_betas(model_beta_corals_taxonomic_present_bayes_nest,
                facet = "Nestedness", time = "Present"),
  
  extract_betas(model_beta_corals_taxonomic_present_bayes_turn,
                facet = "Turnover", time = "Present")
)


## Posterior PAST

post_past <- bind_rows(
  extract_betas(model_beta_corals_taxonomic_past_bayes_all,
                facet = "Total beta", time = "Past"),
  
  extract_betas(model_beta_corals_taxonomic_past_bayes_nest,
                facet = "Nestedness", time = "Past"),
  
  extract_betas(model_beta_corals_taxonomic_past_bayes_turn,
                facet = "Turnover", time = "Past")
)



## Density Plot - PRESENT

corals_taxonomic_pres <- ggplot(
  filter(post_present, parameter != "Intercept"),
  aes(x = value, fill = facet)
) +
  geom_density(alpha = 0.6, color = NA) +
  facet_grid(parameter ~ facet, scales = "free_x") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Posterior estimate",
    y = "Density",
    title = "Coral taxonomic beta diversity – Present"
  ) +
  theme(legend.position = "none")



## Density Plot - PAST

corals_taxonomic_past <- ggplot(
  filter(post_past, parameter != "Intercept"),
  aes(x = value, fill = facet)
) +
  geom_density(alpha = 0.6, color = NA) +
  facet_grid(parameter ~ facet, scales = "free_x") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Posterior estimate",
    y = "Density",
    title = "Coral taxonomic beta diversity – Past"
  ) +
  theme(legend.position = "none")


###

post_all <- bind_rows(post_past, post_present)


##

ggplot(
  filter(post_all, parameter != "Intercept"),
  aes(x = value, fill = time)
) +
  geom_density(alpha = 0.6, color = NA) +
  facet_grid(parameter ~ facet, scales = "free_x") +
  geom_vline(xintercept = 0, linetype = "dashed", inherit.aes = FALSE) +
  scale_fill_manual(values = c("Past" = "#1f78b4", "Present" = "#e31a1c")) +
  theme_minimal() +
  labs(
    x = "Posterior estimate",
    y = "Density",
    fill = "Time",
    title = "Coral taxonomic beta diversity – Past vs Present"
  )


##

corals_taxonomic_past + corals_taxonomic_pres + plot_layout(ncol = 2)




#2. Fish Taxonomic - Past and Pres ----

extract_betas <- function(model, facet, time){
  as_draws_df(model) %>%
    select(starts_with("b_")) %>%
    pivot_longer(
      cols = everything(),
      names_to = "parameter",
      values_to = "value"
    ) %>%
    mutate(
      parameter = str_remove(parameter, "^b_"),
      facet = facet,
      time = time
    )
}


## Posterior PRESENT

post_present <- bind_rows(
  extract_betas(model_beta_fish_taxonomic_present_bayes_all,
                facet = "Total beta", time = "Present"),
  
  extract_betas(model_beta_fish_taxonomic_present_bayes_nest,
                facet = "Nestedness", time = "Present"),
  
  extract_betas(model_beta_fish_taxonomic_present_bayes_turn,
                facet = "Turnover", time = "Present")
)


## Posterior PAST

post_past <- bind_rows(
  extract_betas(model_beta_fish_taxonomic_past_bayes_all,
                facet = "Total beta", time = "Past"),
  
  extract_betas(model_beta_fish_taxonomic_past_bayes_nest,
                facet = "Nestedness", time = "Past"),
  
  extract_betas(model_beta_fish_taxonomic_past_bayes_turn,
                facet = "Turnover", time = "Past")
)



## Density Plot - PRESENT

fish_taxonomic_pres <- ggplot(
  filter(post_present, parameter != "Intercept"),
  aes(x = value, fill = facet)
) +
  geom_density(alpha = 0.6, color = NA) +
  facet_grid(parameter ~ facet, scales = "free_x") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Posterior estimate",
    y = "Density",
    title = "fish taxonomic beta diversity – Present"
  ) +
  theme(legend.position = "none")



## Density Plot - PAST

fish_taxonomic_past <- ggplot(
  filter(post_past, parameter != "Intercept"),
  aes(x = value, fill = facet)
) +
  geom_density(alpha = 0.6, color = NA) +
  facet_grid(parameter ~ facet, scales = "free_x") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Posterior estimate",
    y = "Density",
    title = "fish taxonomic beta diversity – Past"
  ) +
  theme(legend.position = "none")


###

post_all <- bind_rows(post_past, post_present)


##

ggplot(
  filter(post_all, parameter != "Intercept"),
  aes(x = value, fill = time)
) +
  geom_density(alpha = 0.6, color = NA) +
  facet_grid(parameter ~ facet, scales = "free_x") +
  geom_vline(xintercept = 0, linetype = "dashed", inherit.aes = FALSE) +
  scale_fill_manual(values = c("Past" = "#1f78b4", "Present" = "#e31a1c")) +
  theme_minimal() +
  labs(
    x = "Posterior estimate",
    y = "Density",
    fill = "Time",
    title = "fish taxonomic beta diversity – Past vs Present"
  )


##

fish_taxonomic_past + fish_taxonomic_pres + plot_layout(ncol = 2)


###### Plotting corals and fish together 

#3. Corals and fish together ----

extract_betas <- function(model, facet, time){
  as_draws_df(model) %>%
    select(starts_with("b_")) %>%
    pivot_longer(
      cols = everything(),
      names_to = "parameter",
      values_to = "value"
    ) %>%
    mutate(
      parameter = str_remove(parameter, "^b_"),
      parameter = str_replace(parameter, "diff_past_", "diff_"),
      parameter = str_replace(parameter, "diff_present_", "diff_"),
      facet = facet,
      time = time
    )
}


## Corals

# CORALS - Present
post_present_corals <- bind_rows(
  extract_betas(model_beta_corals_taxonomic_present_bayes_all,  "Total beta",  "Present"),
  extract_betas(model_beta_corals_taxonomic_present_bayes_nest, "Nestedness", "Present"),
  extract_betas(model_beta_corals_taxonomic_present_bayes_turn, "Turnover",  "Present")
)

# CORALS - Past
post_past_corals <- bind_rows(
  extract_betas(model_beta_corals_taxonomic_past_bayes_all,  "Total beta",  "Past"),
  extract_betas(model_beta_corals_taxonomic_past_bayes_nest, "Nestedness", "Past"),
  extract_betas(model_beta_corals_taxonomic_past_bayes_turn, "Turnover",  "Past")
)


## Fish

# FISH - Present
post_present_fish <- bind_rows(
  extract_betas(model_beta_fish_taxonomic_present_bayes_all,  "Total beta",  "Present"),
  extract_betas(model_beta_fish_taxonomic_present_bayes_nest, "Nestedness", "Present"),
  extract_betas(model_beta_fish_taxonomic_present_bayes_turn, "Turnover",  "Present")
)

# FISH - Past
post_past_fish <- bind_rows(
  extract_betas(model_beta_fish_taxonomic_past_bayes_all,  "Total beta",  "Past"),
  extract_betas(model_beta_fish_taxonomic_past_bayes_nest, "Nestedness", "Past"),
  extract_betas(model_beta_fish_taxonomic_past_bayes_turn, "Turnover",  "Past")
)


# Order

facet_levels <- c("Total beta", "Nestedness", "Turnover")

post_past_corals <- post_past_corals %>% mutate(facet = factor(facet, levels = facet_levels))
post_present_corals <- post_present_corals %>% mutate(facet = factor(facet, levels = facet_levels))

post_past_fish <- post_past_fish %>% mutate(facet = factor(facet, levels = facet_levels))
post_present_fish <- post_present_fish %>% mutate(facet = factor(facet, levels = facet_levels))


# Final object

post_all <- bind_rows(
  post_past_corals %>% mutate(taxon = "Corals"),
  post_present_corals %>% mutate(taxon = "Corals"),
  post_past_fish %>% mutate(taxon = "Fish"),
  post_present_fish %>% mutate(taxon = "Fish")
) %>%
  mutate(
    taxon = factor(taxon, levels = c("Corals", "Fish")),
    time  = factor(time,  levels = c("Past", "Present"))
  )


####
ggplot(
  filter(post_all, parameter != "Intercept"),
  aes(x = value, fill = taxon)
) +
  geom_density(alpha = .9, color = NA) +
  facet_grid(parameter + taxon ~ time + facet, scales = "free_x") +
  geom_vline(xintercept = 0, linetype = "dashed", inherit.aes = FALSE) +
  scale_fill_manual(values = c("Corals" = "#2364aa", "Fish" = "#a9d3ff")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text.y.left = element_blank(),   
    strip.background.y = element_blank()   
  ) +
  labs(
    x = "Posterior estimate",
    y = "Density",
    title = "Taxonomic beta diversity – Corals vs Fish"
  )


## X axis between -1 and 1.5

marine_tax_plot <- ggplot(
  filter(post_all, parameter != "Intercept"),
  aes(x = value, fill = taxon)
) +
  geom_density(alpha = .9, color = NA) +
  facet_grid(parameter + taxon ~ time + facet, scales = "free_x") +
  geom_vline(xintercept = 0, linetype = "dashed", inherit.aes = FALSE) +
  scale_fill_manual(values = c("Corals" = "#2364aa", "Fish" = "#a9d3ff")) +
  scale_x_continuous(limits = c(-1, 1.5)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text.y.left = element_blank(),
    strip.background.y = element_blank()
  ) +
  labs(
    x = "Posterior estimate",
    y = "Density",
    title = "Taxonomic beta diversity – Corals vs Fish"
  )

# Saving image and creating 'Output' folder

if (!dir.exists("Output")) {
  dir.create("Output")
}

#
ggsave(
  filename = "Output/marine_tax_plot.pdf",
  plot = marine_plot,
  device = "pdf",
  width = 12,
  height = 8
)


## WITHOUT ARCH DIST

marine_tax_plot_without_ArchDis <- ggplot(
  filter(post_all, parameter != "Intercept", parameter != "dist_sc"),
  aes(x = value, fill = taxon)
) +
  geom_density(alpha = .9, color = NA) +
  facet_grid(parameter + taxon ~ time + facet, scales = "free_x") +
  geom_vline(xintercept = 0, linetype = "dashed", inherit.aes = FALSE) +
  scale_fill_manual(values = c("Corals" = "#2364aa", "Fish" = "#a9d3ff")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text.y.left = element_blank(),
    strip.background.y = element_blank()
  ) +
  labs(
    x = "Posterior estimate",
    y = "Density",
    title = "Taxonomic beta diversity – Corals vs Fish (sem dist_sc)"
  )



#
ggsave(
  filename = "Output/marine_tax_plot_without_ArchDis.pdf",
  plot = marine_plot,
  device = "pdf",
  width = 12,
  height = 8
)



#4. Corals Functional - Past and Pres ----


extract_betas <- function(model, facet, time){
  as_draws_df(model) %>%
    select(starts_with("b_")) %>%
    pivot_longer(
      cols = everything(),
      names_to = "parameter",
      values_to = "value"
    ) %>%
    mutate(
      parameter = str_remove(parameter, "^b_"),
      facet = facet,
      time = time
    )
}


## Posterior PRESENT

post_present <- bind_rows(
  extract_betas(model_beta_corals_functional_present_bayes_all,
                facet = "Total beta", time = "Present"),
  
  extract_betas(model_beta_corals_functional_present_bayes_nest,
                facet = "Nestedness", time = "Present"),
  
  extract_betas(model_beta_corals_functional_present_bayes_turn,
                facet = "Turnover", time = "Present")
)


## Posterior PAST

post_past <- bind_rows(
  extract_betas(model_beta_corals_functional_past_bayes_all,
                facet = "Total beta", time = "Past"),
  
  extract_betas(model_beta_corals_functional_past_bayes_nest,
                facet = "Nestedness", time = "Past"),
  
  extract_betas(model_beta_corals_functional_past_bayes_turn,
                facet = "Turnover", time = "Past")
)



## Density Plot - PRESENT

corals_functional_pres <- ggplot(
  filter(post_present, parameter != "Intercept"),
  aes(x = value, fill = facet)
) +
  geom_density(alpha = 0.6, color = NA) +
  facet_grid(parameter ~ facet, scales = "free_x") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Posterior estimate",
    y = "Density",
    title = "Coral functional beta diversity – Present"
  ) +
  theme(legend.position = "none")



## Density Plot - PAST

corals_functional_past <- ggplot(
  filter(post_past, parameter != "Intercept"),
  aes(x = value, fill = facet)
) +
  geom_density(alpha = 0.6, color = NA) +
  facet_grid(parameter ~ facet, scales = "free_x") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Posterior estimate",
    y = "Density",
    title = "Coral functional beta diversity – Past"
  ) +
  theme(legend.position = "none")


###

post_all <- bind_rows(post_past, post_present)


##

ggplot(
  filter(post_all, parameter != "Intercept"),
  aes(x = value, fill = time)
) +
  geom_density(alpha = 0.6, color = NA) +
  facet_grid(parameter ~ facet, scales = "free_x") +
  geom_vline(xintercept = 0, linetype = "dashed", inherit.aes = FALSE) +
  scale_fill_manual(values = c("Past" = "#1f78b4", "Present" = "#e31a1c")) +
  theme_minimal() +
  labs(
    x = "Posterior estimate",
    y = "Density",
    fill = "Time",
    title = "Coral functional beta diversity – Past vs Present"
  )



##

corals_functional_past + corals_functional_pres + plot_layout(ncol = 2)



#5. Fish Functional - Past and Pres ----

extract_betas <- function(model, facet, time){
  as_draws_df(model) %>%
    select(starts_with("b_")) %>%
    pivot_longer(
      cols = everything(),
      names_to = "parameter",
      values_to = "value"
    ) %>%
    mutate(
      parameter = str_remove(parameter, "^b_"),
      facet = facet,
      time = time
    )
}


## Posterior PRESENT

post_present <- bind_rows(
  extract_betas(model_beta_fish_functional_present_bayes_all,
                facet = "Total beta", time = "Present"),
  
  extract_betas(model_beta_fish_functional_present_bayes_nest,
                facet = "Nestedness", time = "Present"),
  
  extract_betas(model_beta_fish_functional_present_bayes_turn,
                facet = "Turnover", time = "Present")
)


## Posterior PAST

post_past <- bind_rows(
  extract_betas(model_beta_fish_functional_past_bayes_all,
                facet = "Total beta", time = "Past"),
  
  extract_betas(model_beta_fish_functional_past_bayes_nest,
                facet = "Nestedness", time = "Past"),
  
  extract_betas(model_beta_fish_functional_past_bayes_turn,
                facet = "Turnover", time = "Past")
)



## Density Plot - PRESENT

fish_functional_pres <- ggplot(
  filter(post_present, parameter != "Intercept"),
  aes(x = value, fill = facet)
) +
  geom_density(alpha = 0.6, color = NA) +
  facet_grid(parameter ~ facet, scales = "free_x") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Posterior estimate",
    y = "Density",
    title = "fish functional beta diversity – Present"
  ) +
  theme(legend.position = "none")



## Density Plot - PAST

fish_functional_past <- ggplot(
  filter(post_past, parameter != "Intercept"),
  aes(x = value, fill = facet)
) +
  geom_density(alpha = 0.6, color = NA) +
  facet_grid(parameter ~ facet, scales = "free_x") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Posterior estimate",
    y = "Density",
    title = "fish functional beta diversity – Past"
  ) +
  theme(legend.position = "none")


###

post_all <- bind_rows(post_past, post_present)


##

ggplot(
  filter(post_all, parameter != "Intercept"),
  aes(x = value, fill = time)
) +
  geom_density(alpha = 0.6, color = NA) +
  facet_grid(parameter ~ facet, scales = "free_x") +
  geom_vline(xintercept = 0, linetype = "dashed", inherit.aes = FALSE) +
  scale_fill_manual(values = c("Past" = "#1f78b4", "Present" = "#e31a1c")) +
  theme_minimal() +
  labs(
    x = "Posterior estimate",
    y = "Density",
    fill = "Time",
    title = "fish functional beta diversity – Past vs Present"
  )


##

fish_functional_past + fish_functional_pres + plot_layout(ncol = 2)


###### Plotting corals and fish together 

#6. Corals and fish together ----

extract_betas <- function(model, facet, time){
  as_draws_df(model) %>%
    select(starts_with("b_")) %>%
    pivot_longer(
      cols = everything(),
      names_to = "parameter",
      values_to = "value"
    ) %>%
    mutate(
      parameter = str_remove(parameter, "^b_"),
      parameter = str_replace(parameter, "diff_past_", "diff_"),
      parameter = str_replace(parameter, "diff_present_", "diff_"),
      facet = facet,
      time = time
    )
}


## Corals

# CORALS - Present
post_present_corals <- bind_rows(
  extract_betas(model_beta_corals_functional_present_bayes_all,  "Total beta",  "Present"),
  extract_betas(model_beta_corals_functional_present_bayes_nest, "Nestedness", "Present"),
  extract_betas(model_beta_corals_functional_present_bayes_turn, "Turnover",  "Present")
)

# CORALS - Past
post_past_corals <- bind_rows(
  extract_betas(model_beta_corals_functional_past_bayes_all,  "Total beta",  "Past"),
  extract_betas(model_beta_corals_functional_past_bayes_nest, "Nestedness", "Past"),
  extract_betas(model_beta_corals_functional_past_bayes_turn, "Turnover",  "Past")
)


## Fish

# FISH - Present
post_present_fish <- bind_rows(
  extract_betas(model_beta_fish_functional_present_bayes_all,  "Total beta",  "Present"),
  extract_betas(model_beta_fish_functional_present_bayes_nest, "Nestedness", "Present"),
  extract_betas(model_beta_fish_functional_present_bayes_turn, "Turnover",  "Present")
)

# FISH - Past
post_past_fish <- bind_rows(
  extract_betas(model_beta_fish_functional_past_bayes_all,  "Total beta",  "Past"),
  extract_betas(model_beta_fish_functional_past_bayes_nest, "Nestedness", "Past"),
  extract_betas(model_beta_fish_functional_past_bayes_turn, "Turnover",  "Past")
)


# Order

facet_levels <- c("Total beta", "Nestedness", "Turnover")

post_past_corals <- post_past_corals %>% mutate(facet = factor(facet, levels = facet_levels))
post_present_corals <- post_present_corals %>% mutate(facet = factor(facet, levels = facet_levels))

post_past_fish <- post_past_fish %>% mutate(facet = factor(facet, levels = facet_levels))
post_present_fish <- post_present_fish %>% mutate(facet = factor(facet, levels = facet_levels))


# Final object

post_all <- bind_rows(
  post_past_corals %>% mutate(taxon = "Corals"),
  post_present_corals %>% mutate(taxon = "Corals"),
  post_past_fish %>% mutate(taxon = "Fish"),
  post_present_fish %>% mutate(taxon = "Fish")
) %>%
  mutate(
    taxon = factor(taxon, levels = c("Corals", "Fish")),
    time  = factor(time,  levels = c("Past", "Present"))
  )


####
ggplot(
  filter(post_all, parameter != "Intercept"),
  aes(x = value, fill = taxon)
) +
  geom_density(alpha = 0.9, color = NA) +
  facet_grid(parameter + taxon ~ time + facet, scales = "free_x") +
  geom_vline(xintercept = 0, linetype = "dashed", inherit.aes = FALSE) +
  scale_fill_manual(values = c("Corals" = "#2364aa", "Fish" = "#a9d3ff")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text.y.left = element_blank(),  
    strip.background.y = element_blank()  
  ) +
  labs(
    x = "Posterior estimate",
    y = "Density",
    title = "Functional beta diversity – Corals vs Fish"
  )


## WITHOUT ARCH DIST

marine_func_plot_without_ArchDis.pdf <- ggplot(
  filter(post_all, parameter != "Intercept", parameter != "dist_sc"),
  aes(x = value, fill = taxon)
) +
  geom_density(alpha = .9, color = NA) +
  facet_grid(parameter + taxon ~ time + facet, scales = "free_x") +
  geom_vline(xintercept = 0, linetype = "dashed", inherit.aes = FALSE) +
  scale_fill_manual(values = c("Corals" = "#2364aa", "Fish" = "#a9d3ff")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text.y.left = element_blank(),
    strip.background.y = element_blank()
  ) +
  labs(
    x = "Posterior estimate",
    y = "Density",
    title = "Taxonomic beta diversity – Corals vs Fish (sem dist_sc)"
  )


#
ggsave(
  filename = "Output/marine_func_plot_without_ArchDis.pdf",
  plot = marine_plot,
  device = "pdf",
  width = 12,
  height = 8
)


## X axis between -1 and 1.5

marine_func_plot <- ggplot(
  filter(post_all, parameter != "Intercept"),
  aes(x = value, fill = taxon)
) +
  geom_density(alpha = .9, color = NA) +
  facet_grid(parameter + taxon ~ time + facet, scales = "free_x") +
  geom_vline(xintercept = 0, linetype = "dashed", inherit.aes = FALSE) +
  scale_fill_manual(values = c("Corals" = "#2364aa", "Fish" = "#a9d3ff")) +
  scale_x_continuous(limits = c(-2, 1)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text.y.left = element_blank(),
    strip.background.y = element_blank()
  ) +
  labs(
    x = "Posterior estimate",
    y = "Density",
    title = "Functional beta diversity – Corals vs Fish"
  )


#
ggsave(
  filename = "Output/marine_func_plot.pdf",
  plot = marine_plot,
  device = "pdf",
  width = 12,
  height = 8
)



#7. Plants Taxonomic - Past and Pres ----


extract_betas <- function(model, facet, time){
  as_draws_df(model) %>%
    select(starts_with("b_")) %>%
    pivot_longer(
      cols = everything(),
      names_to = "parameter",
      values_to = "value"
    ) %>%
    mutate(
      parameter = str_remove(parameter, "^b_"),
      facet = facet,
      time = time
    )
}


## Posterior PRESENT

post_present <- bind_rows(
  extract_betas(model_beta_plants_taxonomic_present,
                facet = "Total beta", time = "Present"),
  
  extract_betas(model_beta_plants_taxonomic_present_nest,
                facet = "Nestedness", time = "Present"),
  
  extract_betas(model_beta_plants_taxonomic_present_turn,
                facet = "Turnover", time = "Present")
)


## Posterior PAST

post_past <- bind_rows(
  extract_betas(model_beta_plants_taxonomic_past,
                facet = "Total beta", time = "Past"),
  
  extract_betas(model_beta_plants_taxonomic_past_nest,
                facet = "Nestedness", time = "Past"),
  
  extract_betas(model_beta_plants_taxonomic_past_turn,
                facet = "Turnover", time = "Past")
)



## Density Plot - PRESENT

plants_taxonomic_pres <- ggplot(
  filter(post_present, parameter != "Intercept"),
  aes(x = value, fill = facet)
) +
  geom_density(alpha = 0.6, color = NA) +
  facet_grid(parameter ~ facet, scales = "free_x") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Posterior estimate",
    y = "Density",
    title = "Plants taxonomic beta diversity – Present"
  ) +
  theme(legend.position = "none")



## Density Plot - PAST

plants_taxonomic_past <- ggplot(
  filter(post_past, parameter != "Intercept"),
  aes(x = value, fill = facet)
) +
  geom_density(alpha = 0.6, color = NA) +
  facet_grid(parameter ~ facet, scales = "free_x") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Posterior estimate",
    y = "Density",
    title = "Plants taxonomic beta diversity – Past"
  ) +
  theme(legend.position = "none")


###

post_all <- bind_rows(post_past, post_present)


##

ggplot(
  filter(post_all, parameter != "Intercept"),
  aes(x = value, fill = time)
) +
  geom_density(alpha = 0.6, color = NA) +
  facet_grid(parameter ~ facet, scales = "free_x") +
  geom_vline(xintercept = 0, linetype = "dashed", inherit.aes = FALSE) +
  scale_fill_manual(values = c("Past" = "#1f78b4", "Present" = "#e31a1c")) +
  theme_minimal() +
  labs(
    x = "Posterior estimate",
    y = "Density",
    fill = "Time",
    title = "Plants taxonomic beta diversity – Past vs Present"
  )



##

plants_taxonomic_past + plants_taxonomic_pres + plot_layout(ncol = 2)



#8. Birds Taxonomic - Past and Pres ----

extract_betas <- function(model, facet, time){
  as_draws_df(model) %>%
    select(starts_with("b_")) %>%
    pivot_longer(
      cols = everything(),
      names_to = "parameter",
      values_to = "value"
    ) %>%
    mutate(
      parameter = str_remove(parameter, "^b_"),
      facet = facet,
      time = time
    )
}


## Posterior PRESENT

post_present <- bind_rows(
  extract_betas(model_beta_birds_taxonomic_present,
                facet = "Total beta", time = "Present"),
  
  extract_betas(model_beta_birds_taxonomic_present,
                facet = "Nestedness", time = "Present"),
  
  extract_betas(model_beta_birds_taxonomic_present,
                facet = "Turnover", time = "Present")
)


## Posterior PAST

post_past <- bind_rows(
  extract_betas(model_beta_birds_taxonomic_past,
                facet = "Total beta", time = "Past"),
  
  extract_betas(model_beta_birds_taxonomic_past,
                facet = "Nestedness", time = "Past"),
  
  extract_betas(model_beta_birds_taxonomic_past,
                facet = "Turnover", time = "Past")
)



## Density Plot - PRESENT

birds_taxonomic_pres <- ggplot(
  filter(post_present, parameter != "Intercept"),
  aes(x = value, fill = facet)
) +
  geom_density(alpha = 0.6, color = NA) +
  facet_grid(parameter ~ facet, scales = "free_x") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Posterior estimate",
    y = "Density",
    title = "birds taxonomic beta diversity – Present"
  ) +
  theme(legend.position = "none")



## Density Plot - PAST

birds_taxonomic_past <- ggplot(
  filter(post_past, parameter != "Intercept"),
  aes(x = value, fill = facet)
) +
  geom_density(alpha = 0.6, color = NA) +
  facet_grid(parameter ~ facet, scales = "free_x") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Posterior estimate",
    y = "Density",
    title = "birds taxonomic beta diversity – Past"
  ) +
  theme(legend.position = "none")


###

post_all <- bind_rows(post_past, post_present)


##

ggplot(
  filter(post_all, parameter != "Intercept"),
  aes(x = value, fill = time)
) +
  geom_density(alpha = 0.6, color = NA) +
  facet_grid(parameter ~ facet, scales = "free_x") +
  geom_vline(xintercept = 0, linetype = "dashed", inherit.aes = FALSE) +
  scale_fill_manual(values = c("Past" = "#1f78b4", "Present" = "#e31a1c")) +
  theme_minimal() +
  labs(
    x = "Posterior estimate",
    y = "Density",
    fill = "Time",
    title = "birds taxonomic beta diversity – Past vs Present"
  )


##

birds_taxonomic_past + birds_taxonomic_pres + plot_layout(ncol = 2)


###### Plotting plants and birds together 

#9. Plants and birds together ----

extract_betas <- function(model, facet, time){
  as_draws_df(model) %>%
    select(starts_with("b_")) %>%
    pivot_longer(
      cols = everything(),
      names_to = "parameter_raw",
      values_to = "value"
    ) %>%
    mutate(
      parameter = str_remove(parameter_raw, "^b_"),
      
      # PADRONIZANDO nomes
      parameter = case_when(
        parameter %in% c("diff_emerse_pres_sc", "diff_emerse_past_sc") ~ "diff_emerse_sc",
        parameter %in% c("distance_km_sc", "d_past_sc") ~ "distance_sc",
        parameter %in% c("diff_present_air", "diff_past_air") ~ "diff_air",
        TRUE ~ parameter
      ),
      
      facet = facet,
      time  = time
    )
}


## Plants

post_present_plants <- bind_rows(
  extract_betas(model_beta_plants_taxonomic_present, "Total beta", "Present"),
  extract_betas(model_beta_plants_taxonomic_present_nest, "Nestedness", "Present"),
  extract_betas(model_beta_plants_taxonomic_present_turn, "Turnover", "Present")
)

post_past_plants <- bind_rows(
  extract_betas(model_beta_plants_taxonomic_past, "Total beta", "Past"),
  extract_betas(model_beta_plants_taxonomic_past_nest, "Nestedness", "Past"),
  extract_betas(model_beta_plants_taxonomic_past_turn, "Turnover", "Past")
)

## Birds
post_present_birds <- bind_rows(
  extract_betas(model_beta_birds_taxonomic_present, "Total beta", "Present"),
  extract_betas(model_beta_birds_taxonomic_present_nest, "Nestedness", "Present"),
  extract_betas(model_beta_birds_taxonomic_present_turn, "Turnover", "Present")
)

post_past_birds <- bind_rows(
  extract_betas(model_beta_birds_taxonomic_past, "Total beta", "Past"),
  extract_betas(model_beta_birds_taxonomic_past_nest, "Nestedness", "Past"),
  extract_betas(model_beta_birds_taxonomic_past_turn, "Turnover", "Past")
)

## both
post_all <- bind_rows(
  post_past_plants %>% mutate(taxon = "Plants"),
  post_present_plants %>% mutate(taxon = "Plants"),
  post_past_birds %>% mutate(taxon = "Birds"),
  post_present_birds %>% mutate(taxon = "Birds")
) %>%
  mutate(
    taxon = factor(taxon, levels = c("Plants", "Birds")),
    time  = factor(time, levels = c("Past", "Present"))
  )

##
terrestrial_tax_plot <- ggplot(
  filter(post_all, parameter != "Intercept"),
  aes(x = value, fill = taxon)
) +
  geom_density(alpha = 0.9, color = NA) +
  facet_grid(parameter + taxon ~ time + facet, scales = "free_x") +
  geom_vline(xintercept = 0, linetype = "dashed", inherit.aes = FALSE) +
  scale_fill_manual(values = c("Plants" = "#9c6644", "Birds" = "#d69f7e")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text.y.left = element_blank(),
    strip.background.y = element_blank()
  ) +
  labs(
    x = "Posterior estimate",
    y = "Density",
    title = "taxonomic beta diversity – Plants vs Birds"
  )



#
ggsave(
  filename = "Output/terrestrial_tax_plot.pdf",
  plot = marine_plot,
  device = "pdf",
  width = 12,
  height = 8
)


#10. Plants Functional - Past and Pres ----


extract_betas <- function(model, facet, time){
  as_draws_df(model) %>%
    select(starts_with("b_")) %>%
    pivot_longer(
      cols = everything(),
      names_to = "parameter",
      values_to = "value"
    ) %>%
    mutate(
      parameter = str_remove(parameter, "^b_"),
      facet = facet,
      time = time
    )
}


## Posterior PRESENT

post_present <- bind_rows(
  extract_betas(model_beta_plants_functional_present,
                facet = "Total beta", time = "Present"),
  
  extract_betas(model_beta_plants_functional_present_nest,
                facet = "Nestedness", time = "Present"),
  
  extract_betas(model_beta_plants_functional_present_turn,
                facet = "Turnover", time = "Present")
)


## Posterior PAST

post_past <- bind_rows(
  extract_betas(model_beta_plants_functional_past,
                facet = "Total beta", time = "Past"),
  
  extract_betas(model_beta_plants_functional_past_nest,
                facet = "Nestedness", time = "Past"),
  
  extract_betas(model_beta_plants_functional_past_turn,
                facet = "Turnover", time = "Past")
)



## Density Plot - PRESENT

plants_functional_pres <- ggplot(
  filter(post_present, parameter != "Intercept"),
  aes(x = value, fill = facet)
) +
  geom_density(alpha = 0.6, color = NA) +
  facet_grid(parameter ~ facet, scales = "free_x") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Posterior estimate",
    y = "Density",
    title = "Plants functional beta diversity – Present"
  ) +
  theme(legend.position = "none")



## Density Plot - PAST

plants_functional_past <- ggplot(
  filter(post_past, parameter != "Intercept"),
  aes(x = value, fill = facet)
) +
  geom_density(alpha = 0.6, color = NA) +
  facet_grid(parameter ~ facet, scales = "free_x") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Posterior estimate",
    y = "Density",
    title = "Plants functional beta diversity – Past"
  ) +
  theme(legend.position = "none")


###

post_all <- bind_rows(post_past, post_present)


##

ggplot(
  filter(post_all, parameter != "Intercept"),
  aes(x = value, fill = time)
) +
  geom_density(alpha = 0.6, color = NA) +
  facet_grid(parameter ~ facet, scales = "free_x") +
  geom_vline(xintercept = 0, linetype = "dashed", inherit.aes = FALSE) +
  scale_fill_manual(values = c("Past" = "#1f78b4", "Present" = "#e31a1c")) +
  theme_minimal() +
  labs(
    x = "Posterior estimate",
    y = "Density",
    fill = "Time",
    title = "Plants functional beta diversity – Past vs Present"
  )



##

plants_functional_past + plants_functional_pres + plot_layout(ncol = 2)



#11. Birds Functional - Past and Pres ----

extract_betas <- function(model, facet, time){
  as_draws_df(model) %>%
    select(starts_with("b_")) %>%
    pivot_longer(
      cols = everything(),
      names_to = "parameter",
      values_to = "value"
    ) %>%
    mutate(
      parameter = str_remove(parameter, "^b_"),
      facet = facet,
      time = time
    )
}


## Posterior PRESENT

post_present <- bind_rows(
  extract_betas(model_beta_birds_functional_present,
                facet = "Total beta", time = "Present"),
  
  extract_betas(model_beta_birds_functional_present,
                facet = "Nestedness", time = "Present"),
  
  extract_betas(model_beta_birds_functional_present,
                facet = "Turnover", time = "Present")
)


## Posterior PAST

post_past <- bind_rows(
  extract_betas(model_beta_birds_functional_past,
                facet = "Total beta", time = "Past"),
  
  extract_betas(model_beta_birds_functional_past,
                facet = "Nestedness", time = "Past"),
  
  extract_betas(model_beta_birds_functional_past,
                facet = "Turnover", time = "Past")
)



## Density Plot - PRESENT

birds_functional_pres <- ggplot(
  filter(post_present, parameter != "Intercept"),
  aes(x = value, fill = facet)
) +
  geom_density(alpha = 0.6, color = NA) +
  facet_grid(parameter ~ facet, scales = "free_x") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Posterior estimate",
    y = "Density",
    title = "birds functional beta diversity – Present"
  ) +
  theme(legend.position = "none")



## Density Plot - PAST

birds_functional_past <- ggplot(
  filter(post_past, parameter != "Intercept"),
  aes(x = value, fill = facet)
) +
  geom_density(alpha = 0.6, color = NA) +
  facet_grid(parameter ~ facet, scales = "free_x") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Posterior estimate",
    y = "Density",
    title = "birds functional beta diversity – Past"
  ) +
  theme(legend.position = "none")


###

post_all <- bind_rows(post_past, post_present)


##

ggplot(
  filter(post_all, parameter != "Intercept"),
  aes(x = value, fill = time)
) +
  geom_density(alpha = 0.6, color = NA) +
  facet_grid(parameter ~ facet, scales = "free_x") +
  geom_vline(xintercept = 0, linetype = "dashed", inherit.aes = FALSE) +
  scale_fill_manual(values = c("Past" = "#1f78b4", "Present" = "#e31a1c")) +
  theme_minimal() +
  labs(
    x = "Posterior estimate",
    y = "Density",
    fill = "Time",
    title = "birds functional beta diversity – Past vs Present"
  )


##

birds_functional_past + birds_functional_pres + plot_layout(ncol = 2)


###### Plotting plants and birds together 

#12. Plants and birds together ----

extract_betas <- function(model, facet, time){
  as_draws_df(model) %>%
    select(starts_with("b_")) %>%
    pivot_longer(
      cols = everything(),
      names_to = "parameter_raw",
      values_to = "value"
    ) %>%
    mutate(
      parameter = str_remove(parameter_raw, "^b_"),
      
      # PADRONIZANDO nomes
      parameter = case_when(
        parameter %in% c("diff_emerse_pres_sc", "diff_emerse_past_sc") ~ "diff_emerse_sc",
        parameter %in% c("distance_km_sc", "d_past_sc") ~ "distance_sc",
        parameter %in% c("diff_present_air", "diff_past_air") ~ "diff_air",
        TRUE ~ parameter
      ),
      
      facet = facet,
      time  = time
    )
}


## Plants

post_present_plants <- bind_rows(
  extract_betas(model_beta_plants_functional_present, "Total beta", "Present"),
  extract_betas(model_beta_plants_functional_present_nest, "Nestedness", "Present"),
  extract_betas(model_beta_plants_functional_present_turn, "Turnover", "Present")
)

post_past_plants <- bind_rows(
  extract_betas(model_beta_plants_functional_past, "Total beta", "Past"),
  extract_betas(model_beta_plants_functional_past_nest, "Nestedness", "Past"),
  extract_betas(model_beta_plants_functional_past_turn, "Turnover", "Past")
)

## Birds
post_present_birds <- bind_rows(
  extract_betas(model_beta_birds_functional_present, "Total beta", "Present"),
  extract_betas(model_beta_birds_functional_present_nest, "Nestedness", "Present"),
  extract_betas(model_beta_birds_functional_present_turn, "Turnover", "Present")
)

post_past_birds <- bind_rows(
  extract_betas(model_beta_birds_functional_past, "Total beta", "Past"),
  extract_betas(model_beta_birds_functional_past_nest, "Nestedness", "Past"),
  extract_betas(model_beta_birds_functional_past_turn, "Turnover", "Past")
)

## both
post_all <- bind_rows(
  post_past_plants %>% mutate(taxon = "Plants"),
  post_present_plants %>% mutate(taxon = "Plants"),
  post_past_birds %>% mutate(taxon = "Birds"),
  post_present_birds %>% mutate(taxon = "Birds")
) %>%
  mutate(
    taxon = factor(taxon, levels = c("Plants", "Birds")),
    time  = factor(time, levels = c("Past", "Present"))
  )

##
terrestrial_func_plot <- ggplot(
  filter(post_all, parameter != "Intercept"),
  aes(x = value, fill = taxon)
) +
  geom_density(alpha = 0.9, color = NA) +
  facet_grid(parameter + taxon ~ time + facet, scales = "free_x") +
  geom_vline(xintercept = 0, linetype = "dashed", inherit.aes = FALSE) +
  scale_fill_manual(values = c("Plants" = "#9c6644", "Birds" = "#d69f7e")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text.y.left = element_blank(),
    strip.background.y = element_blank()
  ) +
  labs(
    x = "Posterior estimate",
    y = "Density",
    title = "Functional beta diversity – Plants vs Birds"
  )


ggsave(
  filename = "Output/terrestrial_func_plot.pdf",
  plot = marine_plot,
  device = "pdf",
  width = 12,
  height = 8
)

