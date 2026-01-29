
## Packages

library(dplyr)

######## Spearman correlations
####
##
#
#------------------------------------------------------------------------

# Corals taxonomic - PRESENT
model_corals_taxonomic$log_distance <- log(model_corals_taxonomic$Distance)
model_corals_taxonomic$log_isolation <- log(model_corals_taxonomic$diff_isolation)

## Scaling variables
coral_vars <- model_corals_taxonomic %>%
  mutate(
    beta_sor_adj = beta_sor - 0.000001,  
    diff_isolation_sc = scale(diff_isolation)[,1],
    diff_reef_area_sc = scale(diff_reef_area)[,1],
    diff_past_isolation_sc = scale(diff_past_isolation)[,1], #quaternary distance
    diff_past_reef_area_sc = scale(diff_past_reef_area)[,1],
    diff_age_sc = scale(diff_age)[,1],
    dist_sc = scale(log_distance)[,1], # log
    diff_present_sst = scale(diff_sst)[,1],
    diff_past_sst = scale(diff_past_sst)
  ) %>%
  select(diff_isolation_sc, diff_reef_area_sc, diff_past_isolation_sc,
         diff_past_reef_area_sc, diff_age_sc, dist_sc, diff_present_sst, diff_past_sst)


# Correlation
cor_matrix <- cor(coral_vars, method = "spearman", use = "pairwise.complete.obs")

cor_matrix

#################################

# Plants taxonomic - PRESENT

plants_vars <- model_plants_taxonomic %>%
  mutate(
    d_past_sc = scale(d_past.y)[,1],
    diff_emerse_past_sc = scale(diff_emerse_past)[,1],
    diff_emerse_pres_sc = scale(diff_emerse_pres)[,1],
    diff_island_age_sc = scale(diff_island_age)[,1],
    distance_km_sc = scale(Distance_km)[,1]
  ) %>%
  select(d_past_sc, diff_emerse_past_sc, diff_emerse_pres_sc,
         diff_island_age_sc, distance_km_sc)


# Correlation
cor_plants_matrix <- cor(plants_vars, method = "spearman", use = "pairwise.complete.obs")

cor_plants_matrix
