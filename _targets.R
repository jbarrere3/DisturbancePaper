#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name _targets.R  
#' @description R script to launch the target pipeline
#' @author Julien BARRERE
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Options and packages ----------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Load targets
library(targets)
# Load functions
lapply(grep("R$", list.files("R"), value = TRUE), function(x) source(file.path("R", x)))
# install if needed and load packages
packages.in <- c("dplyr", "ggplot2", "RCurl", "httr", "tidyr", "data.table", "sp", "R2jags", "rstan", "cowplot",
                 "ggmcmc", "taxize", "rnaturalearth", "ggspatial", "sf", "ggnewscale", "readxl", "scales", 
                 "FactoMineR", "ade4", "factoextra", "xtable", "MASS", "vegan", "lme4", "car", "GGally", "grid")
for(i in 1:length(packages.in)) if(!(packages.in[i] %in% rownames(installed.packages()))) install.packages(packages.in[i])
# Targets options
options(tidyverse.quiet = TRUE, 
        clustermq.scheduler = "multiprocess", 
        dplyr.summarise.inform = FALSE)
tar_option_set(packages = packages.in)
set.seed(2)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Targets workflow --------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

list(
  ## - data files
  tar_target(wood.density_file, "data/traits/GlobalWoodDensityDatabase.xls", format = "file"),
  tar_target(shade.tolerance_file, "data/traits/shade_tolerance_FrenchNFI.csv", format = "file"),
  tar_target(root.depth_file, "data/traits/GRooTAggregateSpeciesVersion.csv", format = "file"),
  tar_target(NFI.traits_file, "data/traits/traitsNFI.csv", format = "file"),
  tar_target(TRY_file, "data/traits/TRY_data_request_21092.txt", format = "file"),
  tar_target(gbif_file, "data/climate/sp_gbif_climate.csv", format = "file"),
  tar_target(gbif_disturbance_file, "data/climate/sp_gbif_disturbance.csv", format = "file"),
  tar_target(sensitivity_file, "data/sensitivity/jags_dominance.Rdata", format = "file"),
  
  ## - Compile traits data
  tar_target(traits_TRY, compile_traits_TRY(TRY_file, NFI.traits_file)),
  tar_target(traits, compile_traits(wood.density_file, shade.tolerance_file, root.depth_file, NFI.traits_file)), 
  
  ## - Get parameters from raw model outputs
  tar_target(param_per_iteration, get_param_from_rdata(sensitivity_file)),
  
  ## - Calculate sensitivity from parameters
  tar_target(sensitivity, get_sensitivity_from_param(param_per_iteration)), 
  
  ## - Plot parameters per species (fig. 2)
  tar_target(fig_param_per_species, plot_param_per_species_ms(
    param_per_iteration, file.in = "output/fig_param_per_species.jpg"), format = "file"), 
  
  ## - Plot the effect of traits on sensitivity (fig. 3)
  tar_target(fig_trait_effect, plot_trait_effect_ms(traits, traits_TRY, sensitivity, "output/fig_trait_effect.jpg"), 
             format = "file"), 
  
  ## Plot traits vs sensitivity to each disturbance (fig. 4)
  tar_target(fig_traits_vs_sensitivity, plot_traits_vs_sensitivity(
    traits, traits_TRY, sensitivity, "output/fig_trait_vs_sensitivity.jpg"), format = "file"), 
  
  ## - Plot the effect of climate on sensitivity (fig. 5)
  tar_target(fig_climate_effect, plot_climate_effect(
    gbif_file, sensitivity, "output/fig_climate_effect.jpg"), format = "file"), 
  
  ## - Plot the effect of climate related disturbance index on sensitivity (fig. 6)
  tar_target(fig_climate_disturbance, plot_disturbance_climate_ms(
    sensitivity, gbif_disturbance_file, "output/fig_climate_disturbance.jpg"), format = "file")
)