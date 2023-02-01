#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions_data.R  
#' @description R script containing all functions relative to data
#               importation and formatting
#' @author Julien BARRERE
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#' Center and Scale variables before fitting a model
#' @param data_model Tree data formated for the IPM. 
#' @param var Variables to scale (character vector)
scale_data <- function(data_model,  var){
  id_var <- which(colnames(data_model) %in% var) 
  out <- data_model
  unscaled = as.matrix(as.data.frame(out)[, id_var])
  scaled = scale(unscaled)
  out[, id_var] <- as.data.frame(scaled)
  colnames(out) <- colnames(data_model)
  return(out)
}

#' Compile all traits data
#' @param bark.thickness_file file containing data on bark thickness
#' @param wood.density_file file containing data on wood density
#' @param FUNDIV_tree original tree dataset (to get height diameter ratio)
#' @param FUNDIV_plot Plot-level data formatted for FUNDIV
#' @param species.in character vector of all the species for which to extract the data
compile_traits <- function(wood.density_file, shade.tolerance_file, 
                           root.depth_file, P50_file, NFI.traits_file){
  
  # Traits from NFI data
  traits.NFI <- fread(NFI.traits_file)
  
  # Species included
  species.in <- traits.NFI$species
  
  # Wood density (Chave et al. 2008 + Dryad to quote)
  wood.density <- read_xls(wood.density_file, sheet = "Data")
  colnames(wood.density) <- c("n", "family", "species", "wood.density_g.cm3", "region", "reference")
  
  # Shade tolerance (Niimenets)
  shade.tolerance <- fread(shade.tolerance_file) %>%
    mutate(shade.tolerance = as.numeric(gsub("\\,", "\\.", shade.tolerance.mean))) %>%
    dplyr::select(species, shade.tolerance)
  
  # Rooting depth (Guerrero-Ramirez et al. 2021 - Groot database)
  root.depth <- fread(root.depth_file) %>%
    mutate(species = paste(genusTNRS, speciesTNRS, sep = " ")) %>%
    filter(traitName %in% c("Rooting_depth", "Root_mass_fraction")) %>%
    filter(species %in% species.in) %>%
    dplyr::select("species", "traitName", "meanSpecies") %>%
    spread(key = "traitName", value = "meanSpecies")
  
  # Read P50 per species
  data_p50 = as.data.frame(read_xlsx(P50_file, sheet = "ALL")) %>%
    dplyr::select(species = Species.binomial, p50 = Pclose) %>%
    # Add data from Loppez et al. 2011  (https://doi.org/10.1093/aob/mct084)
    rbind(data.frame(species = c("Pinus canariensis"), 
                     p50 = mean(c(-3.77, -4.61, -3.16, -4.61, -4.47, -3.13, 
                                  -3.8, -4.32, -5.73, -6.05, -5.44, -4.76, 
                                  -5.29, -4.62, -5.74, -5.65))))
  
  # Global trait dataset
  traits <- traits.NFI %>%
    # Add wood density
    left_join((wood.density %>% 
                 group_by(species) %>%
                 summarize(wood.density_g.cm3 = mean(wood.density_g.cm3))),
              by = "species") %>%
    # Add shade tolerance
    left_join(shade.tolerance, by = "species") %>%
    # Add root depth
    left_join(root.depth, by = "species")  %>%
    # Add p50
    left_join(data_p50, by = "species") 
  
  return(traits)
}



#' Compile all traits data
#' @param TRY_file file containing TRY request
#' @param species.in character vector of all the species for which to extract the data
compile_traits_TRY <- function(TRY_file, NFI.traits_file){
  
  # -- Species name
  species.in <- fread(NFI.traits_file)$species
  
  # -- Translate TRY traits code into abbreviated traits name
  TRY.traits.name <- data.frame(
    TraitID = c(24, 3117, 146, 14, 56, 15, 46, 65, 1111, 2809, 2807, 2808, 159, 30, 318, 
                31, 719, 59, 819, 45, 773, 413, 324, 1229, 153, 865, 837, 3446), 
    trait = c("bark.thickness", "leaf.sla", "leaf.CN.ratio", "leaf.N.mass", "leaf.NP.ratio", 
              "leaf.P.mass", "leaf.thickness", "root.type", "seedbank.density", 
              "seedbank.duration", "seedbank.n.layers", "seedbank.thickness.toplayer", 
              "seedbank.type", "tolerance.drought", "tolerance.fire", "tolerance.frost", 
              "xylem.hydraulic.vulnerability", "plant.lifespan", "plant.resprouting.capacity", 
              "stomata.conductance", "crown.height", "leaf.Chl.content", "crown.length", 
              "wood.Nmass", "budbank.height.distribution", "budbank.seasonality", 
              "bark.structure", "plant.biomass")
  )
  
  
  
  ## -- Compile numeric traits data
  traits.TRY <-data.table::fread(TRY_file) %>%
    filter(!is.na(TraitID)) %>%
    filter(!is.na(StdValue)) %>%
    filter(AccSpeciesName %in% species.in) %>%
    left_join(TRY.traits.name, by = "TraitID") %>%
    mutate(trait = paste("TRY", trait, gsub("\\ ", "", UnitName), sep = "_"), 
           trait = gsub("\\/", ".", trait)) %>%
    rename("species" = "AccSpeciesName") %>%
    group_by(species, trait) %>%
    summarize(trait.value = mean(StdValue, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(trait) %>%
    mutate(n.species.per.trait = n()) %>%
    filter(n.species.per.trait >= 10) %>% 
    dplyr::select(-n.species.per.trait) %>%
    spread(key = trait, value = trait.value)  
  
  return(traits.TRY)
}




#' Function to get parameter values per iteration for each species and disturbance
#' @param rdata.file rdata contianing the model outputs
get_param_from_rdata <- function(rdata.file){
  
  # Load rdata
  load(rdata.file)
  
  # Loop on all disturbances
  for(i in 1:length(names(jags.list))){
    
    # Format the table containing parameter value per species, disturbance and iteration
    param.table.i <- ggs(as.mcmc(jags.list[[i]])) %>%
      # Add disturbance
      mutate(disturbance = names(jags.list)[i]) %>%
      # Extract information on parameter, species and country
      mutate(Param = gsub("\\[.+", "", Parameter), 
             sp = as.integer(ifelse(Param == "a0", gsub(".+\\[", "", gsub("\\,.+", "", Parameter)), 
                                    gsub(".+\\[", "", gsub("\\]", "", Parameter)))), 
             co = as.integer(ifelse(Param == "a0", gsub(".+\\,", "", gsub("\\]", "", Parameter)), NA_integer_))) %>%
      # Remove the estimation of intensity and deviance
      filter(Param != "I") %>%
      filter(Param != "deviance") %>%
      # Add name of the country and species, and weight of each species per country
      left_join(corresp.tables[[i]]$country.table, by = "co") %>%
      left_join(corresp.tables[[i]]$species.table, by = "sp") %>%
      left_join(weight.tables[[i]], by = c("species", "country")) %>%
      # No weight for the parameters that do not rely on the country
      mutate(weight = ifelse(Param == "a0", weight, 1)) %>%
      mutate(weight = ifelse(is.na(weight), 1, weight)) %>%
      # Summarize Parameter value per country (only apply to a0)
      group_by(disturbance, Iteration, Chain, Param, species) %>%
      summarize(val = sum(value*weight, na.rm = TRUE)/sum(weight, na.rm = TRUE)) %>%
      # Format to get one column per parameter
      spread(key = Param, value = val) %>%
      # Set a1 to 0 (dominance effect) if disturbance is not storm or snow
      mutate(a1 = ifelse(disturbance %in% c("storm", "snow"), a1, 0)) %>%
      # Add parameters to scale dbh and logratio
      mutate(dbh.intercept = scale.tables[[i]]$dbh.intercept, 
             dbh.slope = scale.tables[[i]]$dbh.slope, 
             logratio.intercept = scale.tables[[i]]$logratio.intercept, 
             logratio.slope = scale.tables[[i]]$logratio.slope)
    
    # Store table in the final table
    if(i == 1) param.table <- param.table.i
    else param.table <- rbind(param.table, param.table.i)
  }
  
  # return the list
  return(param.table)
  
}



#' Function to get disturbance sensitivity from parameters per iteration
#' @param param_per_iteration value of each parameter at each iteration per species
#' @param dbh.ref reference dbh (in mm) to use to calculate sensitivity
#' @param logratio.ref reference logratio between tree dbh (in mm) and plot quadratic diameter to use to calculate sensitivity
#' @param I.ref reference disturbance intensity (from 0 to 1) to use to calculate sensitivity
#' @param time.ref reference time between two measurements (in years) to use to calculate sensitivity
get_sensitivity_from_param <- function(param_per_iteration, dbh.ref = 250, I.ref = 0.75,
                                       logratio.ref = 0, time.ref = 5){
  
  # Disturbances included in the model
  disturbances.in <- unique(param_per_iteration$disturbance)
  
  # Initialize output list
  list.out <- list()
  
  # Loop on all disturbances
  for(i in 1:length(disturbances.in)){
    
    # Subset and format param_per_iteration for disturbance i
    param_i <- param_per_iteration %>%
      filter(disturbance == disturbances.in[i]) %>%
      mutate(iter = paste0("chain", Chain, "_iter", Iteration)) %>%
      ungroup() %>%
      dplyr::select(-disturbance, -Chain, -Iteration)
    
    # Calculate sensitivity
    data.i <- expand.grid(species = unique(param_i$species),
                          dbh = dbh.ref, logratio = logratio.ref, I = I.ref, 
                          iter = unique(param_i$iter)) %>%
      # Add parameters per species
      left_join(param_i, by = c("species", "iter")) %>%
      # Scale variables when needed
      mutate(dbh.scaled = dbh.intercept + dbh.slope*dbh, 
             logratio.scaled = logratio.intercept + logratio.slope*logratio) %>%
      # Compute probabilities
      mutate(p = 1 - (1 - plogis(a0 + a1*logratio.scaled + b*I*dbh.scaled^c))^time.ref) %>%
      # Keep variables of interest
      dplyr::select(species, iter, p)
    
    # Add to the output list
    eval(parse(text = paste0("list.out$", disturbances.in[i], " <- data.i")))
    
  }
  
  # Return the output list
  return(list.out)
  
}




