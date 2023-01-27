#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions_plot.R  
#' @description R script containing all functions relative to data
#               visualisation
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Function to get the path of a file, and create directories if they don't exist
#' @param file.in character: path of the file, filename included (ex: "plot/plot.png")
create_dir_if_needed <- function(file.in){

  path.in <- strsplit(file.in, "/")[[1]]
  if(length(path.in) > 1){
    for(i in 1:(length(path.in)-1)){
      if(i == 1) path.in_i <- path.in[i]
      else path.in_i <- paste(path.in_i, path.in[i], sep = "/")
      if(!dir.exists(path.in_i)) dir.create(path.in_i)
    }
  }
}


#' Figure of parameter value per species for the manuscript
#' @param param_per_iteration value of each parameter at each iteration per species
#' @param dbh.ref reference dbh (in mm) to use to calculate sensitivity
#' @param logratio.ref reference logratio between tree dbh (in mm) and plot quadratic diameter to use to calculate sensitivity
#' @param I.ref reference disturbance intensity (from 0 to 1) to use to calculate sensitivity
#' @param time.ref reference time between two measurements (in years) to use to calculate sensitivity
#' @param file.in Name and path of the plot to save
plot_param_per_species_ms <- function(param_per_iteration, dbh.ref = 250, I.ref = 0.75,
                                      logratio.ref = 0, time.ref = 5, file.in){
  
  # Create directory if needed
  create_dir_if_needed(file.in)
  
  # Initialize output
  out <- list()
  
  # Identify disturbances 
  disturbances.in <- unique(param_per_iteration$disturbance)
  
  # Create a vector of colors for plotting
  color.vector <- (data.frame(disturbance = disturbances.in) %>%
                     left_join(data.frame(disturbance = c("biotic", "fire", "other", "snow", "storm"), 
                                          color = c("#90A955", "#F77F00", "#5F0F40", "#006D77", "#4361EE")), 
                               by = "disturbance"))$color
  
  # Create a vector of colors for the plot legend
  color.legend <- (data.frame(disturbance = c("storm", "fire", "other", "biotic", "snow")) %>%
                     left_join(data.frame(disturbance = c("biotic", "fire", "other", "snow", "storm"), 
                                          color = c("#90A955", "#F77F00", "#5F0F40", "#006D77", "#4361EE")), 
                               by = "disturbance"))$color
  
  # Create a legend for the plot
  plot.legend <- cowplot::get_legend(
    data.frame(x = c(1:5), y = c(1:5), ymin = c(0:4), ymax = c(2:6), 
               disturbance = factor(disturbances.in, 
                                    levels = c("storm", "fire", "other", "biotic", "snow"))) %>%
      ggplot(aes(x = x, y = y, color = disturbance)) + 
      geom_point() + 
      geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0) + 
      scale_color_manual(values = color.legend) + 
      theme(legend.title = element_blank(), 
            legend.key = element_blank(), 
            legend.text = element_text(size = 19))
  )
  
  # Loop on the three factor levels
  for(j in 1:3){
    
    # Level of factor j
    factor.j <- c("species sensitivity", "dbh effect", "dominance effect")[j]
    
    # Initialize the list that will contain the plots of plot i
    plot.list.j <- list()
    
    # Initialize number of species per disturbance
    n.sp.per.dist <- c()
    
    # Loop on all disturbances to extract data
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
        mutate(pd = 1 - (1 - plogis(a0 + a1*logratio.scaled + b*I*dbh.scaled^c))^time.ref) %>%
        # Average per species
        group_by(species) %>%
        summarise(sensitivity_mean = mean(pd), 
                  sensitivity_ql = as.numeric(quantile(pd, probs = 0.025)), 
                  sensitivity_qh = as.numeric(quantile(pd, probs = 0.975)), 
                  dbh.effect_mean = mean(c), 
                  dbh.effect_ql = as.numeric(quantile(c, probs = 0.025)), 
                  dbh.effect_qh = as.numeric(quantile(c, probs = 0.975)), 
                  dominance.effect_mean = mean(a1, na.rm = TRUE), 
                  dominance.effect_ql = as.numeric(quantile(a1, probs = 0.025, na.rm = TRUE)), 
                  dominance.effect_qh = as.numeric(quantile(a1, probs = 0.975, na.rm = TRUE))) %>%
        # Arrange species by sensitivity value
        ungroup() %>%
        filter(!(species %in% c("Other conifer", "Other broadleaf"))) %>%
        mutate(species = factor(species, levels = .$species[order(.$sensitivity_mean)])) %>%
        # Format for plotting
        gather(key = "variable", value = "value", colnames(.)[which(colnames(.) != "species")]) %>%
        separate(col = "variable", into = c("parameter", "value.type"), sep = "_") %>%
        spread(key = "value.type", value = "value") %>%
        mutate(parameter = gsub("\\.", "\\ ", parameter)) %>%
        mutate(parameter = ifelse(parameter != "sensitivity", parameter, "species sensitivity")) %>%
        mutate(parameter = factor(parameter, levels = c("species sensitivity", "dbh effect", "dominance effect")))
      
      # Case where the plot should be empty
      if(factor.j == "dominance effect" & disturbances.in[i] %in% c("other", "fire", "biotic")){
        
        # Make an empty plot, or use the legend if graph at the top right
        if(j == 3 & disturbances.in[i] == "other") plot.ij <- plot.legend
        else plot.ij <- ggplot() + theme_void() 
        
      }else{
        
        # Build plot ij
        plot.ij <- data.i %>%
          filter(!is.na(species)) %>%
          mutate(parameter = as.character(parameter)) %>%
          filter(parameter == factor.j) %>%
          ggplot(aes(x = species, y = mean, ymin = ql, ymax = qh)) +
          geom_errorbar(width = 0, color = color.vector[i]) +
          geom_point(color = color.vector[i]) + 
          coord_flip() + 
          facet_wrap(~ parameter, scales = "free_x") + 
          xlab("") + ylab("") +
          theme(panel.background = element_rect(color = "black", fill = "white"), 
                panel.grid = element_blank(), 
                strip.background = element_blank(), 
                axis.text.y = element_text(size = 9, face = "italic"), 
                axis.text.x = element_text(size = 7), 
                strip.text = element_text(size = 16), 
                plot.margin = unit(c(0, 0, 0, 0), "cm"))
        
        # If dbh or dominance effect, remove axis ticks and add a vertical line for 0
        if(j > 1) plot.ij <- plot.ij + 
            theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
            geom_hline(yintercept = 0, linetype = "dashed") +
            ylim(-2.3, 2.3)
        
        # If sensitivity, scale between 0 and 1
        if(j == 1) plot.ij <- plot.ij + ylim(0, 1)
        
        # If not storm disturbance (positioned on top), remove strip title
        if(disturbances.in[i] != "storm") plot.ij <- plot.ij + theme(strip.text = element_blank())
        
        # If not snow disturbance (positioned at bottom), remove x-axis text and ticks
        if(disturbances.in[i] != "snow") plot.ij <- plot.ij + 
            theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
        
      }
      
      # Add to the list
      eval(parse(text = paste0("plot.list.j$", disturbances.in[i], " <- plot.ij")))
      
      # Add also the number of species
      n.sp.per.dist <- c(n.sp.per.dist, length(unique(data.i$species)))
      
    }
    
    # Assemble plots ij, add label only if first column
    if(j == 1){
      plot.j <- plot_grid(plotlist = plot.list.j[c("storm", "fire", "other", "biotic", "snow")], 
                          ncol = 1, align = "v", labels = c("(a)", "", "", "(b)", ""),
                          rel_heights = (data.frame(disturbance = c("storm", "fire", "other", "biotic", "snow")) %>%
                                           left_join(data.frame(disturbance = disturbances.in, h = (n.sp.per.dist + 10)), 
                                                     by = "disturbance"))$h)
    }else{
      plot.j <- plot_grid(plotlist = plot.list.j[c("storm", "fire", "other", "biotic", "snow")], 
                          ncol = 1, align = "v", 
                          rel_heights = (data.frame(disturbance = c("storm", "fire", "other", "biotic", "snow")) %>%
                                           left_join(data.frame(disturbance = disturbances.in, h = (n.sp.per.dist + 10)), 
                                                     by = "disturbance"))$h)
    }
    
    # Add to the output list
    eval(parse(text = paste0("out$", gsub("\\ ", "\\.", factor.j), " <- plot.j")))
    
  }
  
  
  # Assemble all plots
  plot.out <- plot_grid(plotlist = out, nrow = 1, align = "h", rel_widths = c(1.5, 1, 1))
  
  ## - save the plot
  ggsave(file.in, plot.out, width = 27, height = 30, units = "cm", dpi = 600, bg = "white")
  return(file.in)
}


#' Plot the effect of traits on disturbance sensitivity
#' @param traits dataframe containing trait values per species
#' @param traits_TRY dataframe containing trait values from TRY per species
#' @param sensivity dataframe containing the sensitivity to each disturbance
#' @param file.in Name of the file to save
plot_trait_effect_ms <- function(traits, traits_TRY, sensitivity, file.in){
  
  # Create dir if needed
  create_dir_if_needed(file.in)
  
  # rename disturbance sensitivity 
  disturbance_sensitivity.in <- sensitivity
  
  # Identify the disturbances
  disturbances.in <- names(disturbance_sensitivity.in)
  
  # Rearrange traits table
  traits.in <- traits %>%
    left_join(traits_TRY, by = "species") %>%
    dplyr::select(
      "species", 
      "Wood dens." = "wood.density_g.cm3", 
      "Shade tol." = "shade.tolerance", 
      "Root mass frac." = "Root_mass_fraction", 
      "Bark thick." = "bark.thickness_mm", 
      "H/dbh ratio" = "height.dbh.ratio", 
      "Lifespan" = "TRY_plant.lifespan_year", 
      "Max. growth" = "growth.max", 
      "Leaf C/N" = "TRY_leaf.CN.ratio_g.cm3", 
      "Leaf Nmass" = "TRY_leaf.N.mass_mg.g", 
      "Leaf thick." = "TRY_leaf.thickness_mm", 
      "p50"
    )
  
  # Center and scale the trait values
  traits.in <- scale_data(traits.in, var = colnames(traits.in)[c(2:dim(traits.in)[2])])
  
  # Loop on all traits
  for(i in 1:(dim(traits.in)[2] - 1)){
    # Identify the name of trait i
    trait.i <- colnames(traits.in)[i+1]
    # Create a table with only species and trait i
    traits.i <- traits.in %>% dplyr::select("species", "trait" = trait.i)
    # Loop on all type of disturbances
    for(j in 1:length(disturbances.in)){
      
      # Create a table with trait i and sensitivity to disturbance j
      data.ij <- disturbance_sensitivity.in[[j]] %>%
        mutate(p.logit = log(p/(1 - p))) %>%
        group_by(species) %>%
        summarize(w = 1/var(p.logit), 
                  p.logit = mean(p.logit)) %>%
        mutate(p = plogis(p.logit)) %>%
        left_join((traits.i), 
                  by = "species") %>%
        drop_na()
      
      # Only perform a test if there is enough data
      if(length(unique(data.ij$species)) > 3){
        
        # Fit a model depending on model type chosen
        model.ij <- lm(p.logit ~ trait, weights = w, data = data.ij)
        
        # Extract results
        table.ij <- data.frame(
          trait = trait.i, 
          disturbance = disturbances.in[j],
          n = dim(data.ij)[1],
          Est.sup = confint.lm(model.ij)[2, 2], 
          Est.inf = confint.lm(model.ij)[2, 1],
          Est = summary(model.ij)$coefficients[2, 1]
        )
        
      }else{table.ij <- data.frame(trait = trait.i, disturbance = disturbances.in[j], n = NA_real_, 
                                   Est.sup = NA_real_, Est.inf = NA_real_, Est = NA_real_)}
      
      
      # Add to the list containing the final results
      if(i == 1 & j == 1) data <- table.ij
      else data <- rbind.data.frame(data, table.ij)
    }
  }
  
  # Add categories for each trait
  data <- data %>%
    mutate(trait.category = case_when(
      trait %in% c("Wood dens.", "Lifespan", "Max. growth") ~ "Growth vs.\n  survival", 
      trait %in% c("Leaf thick.", "Stomata cond.", "p50") ~ "Drought \n traits", 
      trait %in% c("Leaf C/N", "Leaf Nmass") ~ "Growth vs.\n defense", 
      trait %in% c("Shade tol.") ~ "Shade \ntolerance", 
      TRUE ~ "Architectural\ntraits"
    ), 
    significance = ifelse((Est.inf > 0 | Est.sup < 0), "*", ""), 
    label = ifelse(n > 3, paste0("(", n, ") ", significance), "")) %>%
    filter(disturbance %in% c("storm", "fire", "biotic")) %>%
    mutate(disturbance = factor(disturbance, levels = c("storm", "fire", "biotic")))
  
  ## - Make the plot
  plot.out <- data %>%
    group_by(disturbance) %>%
    mutate(label.pos = max(Est.sup, na.rm = TRUE)) %>%
    mutate(label.pos = label.pos + 0.5*(max(Est.sup, na.rm = T) - min(Est.inf, na.rm = T))) %>%
    ggplot(aes(x = trait, y = Est, fill = disturbance)) + 
    geom_hline(yintercept = 0, color = "grey", linetype = "dashed", size = 0.3) + 
    geom_errorbar(aes(ymin = Est.inf, ymax = Est.sup, color = significance), width = 0) + 
    geom_point(size = 1.5, aes(color = significance), shape = 21) +
    scale_color_manual(values = c("#6C757D", "black")) +
    facet_grid(trait.category ~ disturbance, scales = "free", space = "free_y") +
    geom_text(aes(label = label, y = label.pos, color = significance), 
              size = 3, hjust = "inward", show.legend = F) +
    scale_fill_manual(values = c("#4361EE", "#F77F00", "#90A955")) +
    xlab("") + ylab("Trait effect on disturbance sensitivity") +
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 13, angle = 360),
          axis.title = element_text(size = 15),
          legend.position = "none",
          axis.text.x = element_text(size = 8, angle = 360),
          axis.text.y = element_text(size = 14, angle = 360)) + 
    coord_flip() 
  
  ## - Save the plot
  ggsave(file.in, plot.out, width = 20, height = 12, units = "cm", dpi = 600, bg = "white")
  return(file.in)
  
}



#' Plot senstivity to all disturbances against traits on one multipanel
#' @param traits dataset containing trait values per species
#' @param traits_TRY dataset containing trait values per species from TRY
#' @param sensitivity dataset containing disturbance sensitivity per species for each mcmc iteration
#' @param file.in Name of the file to save (including path)
plot_traits_vs_sensitivity <- function(traits, traits_TRY, sensitivity, file.in){
  
  ## - Create directory if needed
  create_dir_if_needed(file.in)
  
  ## - Initialize plot list
  plots.out <- list()
  
  ## - Assemble the two disturbance sensitivity files
  disturbance_sensitivity.in <- sensitivity
  
  ## - Loop on all disturbance types to assemble data sets
  for(i in 1:length(names(disturbance_sensitivity.in))){
    data.in.i <- as.data.frame(disturbance_sensitivity.in[[i]]) %>%
      mutate(disturbance = names(disturbance_sensitivity.in)[i])
    if(i == 1) data.in <- data.in.i
    if(i > 1) data.in <- rbind(data.in, data.in.i)
  }
  
  ## - Rearrange traits table
  traits.in <- traits %>%
    left_join(traits_TRY, by = "species") %>%
    dplyr::select(
      "species", 
      "Wood density" = "wood.density_g.cm3", 
      "Shade tolerance" = "shade.tolerance", 
      "Root mass fraction" = "Root_mass_fraction", 
      "Bark thickness" = "bark.thickness_mm", 
      "H to dbh ratio" = "height.dbh.ratio", 
      "Leaf CN ratio" = "TRY_leaf.CN.ratio_g.cm3", 
      "Leaf NP ratio" = "TRY_leaf.NP.ratio_g.cm3", 
      "Leaf Nmass" = "TRY_leaf.N.mass_mg.g", 
      "Leaf Pmass" = "TRY_leaf.P.mass_mg.g", 
      "SLA" = "TRY_leaf.sla_mm2mg-1", 
      "Leaf thickness" = "TRY_leaf.thickness_mm", 
      "Lifespan" = "TRY_plant.lifespan_year", 
      "Stomata conductance" = "TRY_stomata.conductance_millimolm-2s-1", 
      "Maximum growth" = "growth.max"
    )
  
  
  # Loop on all traits
  for(trait.i in colnames(traits.in)[which(colnames(traits.in) != "species")]){
    
    # Format data for the model
    data.i <- data.in %>%
      # Logit transformation and weight
      mutate(p.logit = log(p/(1 - p))) %>%
      group_by(species, disturbance) %>%
      summarize(w = 1/var(p.logit), 
                p = mean(p),
                p_025.logit = quantile(p.logit, probs = 0.025), 
                p_975.logit = quantile(p.logit, probs = 0.975), 
                p.logit = mean(p.logit)) %>%
      mutate(p_025 = plogis(p_025.logit), 
             p_975 = plogis(p_975.logit)) %>%
      left_join((traits.in %>% dplyr::select("species", "trait" = trait.i)), 
                by = "species") %>%
      # remove NA
      drop_na()
    
    # Scale weight so that the sum equals the number of observations
    data.i$w <- data.i$w*dim(data.i)[1]/sum(data.i$w)
    
    # Fit a model
    model.i <- lmer(p.logit ~ trait + (1|species), weights = w, data = data.i)
    
    # Data with the predictions of the model
    data.fit <- data.frame(
      trait = c(round(min(data.i$trait)*100, digits = 0):round(max(data.i$trait)*100, digits = 0))/100) %>%
      # Extract model coefficients
      mutate(a = summary(model.i)$coefficients[1,1], 
             a.se = summary(model.i)$coefficients[1,2], 
             b = summary(model.i)$coefficients[2,1], 
             b.se = summary(model.i)$coefficients[2,2]) %>%
      # Predict output
      mutate(fit.logit = a + b*trait, 
             fit.logit.inf = a - a.se + (b - b.se)*trait,
             fit.logit.sup = a + a.se + (b + b.se)*trait,
             fit = plogis(fit.logit), 
             fit.inf = plogis(fit.logit.inf), 
             fit.sup = plogis(fit.logit.sup), 
             p = NA_real_)
    
    # Plot predictions 
    plot.i <- data.i %>%
      mutate(disturbance = factor(disturbance, levels = c("biotic", "fire", "other", "snow", "storm"))) %>%
      ggplot(aes(x = trait, y = p)) + 
      geom_errorbar(aes(ymin = p_025, ymax = p_975, color =disturbance), width = 0, alpha = 0.5) +
      geom_point(size = 1, aes(color = disturbance), alpha = 0.5) + 
      scale_color_manual(values = c("#90A955", "#F77F00", "#5F0F40", "#006D77", "#4361EE")) +
      geom_line(data = data.fit, aes(y = fit, group = 1), inherit.aes = TRUE) + 
      geom_line(data = data.fit, aes(y = fit.inf, group = 1), linetype = "dashed", inherit.aes = TRUE) + 
      geom_line(data = data.fit, aes(y = fit.sup, group = 1), linetype = "dashed", inherit.aes = TRUE) + 
      ylab("Disturbance \n sensitivity") + xlab(trait.i) +
      theme(panel.background = element_rect(color = "black", fill = "white"), 
            panel.grid = element_blank(), 
            legend.position = "none", 
            plot.title = element_text(size = 11), 
            axis.title = element_text(size = 13)) + 
      ggtitle(paste0("Chisq = ", round(Anova(model.i)[1, 1], digits = 2), ", ",
                     scales::pvalue(Anova(model.i)[1, 3], add_p = TRUE, accuracy = 0.01))) 
    
    
    # If the model is significant, add to the output list
    if(round(Anova(model.i)[1, 3], digits = 2) <= 0.05){
      eval(parse(text = paste0("plots.out$", gsub("\\ ", "", trait.i), " <- plot.i")))
    }
  }
  
  # Create a legend for the plot
  plot.legend <- cowplot::get_legend(
    data.frame(x = c(1:5), y = c(1:5), ymin = c(0:4), ymax = c(2:6), 
               disturbance = factor(c("storm", "fire", "other", "biotic", "snow"), 
                                    levels = c("storm", "fire", "other", "biotic", "snow"))) %>%
      ggplot(aes(x = x, y = y, color = disturbance)) + 
      geom_point() + 
      geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0) + 
      scale_color_manual(values = c("#4361EE", "#F77F00", "#5F0F40", "#90A955", "#006D77")) + 
      theme(legend.title = element_blank(), 
            legend.key = element_blank(), 
            legend.text = element_text(size = 19)))
  
  # Final plot
  plot.out <- plot_grid(
    plot_grid(plotlist = plots.out, align = "hv", nrow = 1, scale = 0.9), 
    plot.legend, nrow = 1, rel_widths = c(1, 0.25)
  )
  
  # Save the plot
  if(length(names(plots.out)) == 1) ggsave(file.in, plot.out, width = 13, height = 8, units = "cm", dpi = 600, bg = "white")
  if(length(names(plots.out)) == 2) ggsave(file.in, plot.out, width = 21, height = 8, units = "cm", dpi = 600, bg = "white")
  if(length(names(plots.out)) > 3) ggsave(file.in, plot.out, width = 35, height = 6, units = "cm", dpi = 600, bg = "white")
  if(length(names(plots.out)) == 3){
    plots.out$legend <- plot.legend
    plot.out <- plot_grid(plotlist = plots.out, align = "hv", nrow = 2, scale = c(0.95, 0.95, 0.95, 0.5))
    ggsave(file.in, plot.out, width = 17, height = 13, units = "cm", dpi = 600, bg = "white")
  } 
  
  return(file.in)
  
}



#' Function to plot results of the climate analysis for the manuscript
#' @param gbif_file Name of the file containing climate data per species
#' @param sensitivity list of dataset containing disturbance sensitivity per species
#' @param file.in Where to save the plot
plot_climate_effect <- function(gbif_file, sensitivity, file.in){
  
  # Create directory if needed
  create_dir_if_needed(file.in)
  
  
  # - Make PCA 
  pca <- prcomp((fread(gbif_file) %>%
                   dplyr::select(species, mat, tmin, map) %>%
                   drop_na() %>% 
                   dplyr::select(-species)), 
                center = T, scale = T)
  # - Extract the coordinates of the individuals on pca axis
  data.climate <- data.frame(species = (fread(gbif_file) %>%
                                          dplyr::select(species, mat, tmin, map) %>%
                                          drop_na())$species, 
                             pca1 = get_pca_ind(pca)[[1]][, 1], 
                             pca2 = get_pca_ind(pca)[[1]][, 2]) 
  
  # - Extract the coordinates of the individuals on pca axis
  res.ind <- data.frame(species = (fread(gbif_file) %>%
                                     dplyr::select(species, mat, tmin, map) %>%
                                     drop_na())$species, 
                        pca1 = get_pca_ind(pca)[[1]][, 1], 
                        pca2 = get_pca_ind(pca)[[1]][, 2]) %>%
    mutate(sp = paste(substr(gsub("\\ .+", "", species), 1, 2), 
                      substr(gsub(".+\\ ", "", species), 1, 2), 
                      sep = "."))
  
  # - Extract the coordinates of the variables on pca axis
  res.var <- data.frame(var = rownames(get_pca_var(pca)[[1]]), 
                        pca1 = get_pca_var(pca)[[1]][, 1], 
                        pca2 = get_pca_var(pca)[[1]][, 2])
  
  # - Minimum and maximum in each pca axis
  pca.xmin <- -max(abs(res.ind$pca1))
  pca.xmax <- max(abs(res.ind$pca1))
  pca.ymin <- -max(abs(res.ind$pca2))
  pca.ymax <- max(abs(res.ind$pca2))
  
  # Make the plot
  plot.pca <- res.ind %>%
    ggplot(aes(x = pca1, y = pca2)) + 
    geom_point(fill = "#023E8A", color = "black", shape = 21) +
    geom_text(aes(label = sp), nudge_y = 0.1, color = "#023E8A", size = 3) +
    geom_segment(data = (res.var %>% mutate(pca1 = pca1*1.5, pca2 = pca2*1.5)), 
                 aes(x = 0, xend = pca1, y = 0, yend = pca2), 
                 arrow = arrow(length = unit(0.1, "cm")), 
                 type = "closed", color = "#D90429") + 
    geom_text(data = (res.var %>% mutate(pca1 = pca1*1.5, pca2 = pca2*1.5)), 
              aes(label = var), color = "#D90429", size = 5, 
              nudge_x = ifelse(res.var$pca1 < 0, pca.xmin/12, pca.xmax/12)) +
    geom_hline(size = 0.2, yintercept = 0, color = "#6C757D", linetype = "dashed") + 
    geom_vline(size = 0.2, xintercept = 0, color = "#6C757D", linetype = "dashed") + 
    xlim((pca.xmin-0.2), (pca.xmax+0.2)) + 
    ylim((pca.ymin-0.2), (pca.ymax+0.2)) +
    xlab(paste0("PCA1 (", round(summary(pca)$importance[2, 1]*100, digits = 2), "%)")) +
    ylab(paste0("PCA2 (", round(summary(pca)$importance[2, 2]*100, digits = 2), "%)")) +
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          panel.grid = element_blank(), 
          axis.title = element_text(size = 15))
  
  
  ## - Assemble the two disturbance sensitivity files
  disturbance_sensitivity.in <- sensitivity
  
  ## - Names of the disturbances
  disturbances.in <- names(disturbance_sensitivity.in)
  
  # Create a vector of colors for plotting
  color.in <- (data.frame(disturbance = disturbances.in) %>%
                 left_join(data.frame(disturbance = c("biotic", "fire", "other", "snow", "storm"), 
                                      color = c("#90A955", "#F77F00", "#5F0F40", "#006D77", "#4361EE")), 
                           by = "disturbance"))$color
  
  # Compile trait dataset
  data.in <- data.climate %>% rename("var.1" = colnames(.)[2], "var.2" = colnames(.)[3])
  
  
  ## - Initialize output
  plots.out <- list()
  
  ## - Loop on all type of disturbances
  for(i in 1:length(disturbances.in)){
    
    # Data to fit the model for disturbance i
    data.i <- disturbance_sensitivity.in[[i]]  %>%
      mutate(p.logit = log(p/(1 - p))) %>%
      group_by(species) %>%
      summarize(w = 1/var(p.logit), 
                p.logit = mean(p.logit)) %>%
      left_join(data.in, by = "species") %>%
      drop_na()
    
    # Fit model depending on weighting method
    model.i = lm(p.logit ~ scale(var.1) + scale(var.2), weights = w, data = data.i)
    
    # Build result table to plot results
    results.i <- data.frame(disturbance = disturbances.in[i],
                            n = length(unique(data.i$species)),
                            var = colnames(data.climate)[c(2, 3)])
    if(dim(data.i)[1] > 3) results.i <- results.i %>% 
      mutate(est.low = as.numeric(confint(model.i)[c(2, 3), 1]), 
             est.high = as.numeric(confint(model.i)[c(2, 3), 2]))
    else results.i <- results.i %>% mutate(est.low = NA_real_, est.high = NA_real_)
    
    
    # Different extraction of estimate depending on model type
    results.i <- results.i %>% 
      mutate(est = coefficients(summary(model.i))[c(2, 3), 1]) %>% 
      mutate(est = ifelse(is.na(est.low), NA_real_, est))
    
    # Add to the final table
    if(i == 1) results <- results.i
    else results <- rbind(results, results.i)
    
  }
  
  # Label for the plot
  breaks.label <- (data.frame(disturbance = c("storm", "fire", "other", "biotic", "snow")) %>%
                     left_join(results %>%
                                 mutate(label = paste0(disturbance, "\n (n = ", n, ")")) %>%
                                 dplyr::select(disturbance, label) %>%
                                 distinct(), 
                               by = "disturbance"))$label
  
  # Final plot
  plot.effect <- results %>%
    mutate(disturbance = factor(disturbance, levels = rev(c("storm", "fire", "other", "biotic", "snow"))), 
           var = toupper(var)) %>%
    ggplot(aes(x = disturbance, y = est, color = disturbance)) + 
    geom_point() + 
    geom_errorbar(aes(ymin = est.low, ymax = est.high), width = 0) + 
    scale_color_manual(values = rev(c("#4361EE", "#F77F00", "#5F0F40", "#90A955", "#006D77"))) +
    geom_hline(yintercept = 0, linetype = "dashed") + 
    ylab("Effect on disturbance sensitivity") + 
    xlab("") + 
    scale_x_discrete(labels = rev(breaks.label)) +
    facet_wrap( ~ var, nrow = 1) + 
    coord_flip() + 
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(),
          legend.position = "none", 
          axis.text.x = element_text(size = 10), 
          axis.text.y = element_text(size = 14), 
          axis.title = element_text(size = 15), 
          strip.text = element_text(size = 15))
  
  
  # Save the plot
  plot.out <- plot_grid(plot.pca, plot.effect, nrow = 1, scale = 0.9, 
                        labels = c("(a)", "(b)"), rel_widths = c(1, 1.3), align = "h")
  ggsave(file.in, plot.out, width = 25, height = 10, units = "cm", dpi = 600, bg = "white")
  
  # Return file name
  return(file.in)
}





#' Function to plot results of the climate analysis for the manuscript
#' @param gbif_file Name of the file containing species mean climate
#' @param gbif_disturbance_file Name of the file containing species mean disturbance index
#' @param sensitivity list of dataset containing disturbance sensitivity per species
#' @param file.in Where to save the plot
plot_disturbance_climate_ms <- function(sensitivity, gbif_disturbance_file, file.in){
  
  ## - Create directory if needed
  create_dir_if_needed(file.in)
  
  ## - Assemble the two disturbance sensitivity files
  disturbance_sensitivity.in <- sensitivity
  
  ## - Names of the disturbances
  disturbances.in <- names(disturbance_sensitivity.in)
  
  # Create a vector of colors for plotting
  color.in <- (data.frame(disturbance = disturbances.in) %>%
                 left_join(data.frame(disturbance = c("biotic", "fire", "other", "snow", "storm"), 
                                      color = c("#90A955", "#F77F00", "#5F0F40", "#006D77", "#4361EE")), 
                           by = "disturbance"))$color
  
  
  # Initialize plot list
  plots.disturbance <- list()
  
  # Couples disturbance - index
  disturbance.index <- data.frame(disturbance = c("snow", "storm", "fire"), 
                                  index = c("swe", "windspeed", "fwi"), 
                                  name = c("Snow Water \n Equivalent", "Mean windspeed", "Fire Weather \n Index"), 
                                  color = c("#006D77", "#4361EE", "#F77F00"))
  
  # Loop on all couples
  for(j in 1:dim(disturbance.index)[1]){
    # Data for the model
    data.model.j <- as.data.frame(disturbance_sensitivity.in[disturbance.index$disturbance[j]])
    colnames(data.model.j) <- c("species", "iter", "p")
    data.model.j <- data.model.j %>%
      left_join(fread(gbif_disturbance_file) %>% dplyr::select("species", "index" = disturbance.index$index[j]), 
                by = "species") %>%
      drop_na() %>%
      mutate(sensitivity.logit = log(p/(1 - p))) %>%
      group_by(species, index) %>%
      summarize(p_025 = quantile(p, probs = 0.025), 
                p_975 = quantile(p, probs = 0.975), 
                p = mean(p), 
                p.logit = mean(sensitivity.logit), 
                w = 1/var(sensitivity.logit))
    
    # model
    model.j <- lm(p.logit ~ index, weights = w, data = data.model.j)
    # Data with predictions
    data.fit.j <- data.frame(index = c(round(min(data.model.j$index)*100, digits = 0):
                                         round(max(data.model.j$index)*100, digits = 0)/100)) %>%
      mutate(fit.logit = predict(model.j, newdata = .), 
             fit.lwr = predict(model.j, newdata = ., interval = "confidence")[, 2],
             fit.upr = predict(model.j, newdata = ., interval = "confidence")[, 3],
             fit = plogis(fit.logit), 
             fit.inf = plogis(fit.lwr), 
             fit.sup = plogis(fit.upr), 
             p = NA_real_)
    
    
    # Plot
    plot.j <- data.model.j %>%
      ggplot(aes(x = index, y = p, group = 1)) + 
      geom_errorbar(aes(ymin = p_025, ymax = p_975), width = 0, color = "#343A40") +
      geom_point(size = 2, shape = 21, fill = disturbance.index$color[j], color = "#343A40") + 
      xlab(disturbance.index$name[j]) + 
      ylab(paste0("Sensitivity to \n", disturbance.index$disturbance[j])) +
      theme(panel.background = element_rect(color = "black", fill = "white"), 
            panel.grid = element_blank(), 
            plot.title = element_text(size = 17, face = "italic"), 
            axis.title = element_text(size = 17)) + 
      scale_y_continuous(breaks = c(0:5)*0.2) + 
      ggtitle(paste0("F = ", round(anova(model.j)[1, 4], digits = 1), ", ",
                     scales::pvalue(anova(model.j)[1, 5], add_p = TRUE, accuracy = 0.01)))
    
    # Add line and confidence interval only if the regression is significant
    if(anova(model.j)[1, 5] <= 0.05){
      plot.j <- plot.j  + 
        geom_line(data = data.fit.j, aes(y = fit), inherit.aes = TRUE, color = disturbance.index$color[j])  + 
        geom_ribbon(data = data.fit.j, aes(ymin = fit.inf, ymax = fit.sup), 
                    alpha = 0.5, fill = disturbance.index$color[j], inherit.aes = TRUE)
      
    }
    # Add to the output list
    eval(parse(text = paste0("plots.disturbance$", disturbance.index$index[j], " <- plot.j")))
  }
  
  
  # Final plot
  plot.out <- plot_grid(plotlist = plots.disturbance[c("windspeed", "fwi", "swe")], nrow = 1, align = "h", 
                        labels = paste0("(", letters[c(1:3)], ")"), scale = 0.9)
  
  
  # Save the plot
  ggsave(file.in, plot.out, width = 25, height = 8, units = "cm", dpi = 600, bg = "white")
  
  return(file.in)
  
}

