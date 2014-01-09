# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

check_setup()

feature_analyses <- c(feature_analyses, new("FeaturesAnalysis",
id = "interface_overall_metrics",
author = "Jared Adolf-Bryfogle",
brief_description = "Graphs Interface metrics such as packstat and shape complementarity scores",
feature_reporter_dependencies = c("InterfaceFeatures"),
run=function(self, sample_sources, output_dir, output_formats){
  
  sele <- "
    SELECT
      sc_value,
      packstat,
      dSASA,
      dG,
      dG_cross,
      delta_unsatHbonds,
      interface
    FROM
      interfaces
  "
  
  plot_parts <- list(
    geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)),
    scale_y_continuous("Feature Density"),
    theme_bw())
  
  plot_field = function(p, plot_id, grid = NULL){
    
    if (! is.null(grid)){
      p <- p+ facet_grid(facets=grid)
    }
    if(nrow(sample_sources) <= 3){
      p <- p + theme(legend.position="bottom", legend.direction="horizontal")
    }
    save_plots(self, plot_id, sample_sources, output_dir, output_formats)
  }
  
  data = query_sample_sources(sample_sources, sele)
  
  #Basic densities of sc_value and packstat
  fields = c("sc_value", "packstat")
  for(field in fields){
    parts = list(plot_parts, scale_x_continuous("value", limit = c(0, 1.0)))
    
    group = c("sample_source")
    dens <- estimate_density_1d(data,  group, field)
    p <- ggplot(data=dens, na.rm=T) + parts +
      geom_line(aes(x, y, colour=sample_source), size=1.2) +
      ggtitle(field)
    plot_field(p, paste(field, "den", sep="_"))
    
    group = c("sample_source", "interface")
    dens <- estimate_density_1d(data,  group, field)
    p <- ggplot(data=dens, na.rm=T) + parts +
      geom_line(aes(x, y, colour=sample_source), size=1.2) +
      ggtitle(field)
    plot_field(p, paste(field, "den_by_interface", sep="_"), grid=interface ~ .)
  }
  
  #Scatterplots
  #sc_value vs packstat
  parts = list(
    geom_point(size=1.5, pch="o"),
    #stat_smooth(color="grey"),
    stat_smooth(method=lm),
    geom_density2d(),
    #stat_density2d(aes(fill = ..level..), geom="polygon"),
    #stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE),
    theme_bw())
  
  p <- ggplot(data=data, aes(x =sc_value, y=packstat)) + parts + 
    ggtitle("sc_value vs packstat") +
    scale_x_continuous("sc_value", limit = c(0, 1.0)) +
    scale_y_continuous("packstat", limit = c(0, 1.0))
  plot_field(p, "sc_value_vs_packstat", grid = sample_source ~ .)
  plot_field(p, "sc_value_vs_packstat_by_interface", grid=interface ~ sample_source)
  
  #sc_value vs dSASA
  p <- ggplot(data=data, aes(x =sc_value, y=dSASA)) + parts + 
    ggtitle("sc_value vs dSASA") +
    scale_x_continuous("sc_value", limit = c(0, 1.0)) +
    scale_y_continuous("Buried SASA")
  plot_field(p, "sc_value_vs dSASA", grid = sample_source ~ .)
  plot_field(p, "sc_value_vs_dSASA_by_interface", grid=interface ~ sample_source)
  
  #packstat vs dSASA
  p <- ggplot(data = data, aes(x=packstat, y=dSASA)) + parts +
    ggtitle("packstat vs dSASA") +
    scale_x_continuous("packstat", limit=c(0,1.0)) +
    scale_y_continuous("Buried SASA")
  plot_field(p, "packstat_vs_dSASA", grid = sample_source ~.)
  plot_field(p, "packstat_vs_dSASA_by_interface", grid = interface ~ sample_source)
  
  #sc_value vs dG
  p <- ggplot(data = data[data$dG<5000,], aes(x=sc_value, y=dG)) + parts +
    ggtitle("sc_value_vs_dG") +
    scale_x_continuous("sc_value", limit = c(0, 1.0)) +
    scale_y_continuous("REU")
  plot_field(p, "sc_value_vs_dG", grid=sample_source ~ .)
  plot_field(p, "sc_value_vs_dG_by_interface", grid=interface ~ sample_source)
  
  #sc_value vs crossterm
  p <- ggplot(data = data[data$dG<5000,], aes(x=sc_value, y=dG_cross)) + parts +
    ggtitle("sc_value vs dG_cross") +
    scale_x_continuous("sc_value", limit = c(0, 1.0)) +
    scale_y_continuous("REU")
  plot_field(p, "sc_value_vs_dG_cross", grid=sample_source ~ .)
  plot_field(p, "sc_value_vs_dG_cross_by_interface", grid=interface ~ sample_source)
  
  #deltaUnsatHbonds vs sc_value
  p <- ggplot(data = data, aes(x=delta_unsatHbonds, y=sc_value)) + parts +
    ggtitle("Unsatisfied Hbonds at the interface vs sc_value") +
    scale_x_continuous("n") +
    scale_y_continuous("sc_value", limit = c(0, 1.0))
  plot_field(p, "dUnsatHbonds_vs_sc_value", grid=sample_source ~ .)
  plot_field(p, "dUnsatHbonds_vs_sc_value_by_interface", grid=interface ~ sample_source)
  
  #deltaUnsatHbonds vs packstat
  p <- ggplot(data = data, aes(x = delta_unsatHbonds, y=packstat)) + parts +
    ggtitle("Unsatisfied Hbonds at the interface vs packstat") +
    scale_x_continuous("n") +
    scale_y_continuous("packstat", limit = c(0, 1.0))
  plot_field(p, "dUnsatHbonds_vs_packstat", grid=sample_source ~ .)
  plot_field(p, "dUnsatHbonds_vs_packstat_by_interface", grid = interface ~ sample_source)
  #3D Plots
  
  #sc_value vs dG vs dSASA
  
  #Sides:
  
  #sc_value vs interface_energy
  
})) # end FeaturesAnalysis