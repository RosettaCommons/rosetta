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
id = "int_energies-dG_cross_vs",
author = "Jared Adolf-Bryfogle",
brief_description = "Graphs basic interface energy information",
feature_reporter_dependencies = c("InterfaceFeatures/AntibodyFeatures"),
run=function(self, sample_sources, output_dir, output_formats){
  
  sele <- "
  SELECT
    dG,
    dG_cross,
    delta_unsatHbonds,
    hbond_E_fraction,
    dSASA,
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
  #data_rm_out = data[data$dG<=5000 & data$dG>-5000,]#Remove high energy outliers
  
  parts = list(
    geom_point(size=1.2, pch="o"),
    #stat_smooth(color="grey"),
    stat_smooth(method=lm),
    geom_density2d(),
    #stat_density2d(aes(fill = ..level..), geom="polygon"),
    #stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE),
    theme_bw())
  
  #dG_cross vs dG
  p <- ggplot(data = data, aes(y=dG, x=dG_cross, colour=sample_source)) + parts +
    ggtitle("Crossterm dG vs dG") +
    xlab("REU") +
    ylab("REU")
  plot_field(p, "dG_cross_vs_dG_by_all", grid=sample_source ~ .)
  plot_field(p, "dG_cross_vs_dG_by_interface", grid=interface ~ sample_source)
  
  
  #dG_cross vs dSASA
  p <- ggplot(data = data, aes(x=dG_cross, y=dSASA, colour=sample_source)) + parts +
    ggtitle("Crossterm dG vs dSASA") +
    xlab("REU") +
    ylab("SASA")
  plot_field(p, "dg_cross_vs_dSASA")
  plot_field(p, "dG_cross_vs_dSASA_by_interface", grid=~ interface)
  
  
  #dG_cross vs dunsathbonds
  p <- ggplot(data = data, aes(x=dG_cross, y=delta_unsatHbonds, color=sample_source)) + parts +
    ggtitle("Crossterm dG vs Delta unsatisfied Polar Atoms") +
    xlab("REU") +
    ylab("n")
  plot_field(p, "dG_cross_vs_delta_unsat_polars")
  plot_field(p, "dG_cross_vs_delta_unsat_polars", grid=~interface)
  
  
})) # end FeaturesAnalysis