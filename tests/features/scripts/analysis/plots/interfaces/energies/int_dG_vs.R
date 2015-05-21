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
id = "int_energies-dG_vs",
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
  data_rm_out = data[data$dG<=5000 & data$dG>-5000,]#Remove high energy outliers
  
  parts = list(
    geom_point(size=1.2, pch="o"),
    #stat_smooth(color="grey"),
    stat_smooth(method=lm),
    geom_density2d(size=.5),
    #stat_density2d(aes(fill = ..level..), geom="polygon"),
    #stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE),
    theme_bw())
  
  #dG vs dSASA
  p <- ggplot(data=data_rm_out, aes(y = dSASA, x = dG, colour=sample_source)) + parts +
    ggtitle("dG vs dSASA") +
    ylab("SASA") +
    xlab("REU")
  plot_field(p, "dG_vs_dSASA_by_all")
  plot_field(p, "dG_vs_dSASA_by_interface", grid=~ interface)
  
  #dG vs dG_cross
  p <- ggplot(data = data_rm_out, aes(x=dG, y=dG_cross, colour=sample_source)) + parts +
    ggtitle("dG vs Crossterm dG") +
    xlab("REU") +
    ylab("REU")
  plot_field(p, "dG_vs_dG_cross_by_all")
  plot_field(p, "dG_vs_dG_cross_by_interface", grid= ~ interface)
  
  p <- ggplot(data = data_rm_out, aes(x=dG, y=dG_cross, color=dSASA)) + parts +
    ggtitle("dG vs Crossterm dG") +
    xlab("REU") +
    ylab("REU") +
    scale_fill_hue(l=40)
  plot_field(p, "dG_vs_dG_cross_col_by_dSASA_by_all", grid=sample_source ~ .)
  plot_field(p, "dG_vs_dG_cross_col_by_dSASA_by_interface", grid=sample_source ~ interface) 
  
  #dG_cross vs dSASA
  p <- ggplot(data = data_rm_out, aes(x=dG_cross, y=dSASA, colour=sample_source)) + parts +
    ggtitle("dG_cross vs dSASA") +
    xlab("REU") +
    ylab("SASA")
  plot_field(p, "dg_cross_vs_dSASA" )
  plot_field(p, "dG_cross_vs_dSASA_by_interface", grid= ~ interface)
  
  
})) # end FeaturesAnalysis