# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

check_setup()

feature_analyses <- c(feature_analyses, new("FeaturesAnalysis",
id = "int_hbonds-unsat_polars_vs",
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
    geom_point(size=1.5, pch="o"),
    #stat_smooth(color="grey"),
    stat_smooth(method=lm),
    geom_density2d(),
    #stat_density2d(aes(fill = ..level..), geom="polygon"),
    #stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE),
    theme_bw())
  
  #unsat hbonds vs hbond fraction
  p <- ggplot(data = data, aes(y=hbond_E_fraction, x=delta_unsatHbonds)) + parts +
    ggtitle("Delta Unsatisfied Polar Atoms vs Hbond Energy fraction") +
    ylab("fraction") +
    xlab("atoms")
  plot_field(p, "delta_unsat_polars_vs_hbond_E_fraction_by_all", grid=sample_source ~ .)
  plot_field(p, "delta_unsat_polars_vs_hbond_E_fraction_by_interface", grid=interface ~ sample_source)
  
  #unsat hbonds vs dG
  p <- ggplot(data = data, aes(y=dG, x=delta_unsatHbonds)) + parts +
    ggtitle("Delta unsatisfied Polar Atoms vs dG") +
    ylab("REU") +
    xlab("atoms")
  plot_field(p, "delta_unsat_polars_vs_dG_by_interface", grid=sample_source ~ .)
  plot_field(p, "delta_unsat_polars_vs_dG_by_interface", grid=interface ~ sample_source)
  
  #unsat hbonds vs dSASA
  #hbond fraction vs dunsat hbonds
  p <- ggplot(data = data, aes(y=dSASA, x=delta_unsatHbonds)) + parts +
    ggtitle("Delta unsatisfied Polar Atoms vs dSASA") +
    ylab("SASA") +
    xlab("atoms")
  plot_field(p, "delta_unsat_polars_vs_dSASA_by_all", grid=sample_source ~ .)
  plot_field(p, "delta_unsat_polars_vs_dSASA_by_interface", grid=interface ~ sample_source)
  
  #unsat hbonds vs hbond E fraction
  
  #unsat hbonds vs dG_Cross
  p <- ggplot(data = data, aes(y=dG_cross, x=delta_unsatHbonds)) + parts +
    ggtitle("Delta unsatisfied Polar Atoms vs Crossterm dG") +
    ylab("REU") +
    xlab("atoms")
  plot_field(p, "delta_unsat_polars_vs_dG_cross_by_all", grid=sample_source ~ .)
  plot_field(p, "delta_unsat_polars_vs_dG_cross_by_interface", grid=interface ~ sample_source)
  
  
})) # end FeaturesAnalysis