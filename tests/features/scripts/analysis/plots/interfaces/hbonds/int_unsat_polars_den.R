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
id = "int_hbonds-unsat_polars_den",
author = "Jared Adolf-Bryfogle",
brief_description = "Graphs basic interface energy information",
feature_reporter_dependencies = c("InterfaceFeatures/AntibodyFeatures"),
run=function(self, sample_sources, output_dir, output_formats){
  
  sele <- "
  SELECT
    delta_unsatHbonds,
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
  
  #Basic Densities
  
  field = "delta_unsatHbonds"
  parts = list(plot_parts, scale_x_continuous("n"))
  
  group = c("sample_source")
  dens <- estimate_density_1d(data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Delta Unsatisfied Polar Atoms")
  plot_field(p, paste("delta_unsat_polars", "den_by_all", sep="_"), )
  
  group = c("sample_source", "interface")
  dens <- estimate_density_1d(data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Delta Unsatisfied Polar Atoms")
  plot_field(p, paste("delta_unsat_polars", "den_by_interface", sep="_"), grid=interface ~ .)
  
  #Reproduction of Ben Strange's plots within a density.
  #"A comparison of successful and failed protein interface designs highlights the challenges of designing buried hydrogen bonds"
  
  #unsat per 1000 dSASA
  data$unsat_per_thousand = (data$delta_unsatHbonds*1000)/data$dSASA
  
  group = c("sample_source")
  field = "unsat_per_thousand"
  dens <- estimate_density_1d(data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    xlab(expression(paste("No. Unsatisfied /", A^2))) +
    ggtitle("Delta Unsatisfied Polar Atoms per area")
  plot_field(p, paste("delta_unsat_polars_per_1000_dSASA", "den_by_all", sep="_"), )
  
  group = c("sample_source", "interface")
  dens <- estimate_density_1d(data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    xlab(expression(paste("No. Unsatisfied /", A^2))) +
    ggtitle("Delta Unsatisfied Polar Atoms per area")
  plot_field(p, paste("delta_unsat_polars_per_1000_dSASA", "den_by_interface", sep="_"), grid=interface ~ .)
  
  
})) # end FeaturesAnalysis