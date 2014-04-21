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
id = "int_ss_den",
author = "Jared Adolf-Bryfogle",
brief_description = "Graphs basic composition of the interfaces, restypes, etc",
feature_reporter_dependencies = c("InterfaceFeatures/AntibodyFeatures"),


run=function(self, sample_sources, output_dir, output_formats){

  #Aromatic Composition
  sele <- "
  SELECT
    ss_sheet_fraction,
    ss_helix_fraction,
    ss_loop_fraction,
    interface,
    side
  FROM
    interface_sides
  "
  plot_parts <- list(
    geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)),
    scale_y_continuous("Feature Density"),
    theme_bw())

  plot_field = function(p, plot_id, grid = NULL, ssLegend=T){
    
    if (! is.null(grid)){
      p <- p+ facet_grid(facets=grid)
    }
    if(ssLegend){
      if(nrow(sample_sources) <= 3){
        p <- p + theme(legend.position="bottom", legend.direction="horizontal")
      }
    }
    save_plots(self, plot_id, sample_sources, output_dir, output_formats)
  }
  
  fields = c("ss_sheet_fraction", "ss_helix_fraction", "ss_loop_fraction")
  data = query_sample_sources(sample_sources, sele, char_as_factor=F)
  
  for(field in fields){
    fieldSP = unlist(strsplit(field, split="_"))
    parts = list(plot_parts, xlab("fraction"))
    
    group = c("sample_source", "side")
    dens <- estimate_density_1d(data,  group, field)
    p <- ggplot(data=dens, na.rm=T) + parts +
      geom_line(aes(x, y, colour=sample_source), size=1.2) +
      ggtitle(paste(capwords(fieldSP), collapse=" "))
    plot_field(p, paste(field, "den_sides_by_all", sep="_"), grid=side ~ .)
    
    group = c("sample_source", "interface", "side")
    dens <- estimate_density_1d(data,  group, field)
    p <- ggplot(data=dens, na.rm=T) + parts +
      geom_line(aes(x, y, colour=sample_source), size=1.2) +
      ggtitle(paste(capwords(fieldSP), collapse=" "))
    plot_field(p, paste(field, "den_sides","by_interface", sep="_"), grid=side~interface)
  }
  

  
  
})) # end FeaturesAnalysis