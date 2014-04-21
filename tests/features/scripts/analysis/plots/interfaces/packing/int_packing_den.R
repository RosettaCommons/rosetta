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
id = "int_packing-den",
author = "Jared Adolf-Bryfogle",
brief_description = "Graphs Interface metrics such as packstat and shape complementarity scores",
feature_reporter_dependencies = c("InterfaceFeatures/AntibodyFeatures"),
run=function(self, sample_sources, output_dir, output_formats){
  
  sele <- "
    SELECT
      sc_value,
      packstat,
      interface
    FROM
      interfaces
  "
  

  
  data = query_sample_sources(sample_sources, sele, char_as_factor=F)
  
  plot_field = function(p, plot_id, grid = NULL){
    
    if (! is.null(grid)){
      p <- p+ facet_grid(facets=grid)
    }
    if(nrow(sample_sources) <= 3){
      p <- p + theme(legend.position="bottom", legend.direction="horizontal")
    }
    save_plots(self, plot_id, sample_sources, output_dir, output_formats)
  }
  
  plot_parts <- list(
    geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)),
    scale_y_continuous("Feature Density"),
    theme_bw())
  
  #Basic densities of sc_value and packstat
  fields = c("sc_value", "packstat")
  for(field in fields){
    parts = list(plot_parts, scale_x_continuous("value", limit = c(0, 1.0)))
    
    group = c("sample_source")
    dens <- estimate_density_1d(data,  group, field)
    p <- ggplot(data=dens, na.rm=T) + parts +
      geom_line(aes(x, y, colour=sample_source), size=1.2) +
      ggtitle(field)
    plot_field(p, paste(field, "den_by_all", sep="_"))
    
    group = c("sample_source", "interface")
    dens <- estimate_density_1d(data,  group, field)
    p <- ggplot(data=dens, na.rm=T) + parts +
      geom_line(aes(x, y, colour=sample_source), size=1.2) +
      ggtitle(field)
    plot_field(p, paste(field, "den_by_interface", sep="_"), grid=interface ~ .)
  }
  
})) # end FeaturesAnalysis