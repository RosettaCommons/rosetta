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
id = "ab_anchor_dis",
author = "Jared Adolf-Bryfogle",
brief_description = "Graphs basic antibody composition densities",
feature_reporter_dependencies = c("AntibodyFeatures"),
run=function(self, sample_sources, output_dir, output_formats){
  
  #First we run on all the interfaces in the database
  
  if ("FALSE" %in% opt$options$include_cdr4 & "FALSE" %in% opt$options$cdr4_only){
    sele = "
  SELECT
    anchor_CN_distance,
    CDR,
    length
  FROM
    cdr_metrics where CDR NOT LIKE '%Proto%'"
  }
  
  if ("TRUE" %in% opt$options$include_cdr4){
    sele = "
  SELECT
    anchor_CN_distance,
    CDR,
    length
  FROM
    cdr_metrics"
  } 
  
  if ("TRUE" %in% opt$options$cdr4_only){
    sele = "
  SELECT
    anchor_CN_distance,
    CDR,
    length
  FROM
    cdr_metrics where CDR LIKE '%Proto%'"
  }
  
  data = query_sample_sources(sample_sources, sele, char_as_factor=F)
 
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
  
  plot_field_wrap = function(p, plot_id, grid, columns = 3) {
    p <- p + facet_wrap(grid, ncol=columns)
    if(nrow(sample_sources) <= 3){
      p <- p + theme(legend.position="bottom", legend.direction="horizontal")
    }
    save_plots(self, plot_id, sample_sources, output_dir, output_formats)
  }
  
  #C-N distance by CDR and length density
  parts = list(plot_parts, xlab("Angstrom"))
  field = c("anchor_CN_distance")
  group = c("sample_source", "CDR")
  
  dens <- estimate_density_1d(data[data$anchor_CN_distance < 15, ],  group, field)
  group = c("sample_source", "CDR", "length")

  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("C-N distance")
  plot_field_wrap(p, "c_n_dis_by_cdr_den", ~ CDR)
  
  
    
    
})) # end FeaturesAnalysis