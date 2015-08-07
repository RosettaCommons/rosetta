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
id = "ab_SASA-CDR_den",
author = "Jared Adolf-Bryfogle",
brief_description = "CDR Sasas",
feature_reporter_dependencies = c("AntibodyFeatures"),
run=function(self, sample_sources, output_dir, output_formats){
  
  #First we run on all the interfaces in the database

  
  if ("FALSE" %in% opt$options$include_cdr4 & "FALSE" %in% opt$options$cdr4_only){
  sele = "
  SELECT
    SASA,
    CDR,
    length
  FROM
    cdr_metrics
  WHERE
  CDR NOT LIKE '%Proto%'
    "
  }
  
  if ("TRUE" %in% opt$options$include_cdr4){
    sele = "
  SELECT
    SASA,
    CDR,
    length
  FROM
    cdr_metrics"
  }
  
  if ("TRUE" %in% opt$options$cdr4_only){
    sele = "
  SELECT
    SASA,
    CDR,
    length
  FROM
    cdr_metrics
  WHERE
  CDR LIKE '%Proto%'"
  }
  
  data = query_sample_sources(sample_sources, sele, char_as_factor=F)
 
  parts <- list(
    geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)),
    scale_y_continuous("Feature Density"),
    theme_bw())

  plot_field = function(p, plot_id, grid = NULL){
    
    if (! is.null(grid)){
      p <- p+ facet_wrap(facets=grid, ncol=3)
    }
    if(nrow(sample_sources) <= 3){
      p <- p + theme(legend.position="bottom", legend.direction="horizontal")
    }
    save_plots(self, plot_id, sample_sources, output_dir, output_formats)
  }
  
  
  #CDR SASA
  group = c("sample_source", "CDR")
  dens <- estimate_density_1d(data, group, c("SASA"))
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    xlab("SASA") +
    ggtitle("CDR SASA")
  plot_field(p, "cdr_sasa_den", ~CDR)
  
  cdr_avgs = ddply(data, .(sample_source, CDR), function(data){
    data.frame(m=mean(data$SASA))
  })
  
  len_avgs = ddply(data, .(sample_source, CDR, length), function(data){
    data.frame(m=mean(data$SASA))
  })

  p <- ggplot(data=cdr_avgs) +
    geom_bar(aes(x=CDR, y=m, fill=sample_source), position="dodge", stat='identity') +
    xlab("CDR") +
    ylab("SASA") +
    ggtitle("Average CDR SASA")
  plot_field(p, "avg_cdr_sasa_hist")
  
  p <- ggplot(data=len_avgs) +
    geom_bar(aes(x=length, y=m, fill=sample_source), position="dodge", stat='identity') +
    xlab("CDR Length") +
    ylab("SASA") +
    ggtitle("Average CDR SASA")
  plot_field(p, "avg_cdr_sasa_hist_by_length", grid= ~CDR)
  
  
  
})) # end FeaturesAnalysis