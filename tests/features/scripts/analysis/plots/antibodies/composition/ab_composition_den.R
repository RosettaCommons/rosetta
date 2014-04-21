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
id = "ab_composition_den",
author = "Jared Adolf-Bryfogle",
brief_description = "Graphs basic antibody composition densities",
feature_reporter_dependencies = c("AntibodyFeatures"),
run=function(self, sample_sources, output_dir, output_formats){
  
  #First we run on all the interfaces in the database

  

  sele = "
  SELECT
    length,
    CDR,
    aromatic_nres/length as aromatic_makeup
  FROM
    cdr_metrics"
  
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
  
  #Length density
  parts = list(plot_parts, xlab("length"), xlim(5, 20))
  field = c("length")
  group = c("sample_source", "CDR")
  
  dens <- estimate_density_1d(data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    scale_x_continuous("", breaks = seq(min(data$length), max(data$length), 4))
    ggtitle("CDR Lengths")
  plot_field_wrap(p, "cdr_lengths_den",grid= ~ CDR)
  
  perc <- ddply(data, .(sample_source, CDR, length), function(d2){
    perc = nrow(d2)/nrow(data[data$sample_source == d2$sample_source[1] & data$CDR == d2$CDR[1],])
    data.frame(perc = perc)
  })
  
  #Length histogram
  p <- ggplot(data=perc, aes(x=length)) + 
    geom_bar(position="dodge", stat = 'identity', aes(y= perc , fill=sample_source ))+ 
    theme_bw() +
    ggtitle("CDR Lengths") +
    ylab("% of Sample Source")
    #scale_x_continuous("CDR Length", seq(min(perc$length), max(perc$length), 4))
    #scale_x_continuous("restype") +
    #scale_y_continuous("n")
  plot_field_wrap(p, "cdr_lengths_hist", grid= ~ CDR)
  
  avgs <- ddply(data, .(sample_source, CDR), function(d2){
    data.frame(m = mean(d2$length))
  })
  
  p <- ggplot(data=avgs, aes(x=CDR)) + 
    geom_bar(position="dodge", stat = 'identity', aes(y= round(m) , fill=sample_source ))+ 
    theme_bw() +
    ggtitle("Average CDR Lengths ") +
    scale_y_continuous("Avg CDR Length (rounded)", breaks = seq(0, round(max(avgs$m)), 4)) +
    xlab("CDR")
  plot_field(p, "avg_cdr_lengths_hist")
  
  #Aromaticity
  parts = list(plot_parts, xlab("% Aromatic Makeup"))
  field = c("aromatic_makeup")
  dens <- estimate_density_1d(data, group, field)
  p<- ggplot(data = dens, na.rm = T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Aromatic Composition")
  plot_field_wrap(p, "cdr_aromatic_den", grid = ~ CDR)
  
  avgs <- ddply(data, .(sample_source, CDR), function(d2){
    data.frame(m = mean(d2$aromatic_makeup))
  })
  
  p <- ggplot(data=avgs, aes(x=CDR)) + 
    geom_bar(position="dodge", stat = 'identity', aes(y= m , fill=sample_source ))+ 
    theme_bw() +
    ggtitle("Average Aromatic Composition") +
    ylab("avg composition") +
    xlab("CDR")
    #scale_x_continuous("restype") +
    #scale_y_continuous("n")
  plot_field(p, "avg_cdr_aromatic_hist")
  
  avgs <- ddply(data, .(sample_source, CDR, length), function(d2){
    data.frame(m = mean(d2$aromatic_makeup))
  })
  
  p <- ggplot(data=avgs, aes(x=length)) + 
    geom_bar(position="dodge", stat = 'identity', aes(y= m , fill=sample_source ))+ 
    theme_bw() +
    ggtitle("Average Aromatic Composition") +
    ylab("avg composition") +
    scale_x_continuous("CDR Length", breaks = seq(min(avgs$length), max(avgs$length), 4))
    #scale_x_continuous("restype") +
    #scale_y_continuous("n")
  plot_field_wrap(p, "avg_cdr_aromatic_hist_by_length", ~ CDR)
  
  sele = "
  SELECT
    cdr_residues
  FROM
    ab_metrics"
  
  data = query_sample_sources(sample_sources, sele, char_as_factor=F)
  
  #Total CDR residues density
  
  parts = list(plot_parts, xlab("n"))
  field = c("cdr_residues")
  group = c("sample_source")
  
  dens <- estimate_density_1d(data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Total CDR Residues")
  plot_field(p, "total_cdr_res_den")
  
  #Total CDR residues histogram
  perc <- ddply(data, .(sample_source, cdr_residues), function(d2){
    perc = nrow(d2)/nrow(data[data$sample_source == d2$sample_source[1],])
    data.frame(perc = perc)
  })
  p <- ggplot(data=perc, aes(x=cdr_residues)) + 
    geom_bar(position="dodge", stat = 'identity', aes(y= perc , fill=sample_source)) +
    theme_bw() +
    ggtitle("Total CDR Residues") +
    ylab("% of Sample Source") +
    scale_x_continuous("", breaks = seq(min(perc$cdr_residues), max(perc$cdr_residues), 4))
    
    #scale_x_continuous("restype") +
    #scale_y_continuous("n")
  plot_field(p, "total_cdr_res_hist")
})) # end FeaturesAnalysis