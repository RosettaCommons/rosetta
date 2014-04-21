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
id = "ab_composition-length_correlations",
author = "Jared Adolf-Bryfogle",
brief_description = "Graphs basic antibody composition densities",
feature_reporter_dependencies = c("AntibodyFeatures"),
run=function(self, sample_sources, output_dir, output_formats){
  
  #First we run on all the interfaces in the database

  

  sele = "
  SELECT
    length,
    CDR,
    struct_id
  FROM
    cdr_metrics"
  
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
  
  plot_field_wrap = function(p, plot_id, grid, columns = 3) {
    p <- p + facet_wrap(grid, ncol=columns)
    if(nrow(sample_sources) <= 3){
      p <- p + theme(legend.position="bottom", legend.direction="horizontal")
    }
    save_plots(self, plot_id, sample_sources, output_dir, output_formats)
  }
  
  #This does not assume 1 antibody per struct_id, which is currently how it works.  
  #It does assumes all CDRs are present, which is also how AntibodyInfo currently works.
  print("Calculating length pairs..")
  cdrs = c("L1", "L2", "L3", "H1", "H2", "H3")
  #cdrs = c("L1", "L3", "H3") 
  
  data <- ddply(data, .(struct_id, sample_source), function(d2){
    pairs = data.frame()
    for (outer_cdr in cdrs){
      for (inner_cdr in cdrs){
        if (inner_cdr == outer_cdr) {next}
        
        c = data.frame(CDR1 =d2[d2$CDR==outer_cdr,]$CDR, CDR2 = d2[d2$CDR==inner_cdr,]$CDR, CDR1_length = d2[d2$CDR==outer_cdr,]$length, CDR2_length = d2[d2$CDR==inner_cdr,]$length)
        pairs = rbind(pairs, c)
      }
    }
    pairs
  })
  
  print(head(data))
  
  
  perc <- ddply(data, .(sample_source, CDR1, CDR2, CDR1_length, CDR2_length), function(d2){
    total = nrow(data[
      data$sample_source==d2$sample_source &
      data$CDR1 == d2$CDR1 &
      data$CDR2 == d2$CDR2, 
      ])
    perc = nrow(d2)/total * 100
    data.frame(percent_sample_source = perc)
  })
  
  print(head(perc))
  
  heat_map_parts = list(
    geom_tile(aes(x=CDR1_length, y = CDR2_length, fill = percent_sample_source), stat='identity'),
    scale_fill_gradient())
  
  #This will be a bit difficult to see, but its meant to get a rough idea of what is going on.  Add to specific lengths to get more. 
  for (cdr in cdrs){
    
    d = perc[perc$CDR1==cdr,]
    d_all = data[data$CDR1==cdr,]
    
    p <- ggplot(data=d) +
      heat_map_parts +
      scale_x_continuous(paste("CDR", cdr, "length"), breaks = seq(min(d$CDR1_length), max(d$CDR1_length), 4)) +
      scale_y_continuous(paste("CDR", "length"), breaks = seq(min(d$CDR2_length), max(d$CDR2_length), 4))
    plot_field(p, paste("cdr_length_all_heat_map", cdr, "x_axis", sep="_"), sample_source ~ CDR2)
    
    #Points
    p <- ggplot(d_all, aes(x=CDR1_length, y = CDR2_length)) +
      #geom_point(size=1.5, position="jitter", aes(colour = percent_sample_source))+
      
      stat_smooth(data=d_all, method=lm) +
      scale_x_continuous(paste("CDR", cdr, "length"), breaks = seq(min(d$CDR1_length), max(d$CDR1_length), 4)) +
      scale_y_continuous(paste("CDR","length"), breaks = seq(min(d$CDR2_length), max(d$CDR2_length), 4))
    
    p_den = p +
      geom_point(data = d, size=1.2, position="jitter", aes(colour=percent_sample_source))
    
    p_all = p +
      geom_point(data = d_all, size=1.2, position="jitter")
    
    plot_field(p_den, paste("cdr_length_all_regression_jitter_den_points", cdr, "x_axis",sep="_"), sample_source ~ CDR2)
    plot_field(p_all, paste("cdr_length_all_regression_jitter_all_points", cdr, "x_axis", sep="_"), sample_source ~ CDR2)
  
  }
  
  
  
  #Specific length graphs.  Could do all of them (36), but that is too overwhelming.  Only writing these for ones we are interested in
  pairs = list(c("L3", "L1"), c("L3", "H3"), c("H3", "H1"), c("L2", "L3"), c("L2", "H3"), c("H2", "H3"))
  #pairs = list(c("L1", "L3"))
  
  for (pair in pairs){
    d = perc[perc$CDR1==pair[1] & perc$CDR2 == pair[2],]
    p <- ggplot(data = d) +
      heat_map_parts +
      scale_x_continuous(paste("CDR", pair[1], "length"), breaks = seq(min(d$CDR1_length), max(d$CDR1_length), 4)) +
      scale_y_continuous(paste("CDR", pair[2], "length"), breaks = seq(min(d$CDR2_length), max(d$CDR2_length), 4))
    plot_field(p, paste("cdr_length_ind_heat_map", pair[1], pair[2], sep="_"), ~ sample_source)
    
    #Points
    d_all = data[data$CDR1==pair[1] & data$CDR2 == pair[2],]
    p <- ggplot(d_all, aes(x=CDR1_length, y = CDR2_length)) +
      #geom_point(size=1.5, position="jitter", aes(colour = percent_sample_source))+
      
      stat_smooth(data=d_all, method=lm) +
      scale_x_continuous(paste("CDR", pair[1], "length"), breaks = seq(min(d$CDR1_length), max(d$CDR1_length), 4)) +
      scale_y_continuous(paste("CDR", pair[2], "length"), breaks = seq(min(d$CDR2_length), max(d$CDR2_length), 4))
    
    p_den = p +
      geom_point(data = d, size=1.2, position="jitter", aes(colour=percent_sample_source))
    
    p_all = p +
      geom_point(data = d_all, size=1.2, position="jitter")
    
    plot_field(p_den, paste("cdr_length_ind_regression_den_points_jitter", pair[1], pair[2], sep="_"), ~ sample_source)
    plot_field(p_all, paste("cdr_length_ind_regression_all_points_jitter", pair[1], pair[2], sep="_"), ~ sample_source)
  }
  
  
  
})) # end FeaturesAnalysis