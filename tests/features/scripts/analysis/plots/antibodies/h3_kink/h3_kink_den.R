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
id = "H3_kink_den",
author = "Jared Adolf-Bryfogle",
brief_description = "Various statistics on the H3 Kink",
feature_reporter_dependencies = c("AntibodyFeatures"),
run=function(self, sample_sources, output_dir, output_formats){
  
  #First we run on all the interfaces in the database

  

  sele = "
  SELECT
    kink_type,
    RD_Hbond_dis,
    bb_Hbond_dis,
    Trp_Hbond_dis,
    qdis,
    qdih,
    anion_res - cation_res as cation_res_sep
  FROM
    ab_h3_kink_metrics"
  
  data = query_sample_sources(sample_sources, sele, char_as_factor=F)
 
  parts <- list(
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
  
  
  #Kink vs No Kink Histogram
  counts = ddply(data, .(sample_source), function(by_ss){
    result = ddply(by_ss, .(kink_type), function(by_kink) {
      
      df = data.frame(percent = length(by_kink$qdis)/length(by_ss$kink_type))
    })
  })
  
  p <- ggplot(data=counts, na.rm=T) +
    geom_bar(position="dodge", stat="identity", aes(x=kink_type, y= percent, fill=sample_source)) +
    theme_bw() +
    ggtitle("Kink Type Comparison") +
    scale_y_continuous(label = percent) +
    xlab("kink type")
  plot_field(p, "kink_type_hist")
  
  #Anion - Cation res separation
  group = c("sample_source", "kink_type")
  dens <- estimate_density_1d(data, group, c("cation_res_sep"))
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    xlab("Residue Separation") +
    ggtitle("Anion Cation Residue Separation")
  plot_field(p, "anion_cation_res_sep_den", ~kink_type)
  
  #Hbond Distance RD
  group = c("sample_source", "kink_type")
  dens <- estimate_density_1d(data, group, c("RD_Hbond_dis"))
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    xlab("Angstroms") +
    ggtitle("RD Hbond Distance")
  plot_field(p, "hbond_dis_RD_den", ~kink_type)
  
  #Hbond Distance BB
  group = c("sample_source", "kink_type")
  dens <- estimate_density_1d(data, group, c("bb_Hbond_dis"))
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    xlab("Angstroms") +
    ggtitle("BB Hbond Distance")
  plot_field(p, "hbond_dis_BB_den", ~kink_type)
  
  #Hbond Distance Trp
  group = c("sample_source", "kink_type")
  dens <- estimate_density_1d(data, group, c("Trp_Hbond_dis"))
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    xlab("Angstroms") +
    ggtitle("TRP Hbond Distance")
  plot_field(p, "hbond_dis_TRP_den", ~kink_type)
  
  #qDis
  group = c("sample_source", "kink_type")
  dens <- estimate_density_1d(data, group, c("qdis"))
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    xlab("Angstroms") +
    ggtitle("Q Distance")
  plot_field(p, "qdis_den", ~kink_type)
  
  #qdih -  This may be better off on a circular plot however ggplot2 implementation of this does not seem to exist
  group = c("sample_source", "kink_type")
  dens <- estimate_density_1d(data, group, c("qdih"))
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    xlab("Degrees") +
    ggtitle("Q Dihedral")
  plot_field(p, "qdih_den", ~kink_type)
  
})) # end FeaturesAnalysis