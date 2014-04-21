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
id = "ab_charge-clusters_den",
author = "Jared Adolf-Bryfogle",
brief_description = "Graphs basic antibody composition densities",
feature_reporter_dependencies = c("AntibodyFeatures", "CDRClusterFeatures"),
run=function(self, sample_sources, output_dir, output_formats){
  
  #First we run on all the interfaces in the database

  

  sele = "
  SELECT 
    cdr_metrics.charge as charge,
    cdr_metrics.CDR as CDR,
    cdr_metrics.length as length,
    cdr_clusters.fullcluster as cluster
  FROM
    cdr_metrics,
    cdr_clusters
   WHERE
     cdr_clusters.struct_id = cdr_metrics.struct_id AND
     cdr_clusters.CDR = cdr_metrics.CDR
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
  
  plot_field_wrap = function(p, plot_id, grid, columns = 3) {
    p <- p + facet_wrap(grid, ncol=columns)
    if(nrow(sample_sources) <= 3){
      p <- p + theme(legend.position="bottom", legend.direction="horizontal")
    }
    save_plots(self, plot_id, sample_sources, output_dir, output_formats)
  }
  
  
  #CDR Charge Histogram
  data = query_sample_sources(sample_sources, sele, char_as_factor=F)
  data$cdr_length = paste(data$CDR, data$length, sep="_")
  
  avgs = ddply(data, .(sample_source, cdr_length, cluster), function(data) {
    data.frame(m=mean(data$charge))
  })
  
  for (cdr_length in unique(avgs$cdr_length)){
    
    p <- ggplot(data=avgs[avgs$cdr_length==cdr_length,]) +
      geom_bar(position='dodge', stat='identity', aes(x=cluster, y=m, fill=sample_source)) +
      theme_bw() +
      xlab("Cluster") +
      ylab("Charge") +
      ggtitle(paste("Average CDR Charge", cdr_length))
    plot_field(p, paste("cdr_charge_hist_lengths", cdr_length, sep="_"))
    
  }
  
  
  #CDR Charge Density
  for (cluster in unique(data$cluster)) {
    
    clus_data = data[data$cluster==cluster,]
    
      dens = estimate_density_1d(clus_data, c("sample_source"), c("charge"))
      if (nrow(dens)>=1){
        p<- ggplot(data=dens, na.rm=T)+
          geom_line(aes(x, y, colour=sample_source), size = 1.2) +
          geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
          ggtitle(paste("CDR Charge", cluster)) +
          xlab("Charge") +
          scale_y_continuous("Feature Density")
       plot_field(p, paste("cdr_charge_den", cluster, sep="_"))
    }
    
    
      
    p <- ggplot(data=clus_data) +
      theme_bw() +
      ggtitle(paste("CDR Charge", cluster)) +
      scale_y_continuous(label=percent) +
      ylab("% of Sample Source")
  plot_field(p, paste("cdr_charge_hist", cluster, sep="_"))
    
  }
  
})) # end FeaturesAnalysis
  