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
id = "ab_charge_den",
author = "Jared Adolf-Bryfogle",
brief_description = "Graphs basic antibody composition densities",
feature_reporter_dependencies = c("AntibodyFeatures"),
run=function(self, sample_sources, output_dir, output_formats){
  
  #First we run on all the interfaces in the database

  

  sele = "
  SELECT
    net_charge,
    paratope_charge
  FROM
    ab_metrics"
  
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
  
  get_charge_seq <- function(d, spacer = 1){
    d2 = d[! is.nan(d$charge) || is.na(d$charge),]
    r = seq(min(d2$charge), max(d2$charge), spacer)
    r
  }
  #Net Charge density
  parts = list(plot_parts, xlab("charge"))
  field = c("net_charge")
  group = c("sample_source")
  
  dens <- estimate_density_1d(data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    scale_x_continuous("charge", breaks = seq(min(data$net_charge), max(data$net_charge), 2)) +
    ggtitle("Antibody Net Charge")
  plot_field(p, "net_charge_den")
  
  perc <- ddply(data, .(sample_source, net_charge), function(d2){
    perc = nrow(d2)/nrow(data[data$sample_source == d2$sample_source[1],])
    data.frame(perc = perc)
  })
  
  #Net Charge histogram
  p <- ggplot(data=perc ) + 
    geom_bar(position="dodge", stat='identity', aes(x = net_charge, y= perc , fill=sample_source)) +
    theme_bw() +
    scale_x_continuous("charge", breaks = seq(min(perc$net_charge), max(perc$net_charge), 2)) +
    ggtitle("Antibody Net Charge") +
    ylab("% of Sample Source") +
    scale_y_continuous(label=percent)
    #scale_x_continuous("restype") +
    #scale_y_continuous("n")
  plot_field(p, "net_charge_hist")
  

  
  #Avg Net Charge Histogram
  avgs <- ddply(data, .(sample_source), function(d2){
    data.frame(m = mean(d2$net_charge))
  })
  p <- ggplot(data=avgs, ) + 
    geom_bar(position="dodge", stat='identity', aes(x = sample_source, y= m , fill=sample_source)) +
    theme_bw() +
    ggtitle("Average Antibody Net Charge") +
    ylab("charge") +
    xlab("Sample Source")
    #scale_x_continuous("restype") +
    #scale_y_continuous("n")
  plot_field(p, "avg_net_charge_hist")
  
  #Paratope Charge density
  parts = list(plot_parts, xlab("charge"))
  field = c("paratope_charge")
  group = c("sample_source")
  
  dens <- estimate_density_1d(data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    scale_x_continuous("charge", breaks = seq(min(data$paratope_charge), max(data$paratope_charge), 2)) +
    ggtitle("Paratope Net Charge")
  plot_field(p, "paratope_charge_den")
  
  perc <- ddply(data, .(sample_source, paratope_charge), function(d2){
    perc = nrow(d2)/nrow(data[data$sample_source == d2$sample_source[1],])
    data.frame(perc = perc)
  })
  
  #Paratope Charge histogram
  p <- ggplot(data=perc, aes(x=paratope_charge )) + 
    geom_bar(position="dodge", stat='identity', aes( y= perc, fill=sample_source)) +
    theme_bw() +
    ggtitle("Paratope Net Charge") +
    scale_x_continuous("charge", breaks = seq(min(perc$paratope_charge), max(perc$paratope_charge), 2)) +
    ylab("% of Sample Source") +
    scale_y_continuous(label=percent)
    #scale_x_continuous("restype") +
    #scale_y_continuous("n")
  plot_field(p, "paratope_charge_hist")
  
  #Avg Paratope histogram
  avgs <- ddply(data, .(sample_source), function(d2){
    data.frame(m = mean(d2$paratope_charge))
  })
  p <- ggplot(data=avgs ) + 
    geom_bar(position="dodge", stat='identity', aes(x = sample_source, y= m , fill=sample_source)) +
    theme_bw() +
    ggtitle("Average Paratope Net Charge") +
    ylab("charge")
    #scale_x_continuous("restype") +
    #scale_y_continuous("n")
  plot_field(p, "avg_paratope_charge_hist")
  
  sele = "
  SELECT 
    charge,
    CDR,
    length
  FROM
    cdr_metrics
  "
  
  #CDR Charge Density
  data = query_sample_sources(sample_sources, sele, char_as_factor=F)
  dens <- estimate_density_1d(data, c("sample_source", "CDR"), c("charge"))
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size = 1.2) +
    #scale_x_continuous("Charge", breaks = get_charge_seq(data, 1)) +
    ggtitle("CDR Charge")
  plot_field_wrap(p, "cdr_charge_den", ~ CDR)
  
  
  #CDR Charge Histogram
  avgs = ddply(data, .(sample_source, CDR), function(data) {
    data.frame(m=mean(data$charge))
  })
  perc = ddply(data, .(sample_source, CDR, charge), function(d2){
    perc = nrow(d2)/nrow(data[data$CDR == d2$CDR[1] & data$sample_source == d2$sample_source[1],]) *100
    data.frame(perc = perc)
  })
  perc$chargec = as.character(perc$charge)
  
  p <- ggplot(data=perc, aes(x = charge)) +
      geom_bar(position="dodge", stat='identity', aes(y = perc, fill=sample_source))+ 
      theme_bw() +
      #scale_x_continuous("charge",  breaks = get_charge_seq(perc, 1)) +
      ggtitle("CDR Charge") +
      scale_y_continuous("charge", label=percent) +
      ylab("% of Sample Source")
  plot_field_wrap(p, "cdr_charge_hist", ~ CDR)
    
  p <- ggplot(data=avgs, aes(x=charge)) +
    geom_bar(position="dodge", aes(x= CDR, y=m, fill=sample_source), stat = 'identity')+ 
    theme_bw() +
    xlab("CDR") +
    #scale_y_continuous("Charge", breaks = get_charge_seq(avgs, 1)) +
    ylab("charge") +
    ggtitle("Average CDR Charge")
  plot_field(p, "avg_cdr_charge_hist_by_cdr")
  
  avgs = ddply(data, .(sample_source, CDR, length), function(data) {
    data.frame(m=mean(data$charge))
  })
  avgs$lengthc = as.character(avgs$length)
  
  p <- ggplot(data=avgs, aes(x=charge)) +
    geom_bar(position="dodge", aes(x= length, y=m, fill=sample_source), stat = 'identity')+ 
    theme_bw() +
    xlab("CDR Length") +
    ylab("charge") +
    #scale_y_continuous("Charge", breaks = get_charge_seq(avgs, 1)) +
    ggtitle("Average CDR Charge")
  plot_field_wrap(p, "avg_cdr_charge_hist_by_length", ~ CDR)
  
  
})) # end FeaturesAnalysis