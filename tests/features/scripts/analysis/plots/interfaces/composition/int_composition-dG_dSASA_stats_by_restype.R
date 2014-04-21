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
id = "int_composition-dG_dSASA_stats_by_restype",
author = "Jared Adolf-Bryfogle",
brief_description = "Graphs basic composition of the interfaces, restypes, etc",
feature_reporter_dependencies = c("InterfaceFeatures/AntibodyFeatures","ResidueFeatures", "ResidueTypesFeatures", "PdbDataFeatures"),


run=function(self, sample_sources, output_dir, output_formats){

  
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
  
  #Restype composition - Overall interface - this will change to side-chain only contribution.  
  sele <-"
  SELECT
    residues.name3 as restype,
    interface_residues.interface as interface,
    interface_residues.relative_dSASA_fraction as dSASA_fraction,
    interface_residues.dSASA as dSASA,
    interface_residues.dSASA_sc as dSASA_sc,
    (dSASA - dSASA_sc) as dSASA_bb,
    interface_residues.dhSASA as dhSASA,
    interface_residues.dhSASA_sc as dhSASA_sc,
    (dhSASA - dhSASA_sc) as dhSASA_bb,
    interface_residues.dhSASA_rel_by_charge as dhSASA_rel_by_charge,
    interface_residues.dG as dG
  FROM
    residues,
    interface_residues
  WHERE
    interface_residues.struct_id == residues.struct_id and
    interface_residues.resNum == residues.resNum"
  
  res_data = query_sample_sources(sample_sources, sele, char_as_factor=F)
  
  #Scatterplot of dSASA vs dSASA fraction for each restype.
  #All interfaces - Coloring by Sample Source may be too damn confusing, so change this if need be.
  
  #dSASA vs dSASA fraction per restype
  p <- ggplot(data=res_data[res_data$dSASA > 0 & res_data$dSASA_fraction > 0,], aes(x = dSASA_fraction, y=dSASA, color=sample_source)) +
    #geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
    geom_point(size=.5, pch="o") +
    stat_smooth(method=lm) +
    geom_density2d() +
    theme_bw() +
    ggtitle("dSASA fraction vs dSASA per restype") +
    facet_wrap(~ restype, ncol=4) +
    xlab("SASA buried/SASA separated") +
    ylab("dSASA")
  if(nrow(sample_sources) <= 3){
    p <- p + theme(legend.position="bottom", legend.direction="horizontal")
  }
  save_plots(self, "dSASA_fraction_vs_dSASA_per_restype_>0dSASA", sample_sources, output_dir, output_formats)
  

#  #dSASA fractions bins per restype
#  p <- ggplot(data=res_data[res_data$dSASA > 0,], aes(x=dSASA_fraction, fill=factor(restype))) +
#    geom_bar() +
#    ggtitle("dSASA fraction bins") +
#    scale_fill_hue(l=40)+
#    theme_bw()
#  plot_field(p, "dSASA_fraction_by_restype_bins_>0dSASA_by_all", grid=sample_source ~ ., ssLegend=F)
#  plot_field(p, "dSASA_fraction_by_restype_bins_>0dSASA_by_interface", grid=sample_source ~ interface, ssLegend=F)
  
#  #dSASA bins per restype
#  p <- ggplot(data=res_data[res_data$dSASA > 0,], aes(x=dSASA, fill=factor(restype))) +
#    geom_bar() +
#    ggtitle("dSASA bins") +
#    scale_fill_hue(l=40)+
#    theme_bw()
#  plot_field(p, "dSASA_by_restype_bins_>0dSASA_by_all", grid=sample_source ~ ., ssLegend=F)
#  plot_field(p, "dSASA_by_restype_bins_>0dSASA_by_interface", grid=sample_source ~ interface, ssLegend=F)

  #dSASA density per restype
  parts = list(
    geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)),
    ylab("Feature Density"),
    xlab("SASA"),
    theme_bw())

  group = c("sample_source", "restype")
  field = "dSASA_fraction"
  dens <- estimate_density_1d(res_data[res_data$dSASA > 0,],  group, field)
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    xlab("SASA buried/SASA separated") +
    ggtitle("dSASA fraction per restype")
  if(nrow(sample_sources) <= 3){
    p <- p + theme(legend.position="bottom", legend.direction="horizontal")
  }
  p <- p + facet_wrap(~ restype, ncol=4)
  save_plots(self, "dSASA_fraction_by_restype_den_>0dSASA_by_all", sample_sources, output_dir, output_formats)
  
  
  fields = c("dSASA", "dSASA_bb", "dSASA_sc", "dhSASA", "dhSASA_bb", "dhSASA_sc", "dhSASA_rel_by_charge")
  for (field in fields) {
    group = c("sample_source", "restype")
    dens <- estimate_density_1d(res_data[res_data$dSASA > 0,],  group, field)
    p <- ggplot(data=dens, na.rm=T) + parts +
      geom_line(aes(x, y, colour=sample_source), size=1.2) +
     ggtitle( paste(field,"per restype", sep=" "))
    if(nrow(sample_sources) <= 3){
      p <- p + theme(legend.position="bottom", legend.direction="horizontal")
    }
    p <- p + facet_wrap(~ restype, ncol=4)
    save_plots(self, paste(field, "by_restype_den_>0dSASA_by_all", sep="_"), sample_sources, output_dir, output_formats)
  }
  
  #dG by restype
  group = c("sample_source", "restype")
  field = "dG"
  dens <- estimate_density_1d(res_data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("dG per restype") +
    xlab("REU")
  if(nrow(sample_sources) <= 3){
    p <- p + theme(legend.position="bottom", legend.direction="horizontal")
  }
  p <- p + facet_wrap(~ restype, ncol=4)
  save_plots(self, "dG_by_restype_den_by_all", sample_sources, output_dir, output_formats)
  
  group = c("sample_source", "restype")
  field = "dG"
  dens <- estimate_density_1d(res_data[res_data$dSASA > 0,],  group, field)
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("dG per restype") +
    xlab("REU")
  if(nrow(sample_sources) <= 3){
    p <- p + theme(legend.position="bottom", legend.direction="horizontal")
  }
  p <- p + facet_wrap(~ restype, ncol=4)
  save_plots(self, "dG_by_restype_den_>0dSASA_by_all", sample_sources, output_dir, output_formats)
  
  group = c("sample_source", "restype")
  field = "dG"
  dens <- estimate_density_1d(res_data[res_data$dSASA == 0,],  group, field)
  p <- ggplot(data=dens, na.rm=T) + parts +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("dG per restype") +
    xlab("REU")
  if(nrow(sample_sources) <= 3){
    p <- p + theme(legend.position="bottom", legend.direction="horizontal")
  }
  p <- p + facet_wrap(~ restype, ncol=4)
  save_plots(self, "dG_by_restype_den_@0dSASA_by_all", sample_sources, output_dir, output_formats)
  
  #Hard to see:
  #group = c("sample_source", "restype")
  #field = "dSASA_fraction"
  #dens <- estimate_density_1d(res_data[res_data$dSASA_fraction > 0,],  group, field)
  #p <- ggplot(data=dens, na.rm=T) + parts +
  #  geom_line(aes(x, y, color=restype)) +
  #  ggtitle("dSASA fraction per restype") +
  #  facet_grid(sample_source ~ .)
  #save_plots(self,  "dSASA_fraction_per_restype_den_combined", sample_sources, output_dir, output_formats)
  
  #group = c("sample_source", "restype")
  #field = "dSASA"
  #dens <- estimate_density_1d(res_data[res_data$dSASA > 0,],  group, field)
  #p <- ggplot(data=dens, na.rm=T) + parts +
  #  geom_line(aes(x, y, color=restype)) +
  #  ggtitle("dSASA per restype") +
  #  facet_grid(sample_source ~ .)
  #save_plots(self,  "dSASA_per_restype_den_combined", sample_sources, output_dir, output_formats)
  
  
  
  #dG vs dSASA by restype
  p <- ggplot(data=res_data[res_data$dSASA > 0,], aes(x = dG, y=dSASA, color=sample_source)) +
    #geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
    geom_point(size=.5, pch="o") +
    stat_smooth(method=lm) +
    geom_density2d() +
    theme_bw() +
    ggtitle("dG vs dSASA per restype") +
    facet_wrap(~ restype, ncol=4) +
    xlab("REU") +
    ylab("dSASA")
  if(nrow(sample_sources) <= 3){
    p <- p + theme(legend.position="bottom", legend.direction="horizontal")
  }
  save_plots(self, "dSASA_vs_dG_by_restype_>0dSASA_by_all", sample_sources, output_dir, output_formats)
  
  
  
  
  
  #group = c("sample_source", "restype", "interface")
  #field = "dSASA_fraction"
  #dens <- estimate_density_1d(res_data[res_data$dSASA_fraction > 0,],  group, field)
  #p <- ggplot(data=dens, na.rm=T, fill=restype) + parts +
  #  geom_line(aes(x, y, color=restype)) +
  #  ggtitle("Hotspot dSASA density") +
  #  facet_grid(sample_source ~ interface) +
  #save_plots(self, "dSASA_fraction_per_restype_den_combined_by_interface", sample_sources, output_dir, output_formats)
  
  
  #group = c("sample_source", "restype", "interface")
  #field = "dSASA_fraction"
  #dens <- estimate_density_1d(res_data,  group, field)
  #p <- ggplot(data=dens, na.rm=T) + parts +
    #geom_line(aes(x, y, colour=sample_source), size=1.2) +
    #ggtitle("Hotspot dSASA density of Interface residue")
  #plot_field(p, "dSASA_fraction_per_restype_by_interface", grid=restype ~ interface)
  #p <- ggplot(data=res_data[res_data$dSASA_fraction > .05,], fill=sample_source,weight=dSASA_fraction) + 
  #  geom_histogram(aes(x=restype), position="dodge")+ 
  #  theme_bw() +
  #  ggtitle("Interface ResType Composition")
  #scale_x_continuous("restype") +
  #scale_y_continuous("n")
  #plot_field(p, "restype_composition_weighted_by_dSASA_frac") 
  #plot_field(p, "restype_composition_weighted_by_dSASA_frac_by_interface", grid=interface ~ .)
  
  #restype vs avg dSASA
  #p <- ggplot(data=res_data) +
  #geom_histogram(aes(x=mean(restype)), position="dodge") +
  # theme_bw() +
  # ggtitle("Interface ResType Composition")
  #scale_x_continuous("restype") +
  #scale_y_continuous("n")
  #plot_field(p, "restype_vs_avg_dSASA_fraction") 
  #plot_field(p, "restype_vs_avg_dSASA_fraction_by_interface", grid=interface ~ .)
  
})) # end FeaturesAnalysis