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
id = "ab_hbonds-intra_cdr_den",
author = "Jared Adolf-Bryfogle",
brief_description = "Self-CDR Hbonds Excluding those arising from  BB-BB beta sheet",
feature_reporter_dependencies = c("AntibodyFeatures", "HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){
  
  #Checked, works perfectly fine:
  
  if ("FALSE" %in% opt$options$include_cdr4 & "FALSE" %in% opt$options$cdr4_only){
  sele = "
SELECT
    DISTINCT
    hb.energy as energy,
    don.struct_id as struct_id,
    don.resNum as resnum1,
    acc.resNum as resnum2,
    hb_geom.AHdist as distance,
    don_c.CDR as CDR,
    don.atmType as don_atm,
    acc.atmType as acc_atm,
    don.HBChemType as don_type,
    acc.HBChemType as acc_type
  FROM
    hbond_sites AS don,
    hbond_sites AS acc,
    hbonds AS hb,
    hbond_geom_coords as hb_geom,
    cdr_residues as  don_c,
    cdr_residues as  acc_c,
    residue_secondary_structure as don_ss,
    residue_secondary_structure as acc_ss
  WHERE
    acc.struct_id == don.struct_id AND
    don.struct_id == hb.struct_id AND
    hb_geom.struct_id = hb.struct_id AND
    don_c.struct_id = hb.struct_id AND
    acc_c.struct_id = hb.struct_id AND
    don_ss.struct_id = hb.struct_id AND
    acc_ss.struct_id = hb.struct_id AND
    don.site_id =  hb.don_id AND
    acc.site_id = hb.acc_id  AND
    hb_geom.hbond_id = hb.hbond_id AND
    don_c.resNum = don.resNum AND
    acc_c.resNum = acc.resNum AND
    don_c.CDR == acc_c.CDR AND
    acc_ss.resNum = acc.resNum AND
    don_ss.resNum = don.resNum AND
    NOT (don_ss.dssp = acc_ss.dssp AND don_ss.dssp = 'E' AND acc_ss.dssp='E' AND 
    don.HBChemType = 'hbdon_PBA' and acc.HBChemType='hbacc_PBA') AND
    don_c.CDR NOT LIKE '%Proto%' AND
    acc_c.CDR NOT LIKE '%Proto%'
  "
  }
  
  if ("TRUE" %in% opt$options$include_cdr4){
    sele = "
    SELECT
    DISTINCT
    hb.energy as energy,
    don.struct_id as struct_id,
    don.resNum as resnum1,
    acc.resNum as resnum2,
    hb_geom.AHdist as distance,
    don_c.CDR as CDR,
    don.atmType as don_atm,
    acc.atmType as acc_atm,
    don.HBChemType as don_type,
    acc.HBChemType as acc_type
    FROM
    hbond_sites AS don,
    hbond_sites AS acc,
    hbonds AS hb,
    hbond_geom_coords as hb_geom,
    cdr_residues as  don_c,
    cdr_residues as  acc_c,
    residue_secondary_structure as don_ss,
    residue_secondary_structure as acc_ss
    WHERE
    acc.struct_id == don.struct_id AND
    don.struct_id == hb.struct_id AND
    hb_geom.struct_id = hb.struct_id AND
    don_c.struct_id = hb.struct_id AND
    acc_c.struct_id = hb.struct_id AND
    don_ss.struct_id = hb.struct_id AND
    acc_ss.struct_id = hb.struct_id AND
    don.site_id =  hb.don_id AND
    acc.site_id = hb.acc_id  AND
    hb_geom.hbond_id = hb.hbond_id AND
    don_c.resNum = don.resNum AND
    acc_c.resNum = acc.resNum AND
    don_c.CDR == acc_c.CDR AND
    acc_ss.resNum = acc.resNum AND
    don_ss.resNum = don.resNum AND
    NOT (don_ss.dssp = acc_ss.dssp AND don_ss.dssp = 'E' AND acc_ss.dssp='E' AND 
    don.HBChemType = 'hbdon_PBA' and acc.HBChemType='hbacc_PBA')
    "
  }
  
  if ("TRUE" %in% opt$options$cdr4_only){
    sele = "
    SELECT
    DISTINCT
    hb.energy as energy,
    don.struct_id as struct_id,
    don.resNum as resnum1,
    acc.resNum as resnum2,
    hb_geom.AHdist as distance,
    don_c.CDR as CDR,
    don.atmType as don_atm,
    acc.atmType as acc_atm,
    don.HBChemType as don_type,
    acc.HBChemType as acc_type
    FROM
    hbond_sites AS don,
    hbond_sites AS acc,
    hbonds AS hb,
    hbond_geom_coords as hb_geom,
    cdr_residues as  don_c,
    cdr_residues as  acc_c,
    residue_secondary_structure as don_ss,
    residue_secondary_structure as acc_ss
    WHERE
    acc.struct_id == don.struct_id AND
    don.struct_id == hb.struct_id AND
    hb_geom.struct_id = hb.struct_id AND
    don_c.struct_id = hb.struct_id AND
    acc_c.struct_id = hb.struct_id AND
    don_ss.struct_id = hb.struct_id AND
    acc_ss.struct_id = hb.struct_id AND
    don.site_id =  hb.don_id AND
    acc.site_id = hb.acc_id  AND
    hb_geom.hbond_id = hb.hbond_id AND
    don_c.resNum = don.resNum AND
    acc_c.resNum = acc.resNum AND
    don_c.CDR == acc_c.CDR AND
    acc_ss.resNum = acc.resNum AND
    don_ss.resNum = don.resNum AND
    NOT (don_ss.dssp = acc_ss.dssp AND don_ss.dssp = 'E' AND acc_ss.dssp='E' AND 
    don.HBChemType = 'hbdon_PBA' and acc.HBChemType='hbacc_PBA') AND 
    don_c.CDR LIKE '%Proto%' AND
    acc_c.CDR LIKE '%Proto%'
    "
  }
  
  
  data = query_sample_sources(sample_sources, sele, char_as_factor=F)
  
  #print(sum(data$struct_id==1))
  #print(sum(data$struct_id==2))
  #print(sum(data$struct_id==3))
  plot_parts <- list(
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
  
  #Hbond Energy density
  field = "energy"
  group = c("sample_source", "CDR")
  dens <- estimate_density_1d(data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + plot_parts +
    xlab("REU") +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Intra-CDR H-Bond Energies")
  plot_field(p, paste("intra_cdr_hbonds", field, "tog_den_by_cdr", sep="_"), ~ CDR)
  
  #Hbond Distances
  field = "distance"
  group = c("sample_source", "CDR")
  dens <- estimate_density_1d(data,  group, field)
  p <- ggplot(data=dens, na.rm=T) + plot_parts +
    xlab("Angstroms") +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Intra-CDR H-Bond Distances")
  plot_field(p, paste("intra_cdr_hbonds", field, "tog_den_by_cdr", sep="_"), ~ CDR)
  
  #Hbond Counts per CDR
  counts <- ddply(data, .(sample_source, struct_id, CDR), function(data){
    n = length(data$energy > 0)
    df = data.frame(n=n)
  })
  
  field = "n"
  group = c("sample_source", "CDR")
  
  dens <- estimate_density_1d(counts,  group, field)
  p <- ggplot(data=dens, na.rm=T) + plot_parts +
    xlab("hbonds") +
    geom_line(aes(x, y, colour=sample_source), size=1.2) +
    ggtitle("Average Intra-CDR H-Bonds")
  plot_field(p, "intra_cdr_hbond_counts_tog_den_by_cdr", ~ CDR)
    
  
  avgs <- ddply(counts, .(sample_source, CDR), function(data){
    data.frame(sample_source=data$sample_source, CDR=data$CDR, m=mean(data$n))
  })
  
  #Histograms
  p <- ggplot(data=avgs, na.rm=T) + 
    geom_bar(aes(x=CDR, y=m, fill=sample_source), position="dodge", stat='identity') +
    xlab("hbonds") +
    ylab("n")
    ggtitle("Average Intra-CDR H-Bonds")
  plot_field(p, "intra_cdr_hbond_counts_tog_hist_by_cdr")
  
})) # end FeaturesAnalysis