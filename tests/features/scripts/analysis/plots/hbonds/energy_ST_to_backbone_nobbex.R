# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

library(ggplot2)
library(plyr)
library("scales")
library(grid)

# check_setup()

feature_analyses <- c(feature_analyses, new("FeaturesAnalysis",
id = "energy_ST_to_backbone_nobbex",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
  hbond.energy,
  acc_site.HBChemType AS acc_chem_type,
  don_site.HBChemType AS don_chem_type,
  don_site.resType as don_res_type,
	( acc_site.resNum - don_site.resNum ) as seq_sep
FROM
  hbonds AS hbond,
  hbond_sites AS don_site,
  hbond_sites AS acc_site
WHERE
  hbond.struct_id = don_site.struct_id AND hbond.don_id = don_site.site_id AND
  hbond.struct_id = acc_site.struct_id AND hbond.acc_id = acc_site.site_id AND
  don_site.HBChemType = 'hbdon_HXL' AND
  acc_site.HBChemType = 'hbacc_PBA' AND
  abs( acc_site.resNum - don_site.resNum ) < 5
;"

f <- query_sample_sources(sample_sources, sele)

# This is deprecated please use the hbond_chem_types table for the lables instead
# Order the plots better and give more descriptive labels
f$don_chem_type <- factor(f$don_chem_type,
  levels = c("hbdon_IMD", "hbdon_IME", "hbdon_GDE", "hbdon_GDH",
    "hbdon_AHX", "hbdon_HXL", "hbdon_IND", "hbdon_AMO", "hbdon_CXA", "hbdon_PBA"),
  labels = c("dIMD: h", "dIME: h", "dGDE: r", "dGDH: r",
    "dAHX: y", "dHXL: s,t", "dIND: w", "dAMO: k", "dCXA: n,q", "dPBA: bb"))

# This is deprecated please use the hbond_chem_types table for the lables instead
# Order the plots better and give more descriptive labels
f$acc_chem_type <- factor(f$acc_chem_type,
  levels = c("hbacc_IMD", "hbacc_IME", "hbacc_AHX", "hbacc_HXL",
    "hbacc_CXA", "hbacc_CXL", "hbacc_PBA"),
  labels = c("aIMD: h", "aIME: h", "aAHX: y", "aHXL: s,t",
    "aCXA: n,q", "aCXL: d,e", "aPBA: bb"))

dens <- estimate_density_1d(
  f, c("sample_source", "don_res_type", "seq_sep"), "energy")

plot_id <- "hbond_energy_ST_to_BB"
p <- ggplot(data=dens) + theme_bw() +
  geom_line(aes(x=x, y=y, colour=1)) +
  geom_indicator(aes(indicator=counts, colour=1, group=1)) +
  facet_grid( seq_sep ~ don_res_type) +
  ggtitle("BB/SC Hydrogen Bond Energy") +
  scale_x_continuous("Rosetta Energy Units (Unweighted)", limits=c(-1.5,0), breaks=c(-1.5, -1, -.5)) +
  scale_y_continuous("Energy")
p <- p + theme(legend.position="bottom", legend.direction="horizontal")

save_plots(self, plot_id, sample_sources, output_dir, output_formats)

# I think this is slower sele <-"
# I think this is slower CREATE TEMPORARY TABLE bif_hbonds_to_bb AS SELECT
# I think this is slower   hbond.struct_id AS struct_id,
# I think this is slower   hbond.don_id AS don_id,
# I think this is slower   hbond.acc_id AS acc_id,
# I think this is slower   hbond.energy AS energy
# I think this is slower FROM
# I think this is slower   hbonds as hbond,
# I think this is slower   hbond_sites AS acc_site,
# I think this is slower   hbond_site_environment as acc_env
# I think this is slower WHERE
# I think this is slower   hbond.struct_id = acc_site.struct_id AND hbond.acc_id = acc_site.site_id AND
# I think this is slower   acc_site.HBChemType = 'hbacc_PBA' AND
# I think this is slower   acc_env.num_hbonds >= 2;
# I think this is slower
# I think this is slower SELECT
# I think this is slower   hbond.energy,
# I think this is slower   acc_site.HBChemType AS acc_chem_type,
# I think this is slower   don_site.HBChemType AS don_chem_type,
# I think this is slower   don_site.resType as don_res_type,
# I think this is slower   ( acc_site.resNum - don_site.resNum ) as seq_sep
# I think this is slower FROM
# I think this is slower   bif_hbonds_to_bb AS hbond,
# I think this is slower   bif_hbonds_to_bb AS hbond2,
# I think this is slower   hbond_sites AS don_site,
# I think this is slower   hbond_sites AS don_site2
# I think this is slower WHERE
# I think this is slower   hbond.struct_id = don_site.struct_id AND hbond.don_id = don_site.site_id AND
# I think this is slower   hbond2.struct_id = hbond.struct_id AND hbond2.acc_id = hbond.acc_id AND
# I think this is slower   hbond2.struct_id = don_site2.struct_id AND hbond2.don_id = don_site2.site_id AND
# I think this is slower   don_site.HBChemType = 'hbdon_HXL' AND
# I think this is slower   don_site2.HBChemType = 'hbdon_PBA'
# I think this is slower ;"

sele <-"
CREATE TEMPORARY TABLE bif_hbonds_to_bb AS SELECT
  hbond.struct_id AS struct_id,
  hbond.don_id AS don_id,
  hbond.acc_id AS acc_id,
  hbond.energy AS energy
FROM
  hbonds as hbond,
  hbond_sites AS acc_site,
  hbond_site_environment as acc_env
WHERE
  hbond.struct_id = acc_site.struct_id AND hbond.acc_id = acc_site.site_id AND
  acc_site.HBChemType = 'hbacc_PBA' AND
  acc_site.struct_id = acc_env.struct_id AND acc_site.site_id = acc_env.site_id AND
  acc_env.num_hbonds >= 2;

SELECT
  hbond.energy,
  acc_site.HBChemType AS acc_chem_type,
  don_site.HBChemType AS don_chem_type,
  don_site.resType as don_res_type,
  ( acc_site.resNum - don_site.resNum ) as seq_sep
FROM
  bif_hbonds_to_bb AS hbond,
  bif_hbonds_to_bb AS hbond2,
  hbond_sites AS acc_site,
  hbond_sites AS don_site,
  hbond_sites AS don_site2
WHERE
  hbond.struct_id = don_site.struct_id AND hbond.don_id = don_site.site_id AND
  hbond.struct_id = acc_site.struct_id AND hbond.acc_id = acc_site.site_id AND
  hbond2.struct_id = hbond.struct_id AND hbond2.acc_id = hbond.acc_id AND
  hbond2.struct_id = don_site2.struct_id AND hbond2.don_id = don_site2.site_id AND
  don_site.HBChemType = 'hbdon_HXL' AND
  don_site2.HBChemType = 'hbdon_PBA'
;"

f <- query_sample_sources(sample_sources, sele)

# This is deprecated please use the hbond_chem_types table for the lables instead
# Order the plots better and give more descriptive labels
f$don_chem_type <- factor(f$don_chem_type,
  levels = c("hbdon_IMD", "hbdon_IME", "hbdon_GDE", "hbdon_GDH",
    "hbdon_AHX", "hbdon_HXL", "hbdon_IND", "hbdon_AMO", "hbdon_CXA", "hbdon_PBA"),
  labels = c("dIMD: h", "dIME: h", "dGDE: r", "dGDH: r",
    "dAHX: y", "dHXL: s,t", "dIND: w", "dAMO: k", "dCXA: n,q", "dPBA: bb"))

# This is deprecated please use the hbond_chem_types table for the lables instead
# Order the plots better and give more descriptive labels
f$acc_chem_type <- factor(f$acc_chem_type,
  levels = c("hbacc_IMD", "hbacc_IME", "hbacc_AHX", "hbacc_HXL",
    "hbacc_CXA", "hbacc_CXL", "hbacc_PBA"),
  labels = c("aIMD: h", "aIME: h", "aAHX: y", "aHXL: s,t",
    "aCXA: n,q", "aCXL: d,e", "aPBA: bb"))

dens <- estimate_density_1d(
  f, c("sample_source", "don_res_type", "seq_sep"), "energy")

plot_id <- "hbond_energy_ST_to_BB_w_bbbb_hbonds"
p <- ggplot(data=dens) + theme_bw() +
  geom_line(aes(x=x, y=y, colour=1)) +
  geom_indicator(aes(indicator=counts,colour=1,group=1)) +
  facet_grid( seq_sep ~ don_res_type) +
  ggtitle("BB/SC Hydrogen Bond Energy With Competing BB/BB Hbond") +
  scale_x_continuous("Rosetta Energy Units (Unweighted)", limits=c(-1.5,0), breaks=c(-1.5, -1, -.5)) +
  scale_y_continuous("Energy")
p <- p + theme(legend.position="bottom", legend.direction="horizontal")

save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
