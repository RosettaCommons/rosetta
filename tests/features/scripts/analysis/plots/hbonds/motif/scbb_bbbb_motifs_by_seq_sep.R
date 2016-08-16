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
id = "scbb_bbbb_motifs_by_seq_sep",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

source("scripts/analysis/plots/hbonds/hbond_geo_dim_scales.R")

#st motif
sele <-"
SELECT
  aBB2.resNum - aST1.resNum AS seq_sep,
  geom1.AHdist, geom1.cosBAH, geom1.cosAHD, geom1.chi,
  geom1.AHdist, geom1.cosBAH, geom1.cosAHD, geom1.chi,
  aST1.HBChemType AS hb1_acc_chem_type
FROM
  hbonds AS aST1_dBB1,
  hbonds AS aBB2_dBB2,
  hbond_sites AS aST1, hbond_sites AS dBB1,
  hbond_sites AS aBB2, hbond_sites AS dBB2,
  hbond_geom_coords AS geom1,
  hbond_geom_coords AS geom2
WHERE
  aST1_dBB1.struct_id = aBB2_dBB2.struct_id AND
  aST1_dBB1.struct_id = aST1.struct_id AND
  aST1_dBB1.struct_id = dBB1.struct_id AND
  aST1_dBB1.struct_id = aBB2.struct_id AND
  aST1_dBB1.struct_id = dBB2.struct_id AND
  aST1_dBB1.struct_id = geom1.struct_id AND
  aST1_dBB1.struct_id = geom2.struct_id AND
  aST1_dBB1.hbond_id = geom1.hbond_id AND
  aBB2_dBB2.hbond_id = geom2.hbond_id AND
  aST1_dBB1.acc_id = aST1.site_id AND
  aST1_dBB1.don_id = dBB1.site_id AND
  aBB2_dBB2.acc_id = aBB2.site_id AND
  aBB2_dBB2.don_id = dBB2.site_id AND
  aST1.HBChemType != 'hbacc_PBA' AND
  --aST1.HBChemType = 'hbacc_HXL'AND
  dBB1.HBChemType = 'hbdon_PBA' AND
  aBB2.HBChemType = 'hbacc_PBA' AND
  dBB2.HBChemType = 'hbdon_PBA' AND
  aST1.resNum = dBB2.resNum;
"
f <- query_sample_sources(sample_sources, sele)

f <- ddply(f, .(sample_source, hb1_acc_chem_type),
  transform, counts = length(sample_source))

plot_id <- "seq_sep_histogram"
p <- ggplot(f) + theme_bw() +
  stat_bin(aes(x=seq_sep, y=log(..count..), colour=sample_source), geom="line", binwidth=1, position="identity") +
  facet_wrap( ~ hb1_acc_chem_type) +
  geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
  ggtitle("SC-BB + BB-BB HBond Motifs") +
  scale_x_log_pos_neg("Sequence Separation", breaks=c(-200,-25,-4,0,4,25,200)) +
  scale_y_log10("Count")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
