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
id = "lys_rotamer_hbond_interaction",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("SaltBridgeFeatures", "HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
  geom.AHdist, geom.cosBAH, geom.cosAHD, geom.chi,
  lys_conf.chi4 as chi4,
  acc_site.HBChemType AS acc_chem_type,
  don_site.HBChemType AS don_chem_type
FROM
  hbond_sites AS don_site,
  protein_residue_conformation as lys_conf,
  hbonds AS hbond,
  hbond_geom_coords AS geom,
  hbond_sites AS acc_site
WHERE
  don_site.HBChemType == 'hbdon_AMO' AND
  lys_conf.struct_id  == don_site.struct_id AND lys_conf.seqpos  == don_site.resNum AND
  hbond.struct_id     == don_site.struct_id AND hbond.don_id     == don_site.site_id AND
  acc_site.struct_id  == hbond.struct_id    AND acc_site.site_id == hbond.acc_id AND
  geom.struct_id      == hbond.struct_id    AND geom.hbond_id    == hbond.hbond_id;"


f <- query_sample_sources(sample_sources, sele)
#f <- ddply(f, variables=c("sample_source", "acc_chem_type"), function(df){
#	f$counts <- nrow(df)
#	f
#})

plot_parts <- list(
  geom_hex(aes(x=chi4)),
  labs(x="CHI 4 Torsion Angle"),
  facet_grid(acc_chem_type ~ sample_source),
  theme_bw())

plot_id <- "lys_rotamer_hbond_interaction_AHdist"

ggplot(f, aes(y=AHdist)) + plot_parts +
  ggtitle("Lysine Rotamer vs HBond AHdist Length") +
  labs(y=expression(paste('Acceptor -- Proton Distance (', ring(A), ')')))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

plot_id <- "lys_rotamer_hbond_interaction_cosBAH"
ggplot(f, aes(y=acos(cosBAH)*180/pi)) + plot_parts +
  ggtitle("Lysine Rotamer vs HBond Interaction cosBAH Angle") +
  labs(y=expression(paste('Base -- Acceptor -- Hydrogen (degrees)')))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

plot_id <- "lys_rotamer_hbond_interaction_cosAHD"
ggplot(f, aes(y=acos(cosAHD)*180/pi)) + plot_parts +
  ggtitle("Lysine Rotamer vs HBond Interaction cosAHD Angle") +
  labs(y=expression(paste('Acceptor -- Hydrogen -- Donor (degrees)')))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
