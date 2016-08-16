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
id = "beta_mediated_hbonds",
author = "Matthew O'Meara",
brief_description = "",
long_description = "
example of how to condition down to specific types of interactions.
In this case, hydrogen bonds that are forming beta sheet mediated
protein-protein interfaces.",

feature_reporter_dependencies = c("StructureFeatures", "ResidueSecondaryStructureFeatures", "HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


sele <-"
SELECT
  s.tag,
  don_pdb.chain, don_pdb.resNum, don_pdb.iCode, don.atmNum,
  acc_pdb.chain, acc_pdb.resNum, acc_pdb.iCode, don.atmNum,
  hbond.energy,
	don.HBChemType, acc.HBChemType,
  geom.AHdist, geom.cosBAH, geom.cosAHD, geom.chi,
  don_ss.dssp AS don_dssp, acc_ss.dssp AS acc_dssp
FROM
  structures AS s,
  hbonds AS hbond,
  hbond_geom_coords AS geom,
  hbond_sites AS don, hbond_sites AS acc,
	hbond_sites_pdb AS don_pdb, hbond_sites_pdb AS acc_pdb,
  residue_secondary_structure AS don_ss, residue_secondary_structure AS acc_ss
WHERE
  hbond.struct_id = s.struct_id AND
	geom.struct_id = s.struct_id AND hbond.hbond_id = geom.hbond_id AND
  don.struct_id = s.struct_id AND don.site_id = hbond.don_id AND
	acc.struct_id = s.struct_id AND acc.site_id = hbond.acc_id AND
	don_pdb.site_id = don.site_id AND
  acc_pdb.site_id = acc.site_id AND
	don_ss.struct_id = s.struct_id AND don_ss.resNum = don.resNum AND
	acc_ss.struct_id = s.struct_id AND acc_ss.resNum = acc.resNum AND
	don_ss.dssp = 'E' AND acc_ss.dssp = 'E' AND
	don.HBChemType = 'hbdon_PBA' AND acc.HBChemType = 'hbacc_PBA' AND
	don_pdb.chain != acc_pdb.chain';"

f <- query_sample_sources(sample_sources, sele)

dens <- estimate_density_1d(
  data = f,
  ids = c("sample_source"),
  variable = "AHdist",
  weight_fun = radial_3d_normalization)

plot_id <- "beta_mediated_hbonds_AHdist"
ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x, y=log(y+1), colour=sample_source)) +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	ggtitle("Hydrogen Bonds A-H Distance For beta-mediated HBonds\nnormalized for equal weight per unit distance") +
	labs(x=expression(paste('Acceptor -- Proton Distance (', ring(A), ')')),
	     y="log(FeatureDensity + 1)") +
	scale_y_continuous(limits=c(0,2.9), breaks=0:2) +
	scale_x_continuous(limits=c(1.4,2.7), breaks=c(1.6, 1.9, 2.2, 2.6)) +
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

})) # end FeaturesAnalysis
