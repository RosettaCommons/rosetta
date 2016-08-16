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
id = "cosBAH_chem_type_bfac30",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
  geom.cosBAH,
  acc_site.HBChemType AS acc_chem_type,
  don_site.HBChemType AS don_chem_type
FROM
  hbond_geom_coords AS geom,
  hbonds AS hbond,
  hbond_sites AS don_site,
  hbond_sites AS acc_site,
  hbond_sites_pdb AS don_site_pdb,
  hbond_sites_pdb AS acc_site_pdb
WHERE
  hbond.struct_id = geom.struct_id AND
  hbond.hbond_id =  geom.hbond_id AND
  hbond.struct_id = don_site.struct_id AND
  hbond.don_id = don_site.site_id AND
  don_site_pdb.struct_id = hbond.struct_id AND don_site_pdb.site_id = hbond.don_id AND
  acc_site_pdb.struct_id = hbond.struct_id AND acc_site_pdb.site_id = hbond.acc_id AND
  don_site_pdb.heavy_atom_temperature < 30 AND
  acc_site_pdb.heavy_atom_temperature < 30 AND
  hbond.struct_id = acc_site.struct_id AND
	abs( don_site_pdb.resNum - acc_site_pdb.resNum ) > 5 AND
  hbond.acc_id = acc_site.site_id;"

all_geom <- query_sample_sources(sample_sources, sele)

dens <- estimate_density_1d(
  data = all_geom,
  ids = c("sample_source", "acc_chem_type", "don_chem_type"),
  variable = "cosBAH")

# The 'hbdon_' part of the donor labels doesn't fit so strip them out
dens$don_chem_type <- sub("^hbdon_", '', dens$don_chem_type)

plot_id = "cosBAH_chem_type_bfac30"
ggplot(data=dens) +
	geom_line(aes(x=acos(x)*180/pi, y=y, colour=sample_source)) +
	geom_indicator(aes(colour=sample_source, indicator=counts, group=sample_source)) +
	facet_grid(don_chem_type ~ acc_chem_type) +
	ggtitle("Hydrogen Bonds BAH Angle by Chemical Type, SeqSep > 5; BFactors < 30\n(normalized for equal volume per unit distance)") +
	scale_x_continuous(paste('Base -- Acceptor -- Hydrogen (degrees)')) +
	scale_y_continuous("FeatureDensity", limits=c(0,6)) +
	theme_bw()
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
