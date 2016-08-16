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
id = "BAH_chem_type",
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
	hbond_sites AS acc_site
WHERE
	hbond.struct_id = geom.struct_id AND
	hbond.hbond_id = geom.hbond_id AND
	hbond.struct_id = don_site.struct_id AND
	hbond.don_id = don_site.site_id AND
	hbond.struct_id = acc_site.struct_id AND
	hbond.acc_id = acc_site.site_id;"

f <- query_sample_sources(sample_sources, sele)

dens <- estimate_density_1d(
	data = f,
	ids = c("sample_source", "acc_chem_type", "don_chem_type"),
	variable = "cosBAH")

# Order the plots better and give more descriptive labels
f$don_chem_type <- factor(f$don_chem_type,
	levels = c("hbdon_IMD", "hbdon_IME", "hbdon_GDE", "hbdon_GDH",
		"hbdon_AHX", "hbdon_HXL", "hbdon_IND", "hbdon_AMO", "hbdon_CXA", "hbdon_PBA"),
	labels = c("dIMD: h", "dIME: h", "dGDE: r", "dGDH: r",
		"dAHX: y", "dHXL: s,t", "dIND: w", "dAMO: k", "dCXA: n,q", "dPBA: bb"))

# Order the plots better and give more descriptive labels
f$acc_chem_type <- factor(f$acc_chem_type,
	levels = c("hbacc_IMD", "hbacc_IME", "hbacc_AHX", "hbacc_HXL",
		"hbacc_CXA", "hbacc_CXL", "hbacc_PBA"),
	labels = c("aIMD: h", "aIME: h", "aAHX: y", "aHXL: s,t",
		"aCXA: n,q", "aCXL: d,e", "aPBA: bb"))

plot_id = "hbond_BAH_chem_type"
ggplot(data=dens) +
	geom_line(aes(x=acos(x)*180/pi, y=log(y+1), colour=sample_source)) +
	geom_indicator(aes(colour=sample_source, indicator=counts, group=sample_source)) +
	facet_grid(don_chem_type ~ acc_chem_type) +
	ggtitle("HBonds BAH Angle by Chemical Type, B-Factor < 30\n(normalized for equal volume per unit distance)") +
	scale_x_continuous(paste('Base -- Acceptor -- Hydrogen (degrees)')) +
	scale_y_continuous("log(FeatureDensity + 1)")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
