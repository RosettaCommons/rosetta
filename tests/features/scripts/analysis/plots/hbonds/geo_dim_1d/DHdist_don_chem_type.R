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
id = "DHdist_don_chem_type",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT DISTINCT
	don_atoms.base_x AS dx, don_atoms.base_y AS dy, don_atoms.base_z AS dz,
	don_atoms.atm_x  AS hx, don_atoms.atm_y  AS hy, don_atoms.atm_z  AS hz,
	don.HBChemType AS don_chem_type
FROM
	hbond_sites AS don,
	hbond_site_atoms AS don_atoms
WHERE
	don.is_donor = 1 AND
	don_atoms.struct_id = don.struct_id AND don_atoms.site_id = don.site_id;"


f <- query_sample_sources(sample_sources, sele)

# D-H distance is not stored directly in the features database,
# however it can be computed from the coordinates of the hydrogen
# bonding atoms.
f <- transform(f,
	DHdist = vector_distance(cbind(dx, dy, dz), cbind(hx, hy, hz)))

# Order the plots better and give more descriptive labels
f$don_chem_type <- factor(f$don_chem_type,
	levels = c("hbdon_IMD", "hbdon_IME", "hbdon_GDE", "hbdon_GDH",
		"hbdon_AHX", "hbdon_HXL", "hbdon_IND", "hbdon_AMO", "hbdon_CXA", "hbdon_PBA"),
	labels = c("dIMD: h", "dIME: h", "dGDE: r", "dGDH: r",
		"dAHX: y", "dHXL: s,t", "dIND: w", "dAMO: k", "dCXA: n,q", "dPBA: bb"))

f <- na.omit(f, method="r")


plot_id <- "hbond_DHdist_don_chem_type"
p <- ggplot(data=f) + theme_bw() +
	geom_point(aes(x=don_chem_type, y=DHdist, color=sample_source)) +
	ggtitle("HBond D-H Distance by Chemical Type")

if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
