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
id = "AHdist_by_environmental_weight",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
	geom.AHdist,
	hb.envWeight
FROM
	hbonds AS hb,
	hbond_geom_coords AS geom,
	hbond_sites AS don, hbond_sites AS acc,
	hbond_sites_pdb AS don_pdb, hbond_sites_pdb AS acc_pdb
WHERE
	geom.struct_id = hb.struct_id AND geom.hbond_id = hb.hbond_id AND
	don.struct_id = hb.struct_id AND don.site_id = hb.don_id AND
	acc.struct_id = hb.struct_id AND acc.site_id = hb.acc_id AND
	don.HBChemType != 'hbdon_PBA' AND acc.HBChemType != 'hbacc_PBA' AND
	don_pdb.struct_id = hb.struct_id AND don_pdb.site_id = hb.don_id AND
	don_pdb.heavy_atom_temperature < 30 AND
	acc_pdb.struct_id = hb.struct_id AND acc_pdb.site_id = hb.acc_id AND
	acc_pdb.heavy_atom_temperature < 30;"
f <- query_sample_sources(sample_sources, sele)


f$Burial <- cut(f$envWeight, breaks=c(0, .5, .9999, 1), labels=c("Exposed (envWeigth < .5)", "Partial (.5 <= envWeight < 1)", "Buried (1 <= envWeight)"))

dens <- estimate_density_1d(
	f, c("sample_source", "Burial"),
	"AHdist", weight_fun = radial_3d_normalization)

plot_id <- "hbond_AHdist_by_environmental_weight"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x, y=y, colour=Burial), size=1.6) +
	geom_indicator(aes(indicator=counts, colour=Burial, group=Burial)) +
	facet_wrap(~sample_source) +
	ggtitle("Sidechain-Sidechain HBond A-H Distance by Neighbor Count based Solvent Exposure\n B-Factor < 30 normalized for equal weight per unit distance") +
	scale_y_continuous("FeatureDensity", limits=c(0,7), breaks=c(0,1,2,3,4,5,6,7)) +
	scale_x_continuous(expression(paste('Acceptor -- Proton Distance (', ring(A), ')')), limits=c(1.4,2.7), breaks=c(1.6, 1.9, 2.2, 2.6))

if(nrow(sample_sources) < 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
} else if(nrow(sample_sources) == 3){
	p <- p + theme(legend.position=c(.7, .35))
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis

