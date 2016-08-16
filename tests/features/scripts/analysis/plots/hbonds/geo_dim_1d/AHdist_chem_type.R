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
id = "AHdist_chem_type",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
	geom.AHdist,
	don.HBChemType AS don_chem_type, acc.HBChemType AS acc_chem_type
FROM
	hbonds AS hb,
	hbond_geom_coords AS geom,
	hbond_sites AS don, hbond_sites AS acc,
	hbond_sites_pdb AS don_pdb, hbond_sites_pdb AS acc_pdb
WHERE
	geom.struct_id = hb.struct_id AND geom.hbond_id = hb.hbond_id AND
	don.struct_id = hb.struct_id AND don.site_id = hb.don_id AND
	acc.struct_id = hb.struct_id AND acc.site_id = hb.acc_id AND
	don_pdb.struct_id = hb.struct_id AND don_pdb.site_id = hb.don_id AND
	don_pdb.heavy_atom_temperature < 30 AND
	acc_pdb.struct_id = hb.struct_id AND acc_pdb.site_id = hb.acc_id AND
	acc_pdb.heavy_atom_temperature < 30 AND
  abs(don.resNum - acc.resNum ) > 5;"

f <- query_sample_sources(sample_sources, sele)

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

f <- na.omit(f, method="r")

dens <- estimate_density_1d(
	f, c("sample_source", "acc_chem_type", "don_chem_type"),
	"AHdist", weight_fun = radial_3d_normalization)

#library(earthmovdist)
#z <- comparison_statistics(f, c("don_chem_type", "acc_chem_type"), "AHdist", earth_mover_distance_L1)

#z <- smooth_comparison_statistics(dens, c("don_chem_type", "acc_chem_type"), smooth_kl_divergence)

#plot_id <- "emd_sample_size_correlation"
#p <- ggplot(data=z) + theme_bw() +
#	geom_point(aes(x=log(sample.size), y=log(statistic))) +
#	stat_smooth(aes(x=log(sample.size), y=log(statistic)), method="lm", se=F) +
#	ggtitle("Earth Mover's Distance as a function of Sample Size") +
#	scale_x_continuous("log(Sample Size)") +
#	scale_y_continuous("log(Earth Mover's Distance)")
#save_plots(self, plot_id, sample_sources, output_dir, output_formats)
#
#plot_id <- "emd_sample_size_correlation"
#p <- ggplot(data=z[!(z$don_chem_type == "hbdon_PBA" & z$acc_chem_type == "hbacc_PBA"),]) + theme_bw() +
#	geom_point(aes(x=sample.size, y=statistic)) +
#	stat_smooth(aes(x=sample.size, y=statistic), method="lm", se=F) +
#	ggtitle("Earth Mover's Distance as a function of Sample Size") +
#	scale_x_continuous("Sample Size") +
#	scale_y_continuous("Earth Mover's Distance")
#save_plots(self, plot_id, sample_sources, output_dir, output_formats)



plot_id <- "hbond_AHdist_chem_type"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x, y=y, colour=sample_source)) +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
#	geom_indicator(data=z, aes(indicator=round(statistic,2), colour=new_sample_source), xpos="left") +
	facet_grid(don_chem_type ~ acc_chem_type) +
	ggtitle("HBond A-H Distance by Chemical Type, SeqSep > 5, B-Factor < 30\nnormalized for equal weight per unit distance") +
	scale_y_continuous("FeatureDensity", limits=c(0,7.5), breaks=c(1,3,5,7)) +
	scale_x_continuous(expression(paste('Acceptor -- Proton Distance (', ring(A), ')')), limits=c(1.4,2.7), breaks=c(1.6, 1.9, 2.2, 2.5))

if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
