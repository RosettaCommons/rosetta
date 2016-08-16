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
id = "AHdist_by_chem_type_signal",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


sele <-"
SELECT
  geom.AHdist,
  acc_site.HBChemType AS acc_chem_type,
  don_site.HBChemType AS don_chem_type
FROM
  hbond_geom_coords AS geom,
  hbonds AS hbond,
  hbond_sites AS don_site,
  hbond_sites AS acc_site
WHERE
  hbond.struct_id = geom.struct_id AND
  hbond.hbond_id =  geom.hbond_id AND
  hbond.struct_id = don_site.struct_id AND
  hbond.don_id = don_site.site_id AND
  hbond.struct_id = acc_site.struct_id AND
  hbond.acc_id = acc_site.site_id;"
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
  f, c("sample_source", "acc_chem_type", "don_chem_type"),
  "AHdist", weight_fun = radial_3d_normalization)

ref_dens <- dens[dens$sample_source == sample_sources$sample_source[1],]
new_dens <- dens[dens$sample_source != sample_sources$sample_source[1],]


ref_ss_id <- ref_dens$sample_source[1]
d_ply(new_dens, .(sample_source), function(sub_dens){
	new_ss_id <- new_dens$sample_source[1]

	new_dens$ref_y <- ref_dens$y + 1

	plot_id <- "AHdist_chem_type"
	p <- ggplot(data=new_dens) + theme_bw() +
		geom_line(aes(x=x, y=log(y/ref_y)*ref_y, colour=sample_source)) +
		geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
		facet_grid(don_chem_type ~ acc_chem_type) +
		ggtitle(paste("Hydrogen Bonds A-H Distance Signal Ratio by Chemical Type\nNormalized for equal weight per unit distance\nref_ss: ", ref_ss_id, " new_ss: ", new_ss_id, sep="")) +
		scale_y_continuous("log(Signal)") +
		scale_x_continuous(
			expression(paste('Acceptor -- Proton Distance (', ring(A), ')')),
			limits=c(1.4,2.7), breaks=c(1.6, 1.9, 2.2, 2.6))
	if(nrow(sample_sources) <= 3){
		p <- p + theme(legend.position="bottom", legend.direction="horizontal")
	}
	ref_new_sample_sources = sample_sources[
		sample_sources$sample_source == new_ss_id ||
		sample_sources$sample_source == ref_ss_id,]
	save_plots(self, plot_id, ref_new_sample_sources, output_dir, output_formats)
})

})) # end FeaturesAnalysis
