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
id = "angle_dependence_on_AHdist",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
	geom.AHdist, geom.cosBAH, geom.cosAHD, geom.chi,
	don.HBChemType AS don_chem_type, acc.HBChemType AS acc_chem_type
FROM
	hbond_geom_coords AS geom,
	hbonds AS hbond,
	hbond_sites AS don, hbond_sites AS acc
WHERE
	hbond.struct_id = geom.struct_id AND hbond.hbond_id = geom.hbond_id AND
	hbond.struct_id = don.struct_id  AND hbond.don_id   = don.site_id   AND
	hbond.struct_id = acc.struct_id  AND hbond.acc_id   = acc.site_id;"

f <- query_sample_sources(sample_sources, sele)
f$AHdist_quantile <- cut(f$AHdist, quantile(f$AHdist, seq(0,1,by=.5)))

f <- na.omit(f, method="r")

#chi goes from -pi to pi
#rotated chi goes from 0 to 360, but what was pi is now 180
f$chi <- (f$chi*180/pi)%%360

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

plot_parts <- list(
	theme_bw(),
	geom_line(aes(colour=AHdist_quantile)),
	geom_indicator(aes(indicator=counts, colour=AHdist_quantile, group=AHdist_quantile)),
	facet_grid(don_chem_type ~ acc_chem_type),
	scale_y_continuous("log(FeatureDensity + 1)"))

d_ply(f, .variables=("sample_source"), function(sub_f){
	ss_id <- sub_f[1,"sample_source"]

	plot_id <- paste("cosBAH_AHdist_quantiles_", ss_id, "_by_chem_type", sep="")
	dens <- estimate_density_1d(sub_f,
		c("acc_chem_type", "don_chem_type", "AHdist_quantile"), "cosBAH")
	ggplot(data=dens, aes(acos(x)*180/pi, log(y+1))) + plot_parts +
			scale_x_continuous('Base -- Acceptor -- Hydrogen (degrees)') +
			ggtitle(paste("Hydrogen Bond Angle at the Acceptor by Chemical Type and AHdist Quantile\nss_id: ", ss_id, sep=""))
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)

	plot_id <- paste("cosAHD_AHdist_quantiles_", ss_id, "_by_chem_type", sep="")
	dens <- estimate_density_1d_reflect_boundary(sub_f,
		c("acc_chem_type", "don_chem_type", "AHdist_quantile"), "cosAHD",
		reflect_left=TRUE, left_boundary=0)
	ggplot(data=dens, aes(acos(x)*180/pi, log(y+1))) + plot_parts +
		scale_x_continuous('Acceptor -- Hydrogen -- Donor (degrees)') +
		ggtitle(paste("Hydrogen Bond Angle at the Hydrogen by Chemical Type and AHdist Quantile\nss_id: ", ss_id, sep=""))
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)

	plot_id <- paste("chi_AHdist_quantiles_", ss_id, "_by_chem_type", sep="")
	dens <- estimate_density_1d(sub_f,
		c("acc_chem_type", "don_chem_type", "AHdist_quantile"), "chi")
	ggplot(data=dens, aes(x, log(y+1))) + plot_parts +
		scale_x_continuous('Acceptor Base -- Acceptor Torsion (degrees)', breaks=c(90,180, 270)) +
		ggtitle(paste("Hydrogen Bond Base-Acceptor Torsion Angle by Chemical Type and AHdist Quantile\nss_id: ", ss_id, sep=""))
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
})


})) # end FeaturesAnalysis
