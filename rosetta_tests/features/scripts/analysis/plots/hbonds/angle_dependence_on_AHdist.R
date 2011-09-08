# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

check_setup()

sele <-"
SELECT
	geom.AHdist, geom.cosBAH, geom.cosAHD, geom.chi,
	REPLACE(don.HBChemType, 'hbdon_', '') AS d,
	REPLACE(acc.HBChemType, 'hbacc_', '') AS a
FROM
	hbond_geom_coords AS geom,
	hbonds AS hbond,
	hbond_sites AS don,
	hbond_sites AS acc
WHERE
	hbond.struct_id = geom.struct_id AND hbond.hbond_id = geom.hbond_id AND
	hbond.struct_id = don.struct_id  AND hbond.don_id   = don.site_id   AND
	hbond.struct_id = acc.struct_id  AND hbond.acc_id   = acc.site_id;"

f <- query_sample_sources(sample_sources, sele)
f$AHdist_quantile <- cut(f$AHdist, quantile(f$AHdist, seq(0,1,by=.5)))

#chi goes from -pi to pi
#rotated chi goes from 0 to 360, but what was pi is now 180
f$chi <- (f$chi*180/pi)%%360

plot_parts <- list(
	theme_bw(),
	geom_line(aes(colour=AHdist_quantile)),
	geom_indicator(aes(indicator=counts, colour=AHdist_quantile)),
	facet_grid(d ~ a),
	scale_y_continuous("log(FeatureDensity + 1)"))

d_ply(f, .variables=("sample_source"), function(sub_f){
	ss_id <- sub_f[1,"sample_source"]

	plot_id <- paste("cosBAH_AHdist_quantiles_", ss_id, "_by_chem_type", sep="")
	dens <- estimate_density_1d(sub_f, c("a", "d", "AHdist_quantile"), "cosBAH")
	ggplot(data=dens, aes(cos(x)*180/pi, log(y+1))) + plot_parts +
			scale_x_continuous('Base -- Acceptor -- Hydrogen (degrees)') +
			opts(title = paste("Hydrogen Bond Angle at the Acceptor by Chemical Type and AHdist Quantile: ",
				ss_id, "\nnormalized for equal weight per unit distance", sep=""))
	save_plots(plot_id, sample_sources, output_dir, output_formats)

	plot_id <- paste("cosAHD_AHdist_quantiles_", ss_id, "_by_chem_type", sep="")
	dens <- estimate_density_1d(sub_f, c("a", "d", "AHdist_quantile"), "cosAHD")
	ggplot(data=dens, aes(cos(x)*180/pi, log(y+1))) + plot_parts +
		scale_x_continuous('Acceptor -- Hydrogen -- Donor (degrees)') +
		opts(title = paste("Hydrogen Bond Angle at the Hydrogen by Chemical Type and AHdist Quantile: ",
			ss_id, "\nnormalized for equal weight per unit distance", sep=""))
	save_plots(plot_id, sample_sources, output_dir, output_formats)

	plot_id <- paste("chi_AHdist_quantiles_", ss_id, "_by_chem_type", sep="")
	dens <- estimate_density_1d(sub_f, c("a", "d", "AHdist_quantile"), "chi")
	ggplot(data=dens, aes(x, log(y+1))) + plot_parts +
		scale_x_continuous('Acceptor Base -- Acceptor Torsion (degrees)', breaks=c(90,270)) +
		opts(title = paste("Hydrogen Bond Base-Acceptor Torsion Angle by Chemical Type and AHdist Quantile: ",
			ss_id, "\nnormalized for equal weight per unit distance", sep=""))
	save_plots(plot_id, sample_sources, output_dir, output_formats)
})
