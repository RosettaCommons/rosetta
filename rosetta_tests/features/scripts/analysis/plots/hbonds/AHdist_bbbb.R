# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

check_setup()

plot_id <- "AHdist_bbbb"

sele <-"
SELECT
	geom.AHdist,
	CASE don_site.resNum - acc_site.resNum
		WHEN -1 THEN '-1' WHEN -2 THEN '-2' WHEN -3 THEN '-3' WHEN -4 THEN '-4'
		WHEN 1 THEN '1' WHEN 2 THEN '2' WHEN 3 THEN '3' WHEN 4 THEN '4'
		ELSE 'long' END AS seq_sep
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
	hbond.acc_id = acc_site.site_id AND
	acc_site.HBChemType == 'hbacc_PBA' AND
	don_site.HBChemType == 'hbdon_PBA';"
f <- query_sample_sources(sample_sources, sele)

dens <- estimate_density_1d(
	f, c("sample_source", "seq_sep"),
	"AHdist", weight_fun = radial_3d_normalization)

p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x, y=y, colour=sample_source)) +
	geom_indicator(aes(indicator=counts, colour=sample_source)) +
	facet_wrap( ~ seq_sep ) +
	opts(title = "BB/BB Hydrogen Bonds A-H Distance by Sequence Separation\n(donres - accres) normalized for equal weight per unit distance") +
	scale_y_continuous("FeatureDensity)", limits=c(0,6), breaks=c(1,3,5)) +
	scale_x_continuous(expression(paste('Acceptor -- Proton Distance (', ring(A), ')')), limits=c(1.4,2.7), breaks=c(1.6, 1.9, 2.2, 2.6))

if(nrow(sample_sources) <= 3){
	p <- p + opts(legend.position="bottom", legend.direction="horizontal")
}

save_plots(plot_id, sample_sources, output_dir, output_formats)
