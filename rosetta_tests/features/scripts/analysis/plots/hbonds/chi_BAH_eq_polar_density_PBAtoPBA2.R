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
	geom.AHdist,
	geom.cosBAH,
	geom.chi,
  don_site.HBChemType AS don_chem_type,
  acc_site.HBChemType AS acc_chem_type
FROM
	hbond_geom_coords AS geom,
	hbonds AS hbond,
	hbond_sites AS don_site,
	hbond_sites AS acc_site
WHERE
  don_site.HBChemType == 'hbdon_PBA' AND acc_site.HBChemType == 'hbacc_PBA' AND
	hbond.struct_id = geom.struct_id AND hbond.hbond_id = geom.hbond_id AND
	hbond.struct_id = don_site.struct_id AND hbond.don_id = don_site.site_id AND
	hbond.struct_id = acc_site.struct_id AND hbond.acc_id = acc_site.site_id AND
	ABS(don_site.resNum - acc_site.resNum) > 5;";
f <- query_sample_sources(sample_sources, sele)

f <- ddply(f, c("sample_source"),
	transform, counts = length(sample_source))

#equal area projection
f <- transform(f,
	capx = 2*sin(acos(cosBAH)/2)*cos(chi),
	capy = 2*sin(acos(cosBAH)/2)*sin(chi))

##orthographic projection
#f <- transform(f,
#	capx = sin(acos(cosBAH))*cos(chi),
#	capy = sin(acos(cosBAH))*sin(chi))

#capx_limits <- range(f$capx);
#capy_limits <- range(f$capy)
capx_limits <- c(-1.5,1.5)
capy_limits <- capx_limits

narrow_output_formats <- transform(output_formats, width=height)

d_ply(f, .(sample_source), function(sub_f){
	ss_id <- sub_f$sample_source[1]
	plot_id = paste("chi_BAH_eq_polar_density_lr_bbbb", ss_id, sep="_")
	ggplot(data=sub_f) + theme_bw() +
		stat_density2d(aes(x=capx,y=capy, fill=..density..), geom="tile", contour=FALSE) +
		geom_indicator(aes(indicator=counts)) +
		facet_grid(acc_chem_type ~ don_chem_type) +
		opts(title = paste("Hydrogen Bonds chi vs BAH Angles with Sequence Separation > 5\nBackbone/Backbone Hydrogen Bonds\nEqual Coordinate Projection   Sample Source: ", ss_id, sep="")) +
		scale_x_continuous('2*sin(BAH/2) * cos(CHI)', limits=capx_limits, breaks=c(-1, 0, 1)) +
		scale_y_continuous('2*sin(BAH/2) * sin(CHI)', limits=capy_limits, breaks=c(-1, 0, 1)) +
		polar_equal_area_grids_bw() +
		coord_fixed(ratio = 1) +
		scale_fill_gradientn('Density', colour=jet.colors(10)) +
#        	opts(legend.position="bottom", legend.direction="horizontal")
	save_plots(plot_id, sample_sources[sample_sources$sample_source == ss_id,],
		output_dir, narrow_output_formats)
})
