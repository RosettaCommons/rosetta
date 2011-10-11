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
	acc_site.HBChemType AS acc_chem_type,
	don_site.HBChemType AS don_chem_type,

  CASE acc_site.HBChemType
		WHEN 'hbacc_IMD' THEN 'ring' WHEN 'hbacc_IME' THEN 'ring'
		WHEN 'hbacc_AHX' THEN 'sp3'  WHEN 'hbacc_HXL' THEN 'sp3'
		WHEN 'hbacc_CXA' THEN 'sp2'  WHEN 'hbacc_CXL' THEN 'sp2'
		WHEN 'hbacc_PBA' THEN 'bb_sp2'  END AS hybrid
FROM
	hbond_geom_coords AS geom,
	hbonds AS hbond,
	hbond_sites AS don_site,
	hbond_sites AS acc_site
WHERE
	hbond.struct_id = geom.struct_id AND hbond.hbond_id = geom.hbond_id AND
	hbond.struct_id = don_site.struct_id AND hbond.don_id = don_site.site_id AND
	hbond.struct_id = acc_site.struct_id AND hbond.acc_id = acc_site.site_id AND
	ABS(don_site.resNum - acc_site.resNum) > 5;";
f <- query_sample_sources(sample_sources, sele)

#equal area projection
f <- transform(f,
	capx = 2*sin(acos(cosBAH)/2)*cos(chi),
	capy = 2*sin(acos(cosBAH)/2)*sin(chi))

##orthographic projection
#f <- transform(f,
#	capx = sin(acos(cosBAH))*cos(chi),
#	capy = sin(acos(cosBAH))*sin(chi))

capx_limits <- range(f$capx); capy_limits <- range(f$capy)

narrow_output_formats <- transform(output_formats, width=height)

f <- ddply(f, c("sample_source", "hybrid"),
	transform, counts = length(sample_source))

plot_id = "chi_BAH_eq_polar_density"
d_ply(f, .(sample_source, hybrid), function(sub_f){
	ss_id <- sub_f$sample_source[1]
	hybrid <- sub_f$hybrid[1]
	sub_plot_id <- paste(plot_id, hybrid, ss_id, sep="_")
	ggplot(data=sub_f) + theme_bw() +
		polar_equal_area_grids_bw() +
		geom_indicator(aes(indicator=counts)) +
		stat_bin2d(aes(x=capx, y=capy, fill=log(..density..)), binwidth=c(.01, .01)) +
		opts(title = paste(hybrid, " Acceptor Hydrogen Bonds seq_sep > 5: chi vs BAH\nEqual Coordinate Projection  ss_id: ", ss_id, sep="")) +
		coord_equal() +
		scale_x_continuous('2*sin(BAH/2) * cos(CHI)', limits=capx_limits, breaks=c(-1, 0, 1)) +
		scale_y_continuous('2*sin(BAH/2) * sin(CHI)', limits=capy_limits, breaks=c(-1, 0, 1)) +
		scale_fill_gradientn('log(Normalized\nDensity)', colour=jet.colors(15))
	save_plots(sub_plot_id, sample_sources[sample_sources$sample_source == ss_id,],
		output_dir, narrow_output_formats)
})


f <- ddply(f, c("sample_source", "acc_chem_type", "don_chem_type"),
	transform, counts = length(sample_source))

plot_id = "chi_BAH_eq_polar_density_by_chem_type"
l_ply(levels(f$sample_source), function(ss){
	ggplot(data=f[f$sample_source == ss,]) + theme_bw() +
		polar_equal_area_grids_bw() +
		geom_indicator(aes(indicator=counts)) +
		geom_bin2d(aes(x=capx, y=capy, fill=log(..density..)), binwidth=c(.06, .06)) +
		facet_grid(acc_chem_type ~ don_chem_type) +
		opts(title = paste("Hydrogen Bonds chi vs BAH Angles by Chemical Type Sequence Separation > 5\nEqual Coordinate Projection   Sample Source: ", ss, sep="")) +
		scale_x_continuous('2*sin(BAH/2) * cos(CHI)', limits=capx_limits, breaks=c(-1, 0, 1)) +
		scale_y_continuous('2*sin(BAH/2) * sin(CHI)', limits=capy_limits, breaks=c(-1, 0, 1)) +
		scale_fill_gradientn('log(Normalized\nDensity)', colour=jet.colors(15)) +
#		opts(legend.position="bottom", legend.direction="horizontal")
	save_plots(plot_id, sample_sources[sample_sources$sample_source == ss,], output_dir, output_formats)
})
