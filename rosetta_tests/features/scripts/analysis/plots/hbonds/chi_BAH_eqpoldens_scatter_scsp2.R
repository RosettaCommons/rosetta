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
	don_site.HBChemType AS don_chem_type
FROM
	hbond_geom_coords AS geom,
	hbonds AS hbond,
	hbond_sites AS don_site,
	hbond_sites AS acc_site
WHERE
  don_site.HBChemType != 'hbdon_PBA' AND
	( acc_site.HBChemType == 'hbacc_CXL' OR acc_site.HBChemType == 'hbacc_CXA' ) AND
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

#capx_limits <- range(f$capx);
#capy_limits <- range(f$capy)
capx_limits <- c(-1.5,1.5)
capy_limits <- capx_limits

f$weight <- radial_3d_normalization(f$AHdist)

jet.colors <-
  colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

plot_id = "chi_BAH_eqpoldens_and_scatter_scsc_to_sp2"
#l_ply(levels(f$sample_source), function(ss){
#	ggplot(data=f) + theme_bw() +
#		polar_equal_area_grids_bw +
#		#geom_bin2d(aes(x=capx, y=capy, fill=log(..count..), weight=weight), binwidth=c(.06, .06)) +
#    stat_density2d(aes(x=capx,y=capy, fill=..density..), geom="tile", contour=FALSE ) +
#		#facet_grid(acc_chem_type ~ don_chem_type) +
#		opts(title = paste("Hydrogen Bonds chi vs BAH Angles with Sequence Separation > 5\nSidechain Hydroxyl Donors to Sidechain Carboxyl Acceptors\nEqual Coordinate Projection   Sample Source: ", ss, sep="")) +
#		scale_x_continuous('2*sin(BAH/2) * cos(CHI)', limits=capx_limits, breaks=c(-1, 0, 1)) +
#		scale_y_continuous('2*sin(BAH/2) * sin(CHI)', limits=capy_limits, breaks=c(-1, 0, 1)) +
#		coord_fixed(ratio = 1) +
#		scale_fill_gradientn('Density', colour=jet.colors(10)) +
##        	opts(legend.position="bottom", legend.direction="horizontal")
#	save_plots(plot_id, sample_sources[sample_sources$sample_source == ss,], output_dir, output_formats)
#})

f_first <- f[ f$sample_source == levels(sample_sources$sample_source)[1], ]
f_second <- f[ f$sample_source == levels(sample_sources$sample_source)[2], ]

ggplot(data=f_first) + theme_bw() +
		polar_equal_area_grids_bw +
		#geom_bin2d(aes(x=capx, y=capy, fill=log(..count..), weight=weight), binwidth=c(.06, .06)) +
    stat_density2d(aes(x=capx,y=capy, fill=..density..), geom="tile", contour=FALSE ) +
		#facet_grid(acc_chem_type ~ don_chem_type) +
		opts(title = paste("Hydrogen Bonds chi vs BAH Angles with Sequence Separation > 5\nSidechain Donors to Sidechain sp2 Acceptors, Equal Coordinate Projection\nNatives (density) vs Interface Designs (white circles)")) +
		scale_x_continuous('2*sin(BAH/2) * cos(CHI)', limits=capx_limits, breaks=c(-1, 0, 1)) +
		scale_y_continuous('2*sin(BAH/2) * sin(CHI)', limits=capy_limits, breaks=c(-1, 0, 1)) +
		coord_fixed(ratio = 1) +
		scale_fill_gradientn('Density', colour=jet.colors(10)) +
		geom_point( data=f_second,aes(x=capx,y=capy),colour="white",size=2)
#        	opts(legend.position="bottom", legend.direction="horizontal")
save_plots(plot_id, sample_sources, output_dir, output_formats)

