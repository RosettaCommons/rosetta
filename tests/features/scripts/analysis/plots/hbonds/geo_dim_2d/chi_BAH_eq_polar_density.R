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
id = "chi_BAH_eq_polar_density",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
	geom.AHdist,
	geom.chi,
	acc_atoms.base_x AS bx, acc_atoms.base_y AS by, acc_atoms.base_z AS bz,
	acc_atoms.atm_x  AS ax, acc_atoms.atm_y  AS ay, acc_atoms.atm_z  AS az,
	don_atoms.atm_x  AS hx, don_atoms.atm_y  AS hy, don_atoms.atm_z  AS hz,
	acc_site.HBChemType AS acc_chem_type,
	don_site.HBChemType AS don_chem_type,
	CASE acc_site.HBChemType
		WHEN 'hbacc_IMD' THEN 'ring' WHEN 'hbacc_IME' THEN 'ring'
		WHEN 'hbacc_AHX' THEN 'sp3'  WHEN 'hbacc_HXL' THEN 'sp3'
		WHEN 'hbacc_CXA' THEN 'sp2'  WHEN 'hbacc_CXL' THEN 'sp2'
		WHEN 'hbacc_PBA' THEN 'bb_sp2'  END AS hybrid
FROM
	hbond_geom_coords AS geom,
	hbonds AS hb,
	hbond_sites AS don_site,
	hbond_sites AS acc_site,
	hbond_site_atoms AS don_atoms,
	hbond_site_atoms AS acc_atoms
WHERE
	hb.struct_id = geom.struct_id AND hb.hbond_id = geom.hbond_id AND
	hb.struct_id = don_site.struct_id AND hb.don_id = don_site.site_id AND
	hb.struct_id = acc_site.struct_id AND hb.acc_id = acc_site.site_id AND
	don_atoms.struct_id = hb.struct_id AND don_atoms.site_id = hb.don_id AND
	acc_atoms.struct_id = hb.struct_id AND acc_atoms.site_id = hb.acc_id AND
	ABS(don_site.resNum - acc_site.resNum) > 5;"

f <- query_sample_sources(sample_sources, sele)

f$don_chem_type_name <- factor(f$don_chem_type,
	levels = c("hbdon_IMD", "hbdon_IME", "hbdon_GDE", "hbdon_GDH",
		"hbdon_AHX", "hbdon_HXL", "hbdon_IND", "hbdon_AMO", "hbdon_CXA", "hbdon_PBA"),
	labels = c("dIMD: h", "dIME: h", "dGDE: r", "dGDH: r",
		"dAHX: y", "dHXL: s,t", "dIND: w", "dAMO: k", "dCXA: n,q", "dPBA: bb"))

f$acc_chem_type_name <- factor(f$acc_chem_type,
	levels = c("hbacc_IMD", "hbacc_IME", "hbacc_AHX", "hbacc_HXL",
		"hbacc_CXA", "hbacc_CXL", "hbacc_PBA"),
	labels = c("aIMD: h", "aIME: h", "aAHX: y", "aHXL: s,t",
		"aCXA: n,q", "aCXL: d,e", "aPBA: bb"))

f <- transform(f,
	cosBAH = vector_dotprod(
		vector_normalize(cbind(ax-bx, ay-by, az-bz)),
		vector_normalize(cbind(hx-ax, hy-ay, hz-az))))

#equal area projection
f <- transform(f,
	capx = 2*sin(acos(cosBAH)/2)*cos(chi),
	capy = 2*sin(acos(cosBAH)/2)*sin(chi))

capx_limits <- c(-1.5,1.5)
capy_limits <- c(-1.5,1.5)

max_BAH_angle = 110 # in degrees
abs_cap_limit <- 2*sin((max_BAH_angle*pi/180)/2)
capx_limits <- c(-abs_cap_limit, abs_cap_limit)
capy_limits <- c(-abs_cap_limit, abs_cap_limit)

narrow_output_formats <- transform(output_formats, width=height)

plot_parts <- list(
	theme_bw(),
	theme(panel.background=element_rect(fill="#00007F", colour="#00007F")),
	stat_density2d(
		aes(x=capx,y=capy, fill=..density..), geom="tile", contour=FALSE),
	geom_indicator(aes(indicator=counts), color="white"),
	polar_equal_area_grids_bw(),
	coord_equal(ratio=1),
	scale_fill_gradientn('Density', colours=jet.colors(15)),
	scale_x_continuous(limits=capx_limits),
	scale_y_continuous(limits=capy_limits),
	theme(
		axis.text.x=theme_blank(),
		axis.text.y=theme_blank(),
		axis.title.x=theme_blank(),
		axis.title.y=theme_blank(),
		axis.ticks.x = theme_blank(),
		axis.ticks.y = theme_blank()))


f <- ddply(f, c("sample_source", "hybrid"),
	transform, counts = length(sample_source))

plot_id <- "hbond_chi_BAH_eq_polar_density"
d_ply(f, .(sample_source, hybrid), function(sub_f){
	ss_id <- sub_f$sample_source[1]
	hybrid <- sub_f$hybrid[1]
	sub_plot_id <- paste(plot_id, hybrid, sep="_")
	ggplot(data=sub_f) + plot_parts +
		ggtitle(paste(hybrid, " Acceptor Hydrogen Bonds seq_sep > 5: chi vs BAH\nEqual Coordinate Projection  ss_id: ", ss_id, sep="")) +
		coord_equal(ratio=1) +
		polar_equal_area_grids_bw() +
		geom_indicator(aes(indicator=counts), color="white") +
		scale_x_continuous('', limits=capx_limits, breaks=c()) +
		scale_y_continuous('', limits=capy_limits, breaks=c()) +
		scale_fill_gradientn('log(Normalized\nDensity)', colours=jet.colors(15))
	save_plots(self, sub_plot_id, sample_sources[sample_sources$sample_source == ss_id,],
		output_dir, narrow_output_formats)
})


f <- ddply(f, c("sample_source", "acc_chem_type", "don_chem_type"),
	transform, counts = length(sample_source))

plot_id = "hbond_chi_BAH_eq_polar_density_by_chem_type"
l_ply(levels(f$sample_source), function(ss){
	ggplot(data=f[f$sample_source == ss,]) + plot_parts +
		facet_grid(acc_chem_type_name ~ don_chem_type_name) +
		ggtitle(paste("Hydrogen Bonds chi vs BAH Angles by Chemical Type Sequence Separation > 5\nEqual Coordinate Projection   Sample Source: ", ss, sep="")) +
		coord_equal(ratio=1) +
		geom_indicator(aes(indicator=counts), color="white") +
		scale_x_continuous('', limits=capx_limits, breaks=c()) +
		scale_y_continuous('', limits=capy_limits, breaks=c()) +
		scale_fill_gradientn('log(Normalized\nDensity)', colours=jet.colors(15))
	save_plots(self, plot_id, sample_sources[sample_sources$sample_source == ss,], output_dir, output_formats)
})

plot_id = "hbond_chi_BAH_eq_polar_density"
d_ply(f, .(sample_source, don_chem_type, acc_chem_type), function(sub_f){
	ss_id <- sub_f$sample_source[1]

	sub_plot_id <- paste(
		plot_id, sub_f$don_chem_type[1], sub_f$acc_chem_type[1], sep="_")

	ggplot(data=sub_f) + plot_parts +
		facet_grid(acc_chem_type_name ~ don_chem_type_name) +
		ggtitle(paste(
			"HBonds chi vs BAH Angles Sequence Separation > 5 ",
			"Don:", sub_f$don_chem_type_name[1], " ",
			"Acc:", sub_f$acc_chem_type_name[1], "\n",
			"Equal Coordinate Projection   Sample Source: ", ss_id, sep=""))
	save_plots(self, sub_plot_id, sample_sources[sample_sources$sample_source == ss_id,],
		output_dir, output_formats)
})


})) # end FeaturesAnalysis
