# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

feature_analyses <- c(feature_analyses, new("FeaturesAnalysis",
id = "AHchi_AHD_eq_polar_density_salt_bridge_geometries",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

source("scripts/analysis/plots/hbonds/hbond_geo_dim_scales.R")

sele <-"
SELECT
	geom.cosAHD, geom.chi,
	don_site.HBChemType AS don_chem_type,
	acc_site.HBChemType AS acc_chem_type,
 	acc_atoms.base_x AS abx, acc_atoms.base_y AS aby, acc_atoms.base_z AS abz, -- acceptor base atom
	acc_atoms.atm_x  AS ax,  acc_atoms.atm_y  AS ay,  acc_atoms.atm_z  AS az,  -- acceptor atom
	don_atoms.atm_x  AS hx,  don_atoms.atm_y  AS hy,  don_atoms.atm_z  AS hz,  -- hydrogen atom
	don_atoms.base_x AS dx,  don_atoms.base_y AS dy,  don_atoms.base_z AS dz  -- donor atom
FROM
	hbond_geom_coords AS geom,
	hbonds AS hbond,
	hbond_sites AS don_site,
	hbond_sites AS acc_site,
	hbond_site_atoms AS don_atoms, hbond_site_atoms AS acc_atoms,
	hbond_sites_pdb AS don_pdb, hbond_sites_pdb AS acc_pdb
WHERE
	hbond.struct_id = geom.struct_id AND hbond.hbond_id = geom.hbond_id AND
	hbond.struct_id = don_site.struct_id AND hbond.don_id = don_site.site_id AND
	hbond.struct_id = acc_site.struct_id AND hbond.acc_id = acc_site.site_id AND
	don_atoms.struct_id = hbond.struct_id AND don_atoms.site_id = hbond.don_id AND
	acc_atoms.struct_id = hbond.struct_id AND acc_atoms.site_id = hbond.acc_id AND
  don_pdb.struct_id = hbond.struct_id AND don_pdb.site_id = don_site.site_id AND
  acc_pdb.struct_id = hbond.struct_id AND acc_pdb.site_id = acc_site.site_id AND
	don_pdb.heavy_atom_temperature < 30 AND acc_pdb.heavy_atom_temperature < 30 AND
	ABS(don_site.resNum - acc_site.resNum) > 5 AND
	(don_site.HBChemType == 'hbdon_GDE' OR don_site.HBChemType == 'hbdon_GDH' OR
	don_site.HBChemType == 'hbdon_AMO' OR	don_site.HBChemType == 'hbdon_IMD' OR
	don_site.HBChemType == 'hbdon_IME') AND
	(acc_site.HBChemType == 'hbacc_CXA' OR acc_site.HBChemType == 'hbacc_CXL');"
f <- query_sample_sources(sample_sources, sele)

f$don_chem_type_name <- don_chem_type_name_wrap(f$don_chem_type)
f$acc_chem_type_name <- acc_chem_type_name_wrap(f$acc_chem_type)

f <- transform(f,
	AHchi = vector_dihedral(
		cbind(abx, aby, abz), cbind(ax, ay, az),
		cbind(hx, hy, hz), cbind(dx, dy, dz)))

#equal area projection
f <- transform(f,
	capx = 2*sin(acos(cosAHD)/2)*cos(AHchi),
	capy = 2*sin(acos(cosAHD)/2)*sin(AHchi))

capx_limits <- c(-1.5,1.5)
capy_limits <- capx_limits

plot_parts <- list(
	theme_bw(),
	theme(panel.background=element_rect(fill="#00007F", colour="#00007F")),
	stat_density2d(
		aes(x=capx,y=capy, fill=..density..), geom="tile", contour=FALSE),
	geom_indicator(aes(indicator=counts), color="white"),
	polar_equal_area_grids_bw(),
	coord_fixed(ratio = 1),
	scale_fill_gradientn('Density', colour=jet.colors(15)),
	scale_x_continuous(limits=capx_limits),
	scale_y_continuous(limits=capy_limits),
	theme(
		axis.text.x=theme_blank(),
		axis.text.y=theme_blank(),
		axis.title.x=theme_blank(),
		axis.title.y=theme_blank(),
		axis.ticks.x = theme_blank(),
		axis.ticks.y = theme_blank()))

narrow_output_formats <- transform(output_formats, width=height)

d_ply(f, .(sample_source), function(sub_f){
	ss_id <- sub_f$sample_source[1]

	d_ply(sub_f, .(don_chem_type), function(sub_sub_f){
		don_chem_type <- sub_sub_f$don_chem_type[1]
		don_chem_type_name <- sub_sub_f$don_chem_type_name[1]
		sub_sub_f$counts <- nrow(sub_sub_f)
		plot_id = paste("AHchi_AHD_eq_polar_density", don_chem_type, ss_id, sep="_")
		ggplot(data=sub_sub_f) + plot_parts +
			ggtitle(paste(
				"Hydrogen Bonds AHchi vs AHD Angle, ", don_chem_type_name, "; SeqSep > 5, BFact < 30\n",
				"Equal Coordinate Projection   Sample Source: ", ss_id, sep=""))
		save_plots(self, plot_id, sample_sources, output_dir, output_formats)

	})

	d_ply(sub_f, .(don_chem_type, acc_chem_type), function(sub_sub_f){
		don_chem_type <- sub_sub_f$don_chem_type[1]
		don_chem_type_name <- sub_sub_f$don_chem_type_name[1]
		acc_chem_type <- sub_sub_f$acc_chem_type[1]
		acc_chem_type_name <- sub_sub_f$acc_chem_type_name[1]
		sub_sub_f$counts <- nrow(sub_sub_f)
		plot_id = paste("AHchi_AHD_eq_polar_density", don_chem_type, acc_chem_type, ss_id, sep="_")
		ggplot(data=sub_sub_f) + plot_parts +
			ggtitle(paste(
				"Hydrogen Bonds AHchi vs AHD Angle, ", don_chem_type_name, " ", acc_chem_type_name, ";",
				"SeqSep > 5, BFact < 30\n",
				"Equal Coordinate Projection   Sample Source: ", ss_id, sep=""))
		save_plots(self, plot_id, sample_sources, output_dir, narrow_output_formats)

	})

})


})) # end FeaturesAnalysis
