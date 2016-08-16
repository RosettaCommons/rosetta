# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

feature_analyses <- c(feature_analyses, new("FeaturesAnalysis",
id = "AHchi_AHD_eq_polar_density_PBAtoPBA",
author = "Matthew O'Meara",

brief_description = "",

feature_reporter_dependencies = c("HBondFeatures"),

run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
	geom.cosAHD,
	don_site.HBChemType AS don_chem_type,
	acc_site.HBChemType  acc_chem_type,
  don_ss.dssp AS don_dssp,
 	acc_atoms.base_x AS abx, acc_atoms.base_y AS aby, acc_atoms.base_z AS abz, -- acceptor base atom
	acc_atoms.atm_x  AS ax,  acc_atoms.atm_y  AS ay,  acc_atoms.atm_z  AS az,  -- acceptor atom
	don_atoms.atm_x  AS hx,  don_atoms.atm_y  AS hy,  don_atoms.atm_z  AS hz,  -- hydrogen atom
	don_atoms.base_x AS dx,  don_atoms.base_y AS dy,  don_atoms.base_z AS dz   -- donor atom
FROM
	hbond_geom_coords AS geom,
	hbonds AS hbond,
	hbond_sites AS don_site,
	hbond_sites AS acc_site,
	residue_secondary_structure AS don_ss,
	hbond_site_atoms AS don_atoms, hbond_site_atoms AS acc_atoms
WHERE
	don_site.HBChemType == 'hbdon_PBA' AND acc_site.HBChemType == 'hbacc_PBA' AND
	hbond.struct_id = geom.struct_id AND hbond.hbond_id = geom.hbond_id AND
	hbond.struct_id = don_site.struct_id AND hbond.don_id = don_site.site_id AND
	hbond.struct_id = acc_site.struct_id AND hbond.acc_id = acc_site.site_id AND
	don_atoms.struct_id = hbond.struct_id AND don_atoms.site_id = hbond.don_id AND
	acc_atoms.struct_id = hbond.struct_id AND acc_atoms.site_id = hbond.acc_id AND
  don_ss.struct_id = hbond.struct_id AND don_ss.resNum = don_site.resNum;";

f <- query_sample_sources(sample_sources, sele)

f <- ddply(f, .(sample_source, don_dssp), transform, counts = length(sample_source))

# reorder factors for better layout
f$don_dssp <- factor(f$don_dssp,
	levels = c(
		"H", "E", "T",
		"G", "B", "S",
		"I", " "),
	labels = c(
		'H: a-Helix',    'E: b-Sheet',  'T: HB Turn',
		'G: 3/10 Helix', 'B: b-Bridge', 'S: Bend',
		'I: pi-Helix',   'Irregular'))

# remove non-canonical amino acids
f <- na.omit(f, method="r")

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

narrow_output_formats <- transform(output_formats, width=height*.8)

d_ply(f, .(sample_source), function(sub_f){
	ss_id <- sub_f$sample_source[1]
	plot_id = paste("AHchi_AHD_eq_polar_density_bbbb_by_sec_struct", ss_id, sep="_")
	ggplot(data=sub_f) + theme_bw() +
		theme(panel.background=element_rect(fill="#00007F", colour="#00007F")) +
		stat_density2d(
			aes(x=capx,y=capy, fill=..density..), geom="tile", contour=FALSE) +
		polar_equal_area_grids_bw() +
		geom_indicator(aes(indicator=counts)) +
		facet_wrap(~don_dssp) +
		ggtitle =
			paste("Backbone Backbone Hydrogen Bonds AHchi vs AHD Angles\n",
			"Equal Coordinate Projection   Sample Source: ", ss_id, sep="")) +
		scale_x_continuous('', limits=capx_limits, breaks=c(-1, 0, 1)) +
		scale_y_continuous('', limits=capy_limits, breaks=c(-1, 0, 1)) +
		coord_fixed(ratio = 1) +
		scale_fill_gradientn('Density', colours=jet.colors(10))
	save_plots(self, plot_id, sample_sources[sample_sources$sample_source == ss_id,],
		output_dir, narrow_output_formats)
})

})) # end FeaturesAnalysis
