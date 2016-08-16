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
id = "OHacceptor_chi",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
	geom.cosBAH,
	geom.chi,
  CASE acc_site.HBChemType
		WHEN 'hbacc_IMD' THEN 'ring' WHEN 'hbacc_IME' THEN 'ring'
		WHEN 'hbacc_AHX' THEN 'sp3'  WHEN 'hbacc_HXL' THEN 'sp3'
		WHEN 'hbacc_CXA' THEN 'sp2'  WHEN 'hbacc_CXL' THEN 'sp2'
		WHEN 'hbacc_PBA' THEN 'bb_sp2'  END AS hybrid,
	acc_site.HBChemType AS acc_chem_type, don_site.HBChemType AS don_chem_type,
	acc_atoms.base2_x AS ab2x, acc_atoms.base2_y AS ab2y, acc_atoms.base2_z AS ab2z, -- acceptor base 2 atom
 	acc_atoms.base_x  AS abx,  acc_atoms.base_y  AS aby,  acc_atoms.base_z  AS abz,  -- acceptor base atom
	acc_atoms.atm_x   AS ax,   acc_atoms.atm_y   AS ay,   acc_atoms.atm_z   AS az,   -- acceptor atom
	don_atoms.atm_x   AS hx,   don_atoms.atm_y   AS hy,   don_atoms.atm_z   AS hz    -- hydrogen atom
FROM
	hbond_geom_coords AS geom,
	hbonds AS hbond,
	hbond_sites AS don_site, hbond_sites AS acc_site,
	hbond_site_atoms AS don_atoms, hbond_site_atoms AS acc_atoms
WHERE
	hbond.struct_id = geom.struct_id AND hbond.hbond_id = geom.hbond_id AND
	hbond.struct_id = don_site.struct_id AND hbond.don_id = don_site.site_id AND
	hbond.struct_id = acc_site.struct_id AND hbond.acc_id = acc_site.site_id AND
	don_atoms.struct_id = hbond.struct_id AND don_atoms.site_id = hbond.don_id AND
	acc_atoms.struct_id = hbond.struct_id AND acc_atoms.site_id = hbond.acc_id AND
	(don_site.HBChemType = 'hbdon_AXL' OR don_site.HBChemType = 'hbdon_HXL');"

f <- query_sample_sources(sample_sources, sele)

alt_chi_dihedral_angle <- function(ab2, ab, a, h){
	alt_ab <- (ab + ab2)/2
	alt_ab2 <- vector_crossprod(ab - ab2, a - ab) - alt_ab
	vector_dihedral(alt_ab2, alt_ab, a, h)
}

f[f$hybrid %in% c("sp3", "ring"), "chi"] <-
	with(f[f$hybrid %in% c("sp3", "ring"),], alt_chi_dihedral_angle(
		cbind(ab2x, ab2y, ab2z), cbind(abx, aby, abz),
		cbind(ax, ay, az), cbind(hx, hy, hz)))

alt_sp3_cosBAH <- function(ab2, ab, a, h){
	alt_ab <- (ab + ab2)/2
	vector_dotprod(vector_normalize(a-alt_ab), vector_normalize(h-a))
}

f[f$hybrid == "sp3", "cosBAH"] <- with(f[f$hybrid == "sp3",], alt_sp3_cosBAH(
		cbind(ab2x, ab2y, ab2z), cbind(abx, aby, abz),
		cbind(ax, ay, az), cbind(hx, hy, hz)))

#equal area projection
f <- transform(f,
	capx = 2*sin(acos(cosBAH)/2)*cos(chi),
	capy = 2*sin(acos(cosBAH)/2)*sin(chi))
#capx_limits <- range(f$capx); capy_limits <- range(f$capy)
capx_limits <- c(-1.5,1.5); capy_limits <- capx_limits;

l_ply(levels(f$hybrid), function(hybrid){
	d_ply(sample_sources, .("sample_sources"), function(sample_source){
		ss_id <- sample_source$sample_source[1]
		plot_id = paste("hbond_chi_BAH_polar_density_", hybrid, "_by_don_chem_type_", ss_id, sep="")
		ggplot(data=subset(f, sample_source == ss_id & hybrid==hybrid)) + theme_bw() +
			polar_equal_area_grids_bw() +
			stat_bin2d(aes(x=capx, y=capy, fill=log(..density..)), binwidth=c(.06,.06)) +
			facet_grid(don_ss ~ acc_ss) +
			ggtitle(paste("Backbone-Backbone Hydrogen Bonds chi vs sinBAH Angles by Secondary Structure\nEqual Coordinate Projection   Sample Source: ", ss_id, sep="")) +
			scale_x_continuous('2*sin(BAH/2) * cos(CHI)', limits=capx_limits, breaks=c(-1, 0, 1)) +
			scale_y_continuous('2*sin(BAH/2) * sin(CHI)', limits=capy_limits, breaks=c(-1, 0, 1)) +
			scale_fill_gradientn('log(Normalized\nDensity)', colours=jet.colors(15))
		save_plots(self, plot_id, sample_source, output_dir, output_formats)
	})
})


})) # end FeaturesAnalysis
