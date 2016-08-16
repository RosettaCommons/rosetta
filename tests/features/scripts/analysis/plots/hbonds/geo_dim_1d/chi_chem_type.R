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
id = "chi_chem_type",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
	geom.chi,
	acc.HBChemType AS acc_chem_type,
	don.HBChemType AS don_chem_type,
  CASE acc.HBChemType
		WHEN 'hbacc_IMD' THEN 'ring' WHEN 'hbacc_IME' THEN 'ring'
		WHEN 'hbacc_AHX' THEN 'sp3'  WHEN 'hbacc_HXL' THEN 'sp3'
		WHEN 'hbacc_CXA' THEN 'sp2'  WHEN 'hbacc_CXL' THEN 'sp2'
		WHEN 'hbacc_PBA' THEN 'bb_sp2' END AS hybrid,
	acc_atoms.base2_x AS ab2x, acc_atoms.base2_y AS ab2y, acc_atoms.base2_z AS ab2z, -- acceptor base 2 atom
 	acc_atoms.base_x  AS abx,  acc_atoms.base_y  AS aby,  acc_atoms.base_z  AS abz,  -- acceptor base atom
	acc_atoms.atm_x   AS ax,   acc_atoms.atm_y   AS ay,   acc_atoms.atm_z   AS az,   -- acceptor atom
	don_atoms.atm_x   AS hx,   don_atoms.atm_y   AS hy,   don_atoms.atm_z   AS hz    -- hydrogen atom
FROM
	hbonds AS hb,
	hbond_geom_coords AS geom,
	hbond_sites AS don, hbond_sites AS acc,
	hbond_sites_pdb AS don_pdb, hbond_sites_pdb AS acc_pdb,
	hbond_site_atoms AS don_atoms, hbond_site_atoms AS acc_atoms
WHERE
	geom.struct_id = hb.struct_id AND geom.hbond_id = hb.hbond_id AND
	don.struct_id = hb.struct_id AND don.site_id = hb.don_id AND
	acc.struct_id = hb.struct_id AND acc.site_id = hb.acc_id AND
	don_pdb.struct_id = hb.struct_id AND don_pdb.site_id = hb.don_id AND
	don_pdb.heavy_atom_temperature < 30 AND
	acc_pdb.struct_id = hb.struct_id AND acc_pdb.site_id = hb.acc_id AND
	acc_pdb.heavy_atom_temperature < 30 AND
	don_atoms.struct_id = hb.struct_id AND don_atoms.site_id = hb.don_id AND
	acc_atoms.struct_id = hb.struct_id AND acc_atoms.site_id = hb.acc_id;"
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

#Convert from radians to degrees
f$chi <- (f$chi*180/pi)

# Order the plots better and give more descriptive labels
f$don_chem_type <- factor(f$don_chem_type,
	levels=c("hbdon_IMD", "hbdon_GDE", "hbdon_HXL", "hbdon_PBA",
		"hbdon_IME", "hbdon_GDH", "hbdon_AHX", "hbdon_IND",
		"hbdon_CXA", "hbdon_AMO"),
	labels = c("dIMD: h", "dGDE: r", "dHXL: s,t", "dPBA: bb",
		"dIME: h", "dGDH: r", "dAHX: y", "dIND: w",
		"dCXA: n,q", "dAMO: k"))

dens <- estimate_density_1d_wrap(
	f, c("sample_source", "hybrid", "don_chem_type"), "chi")

# This is a prototype for a multiplicative factor of the hbond energy
# to account for the chi angle preference.
cos_potential <- function(chi) -1*(cos( 2 * chi * pi/180 ) + 2 ) / 3

# Transform the potential energy to probability via the Boltzmann
# transformation where the partition function is over the chi angle
potential <- data.frame( x=seq(0,360,length.out=200))
potential$y <- exp(-1*cos_potential(potential$x))/(12.5803 *180/pi)

l_ply(levels(dens$hybrid), function(hybrid){
	plot_id = paste("hbond_chi_", hybrid, "_acceptor_by_don_chem_type", sep="")
	ggplot(data=dens[dens$hybrid==hybrid,]) + theme_bw() +
		geom_line(aes(x=x, y=y, colour=sample_source)) +
		geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
		facet_wrap( ~ don_chem_type) +
		ggtitle(paste("HBonds CHI Angle for ", hybrid, " acceptors and B-Factor < 30\n(Normalized for Equal Volume per Unit Distance)")) +
		scale_x_continuous('Acceptor Base -- Acceptor Torsion (degrees)', breaks=c(90,270)) +
		scale_y_continuous('Feature Density') +
		theme(legend.position=c(.58,.35)) +
		theme(legend.justification=c("left", "top"))
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
})


# Order the plots better and give more descriptive labels
f$don_chem_type <- factor(f$don_chem_type,
	levels = c("dIMD: h", "dIME: h", "dGDE: r", "dGDH: r",
		"dAHX: y", "dHXL: s,t", "dIND: w", "dAMO: k", "dCXA: n,q", "dPBA: bb"))

# Order the plots better and give more descriptive labels
f$acc_chem_type <- factor(f$acc_chem_type,
	levels=c("hbacc_IMD", "hbacc_IME", "hbacc_AHX", "hbacc_HXL",
		"hbacc_CXA", "hbacc_CXL", "hbacc_PBA"),
	labels = c("aIMD: h", "aIME: h", "aAHX: y", "aHXL: s,t",
		"aCXA: n,q", "aCXL: d,e", "aPBA: bb"))

plot_id = "hbond_chi_chem_type"
dens <- estimate_density_1d_wrap(
	f, c("sample_source", "acc_chem_type", "don_chem_type"), "chi")
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x, y=y, colour=sample_source)) +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	facet_grid(don_chem_type ~ acc_chem_type) +
	ggtitle("HBonds CHI Angle by Chemical Types, B-Factor < 30\n(normalized for equal volume per unit distance)") +
	scale_x_continuous('Acceptor Base -- Acceptor Torsion (degrees)', breaks=c(90,270)) +
	scale_y_continuous('Feature Density')

if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
