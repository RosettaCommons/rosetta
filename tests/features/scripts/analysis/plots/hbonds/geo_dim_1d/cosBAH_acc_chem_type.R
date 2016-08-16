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
id = "cosBAH_acc_chem_type",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
	acc_atoms.base_x AS bx, acc_atoms.base_y AS by, acc_atoms.base_z AS bz,
	acc_atoms.atm_x  AS ax, acc_atoms.atm_y  AS ay, acc_atoms.atm_z  AS az,
	don_atoms.atm_x  AS hx, don_atoms.atm_y  AS hy, don_atoms.atm_z  AS hz,
	acc.HBChemType AS acc_chem_type,
	don.HBChemType AS don_chem_type,
	abs(don.resNum - acc.resNum) AS seq_sep,
	don_pdb.heavy_atom_temperature AS don_temp,
	acc_pdb.heavy_atom_temperature AS acc_temp
FROM
	hbonds AS hb,
	hbond_geom_coords AS geom,
	hbond_sites AS don, hbond_sites AS acc,
	hbond_site_atoms AS don_atoms, hbond_site_atoms AS acc_atoms,
	hbond_sites_pdb AS don_pdb, hbond_sites_pdb AS acc_pdb
WHERE
	geom.struct_id = hb.struct_id AND geom.hbond_id = hb.hbond_id AND
	don.struct_id = hb.struct_id AND don.site_id = hb.don_id AND
	acc.struct_id = hb.struct_id AND acc.site_id = hb.acc_id AND
	don_atoms.struct_id = hb.struct_id AND don_atoms.site_id = hb.don_id AND
	acc_atoms.struct_id = hb.struct_id AND acc_atoms.site_id = hb.acc_id AND
	don_pdb.struct_id = hb.struct_id AND don_pdb.site_id = hb.don_id AND
	acc_pdb.struct_id = hb.struct_id AND acc_pdb.site_id = hb.acc_id;"

f <- query_sample_sources(sample_sources, sele)

f$acc_chem_type_label <- factor(f$acc_chem_type,
	levels =
		c("hbacc_CXA", "hbacc_AHX", "hbacc_IMD",
			"hbacc_CXL", "hbacc_HXL", "hbacc_IME",
			"hbacc_PBA"),
	labels =
		c("aCXA: n,q", "aAHX: y",   "aIMD: h",
			"aCXL: d,e", "aHXL: s,t", "aIME: h",
			"aPBA: bb"))

f <- transform(f,
	cosBAH = vector_dotprod(
		vector_normalize(cbind(ax-bx, ay-by, az-bz)),
		vector_normalize(cbind(hx-ax, hy-ay, hz-az))))

dens <- estimate_density_1d(
	f, c("sample_source", "acc_chem_type_label" ), "cosBAH")

plot_id = "cosBAH_acc_chem_type"
ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=acos(x)*180/pi, y=-log(y), colour=sample_source)) +
	geom_indicator(aes(colour=sample_source, indicator=counts, group=sample_source)) +
	facet_wrap( ~ acc_chem_type_label) +
	ggtitle("Hydrogen Bonds BAH Angle by Acceptor Chemical Type\n(normalized for equal volume per unit distance)") +
	scale_x_continuous(paste('Base -- Acceptor -- Hydrogen (degrees)')) +
	scale_y_continuous("-log(FeatureDensity)", limits=c(-2.3,6)) +
	theme(legend.position=c(.58,.35)) +
	theme(legend.justification=c("left", "top"))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

filt_f <- f[
	f$don_temp < 30 &
	f$acc_temp < 30 &
	f$don_chem_type != "hbdon_PBA" &
	f$seq_sep > 5, ]

dens <- estimate_density_1d(
	filt_f, c("sample_source", "acc_chem_type_label" ), "cosBAH")

plot_id = "cosBAH_acc_chem_type_filtered"
ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=acos(x)*180/pi, y=-log(y), colour=sample_source)) +
	geom_indicator(aes(colour=sample_source, indicator=counts, group=sample_source)) +
	facet_wrap( ~ acc_chem_type_label) +
	ggtitle("HBonds: Sidechain Donor, BFact < 30, SeqSep > 5   BAH Angle by Acceptor Chemical Type\n(normalized for equal volume per unit distance)") +
	scale_x_continuous(paste('Base -- Acceptor -- Hydrogen (degrees)')) +
	scale_y_continuous("-log(FeatureDensity)", limits=c(-2.3,6)) +
	theme(legend.position=c(.58,.35)) +
	theme(legend.justification=c("left", "top"))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
