# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

check_setup()
feature_analyses <- c(feature_analyses, new("FeaturesAnalysis",
id = "AHdist_mode_differences",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
	geom.AHdist,
	acc.HBChemType AS acc_chem_type,
	don.HBChemType AS don_chem_type
FROM
	hbond_geom_coords AS geom,
	hbonds AS hb,
	hbond_sites_pdb AS don_pdb,
	hbond_sites_pdb AS acc_pdb,
	hbond_sites AS don,
	hbond_sites AS acc
WHERE
	hb.struct_id = geom.struct_id AND
	hb.hbond_id = geom.hbond_id AND
	hb.struct_id = don.struct_id AND
	hb.don_id = don.site_id AND
	hb.struct_id = acc.struct_id AND
	hb.acc_id = acc.site_id AND
	don_pdb.struct_id = hb.struct_id AND don_pdb.site_id = hb.don_id AND
	don_pdb.heavy_atom_temperature < 30 AND
	acc_pdb.struct_id = hb.struct_id AND acc_pdb.site_id = hb.acc_id AND
	acc_pdb.heavy_atom_temperature < 30;";

f <- query_sample_sources(sample_sources, sele)

print(summary(f))

all_modes <- estimate_primary_modes_1d(f,
  c("sample_source"), "AHdist")
print(all_modes)

acc_modes <- estimate_primary_modes_1d(f,
  c("sample_source", "acc_chem_type"), "AHdist")
print(acc_modes)

don_modes <- estimate_primary_modes_1d(f,
  c("sample_source", "don_chem_type"), "AHdist")
print(don_modes)

each_modes <- estimate_primary_modes_1d(f,
  c("sample_source", "don_chem_type", "acc_chem_type"), "AHdist")
print(each_modes)

z <- data.frame(
	AHdist=each_modes$mode_location,
	acc=each_modes$acc_chem_type,
	don=each_modes$don_chem_type,
	each_counts=each_modes$counts)

z <- merge(z,
	data.frame(
		acc=acc_modes$acc_chem_type,
		acc_AHdist=acc_modes$mode_location,
		acc_counts=acc_modes$counts))

z <- merge(z,
	data.frame(
		don=don_modes$don_chem_type,
		don_AHdist=don_modes$mode_location,
		don_counts=don_modes$counts))

z <- merge(z,
	data.frame(
		all_AHdist=all_modes$mode_location))

print(summary(z))

})) # end FeaturesAnalysis
