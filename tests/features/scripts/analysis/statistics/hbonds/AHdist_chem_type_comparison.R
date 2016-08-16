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
id = "AHdist_chem_type_comparison",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){



source("scripts/analysis/plots/hbonds/hbond_geo_dim_scales.R")

sele <-"
SELECT
	geom.AHdist,
	acc.HBChemType AS acc_chem_type,
	don.HBChemType AS don_chem_type,
	CASE acc.HBChemType
			WHEN 'hbacc_IMD' THEN 'ring' WHEN 'hbacc_IME' THEN 'ring'
			WHEN 'hbacc_AHX' THEN 'sp3'  WHEN 'hbacc_HXL' THEN 'sp3'
			WHEN 'hbacc_CXA' THEN 'sp2'  WHEN 'hbacc_CXL' THEN 'sp2'
			WHEN 'hbacc_PBA' THEN 'sp2'  END AS acc_hybrid
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

f$don_chem_type_name <- don_chem_type_name_linear(f$don_chem_type)
f$acc_chem_type_name <- acc_chem_type_name_linear(f$acc_chem_type)

f <- na.omit(f, method="r")

tests <- c("kolmogorov_smirnov_test", "histogram_kl_divergence", "primary_modes_radial_3d_norm_diff")

comp_stats <- comparison_statistics(
	sample_sources, f, c(), "AHdist", tests)
table_id <- "AHdist_chem_type_comparison"
table_title <- "H-Bond Distance Distribution Comparison, B-Factor < 30"
save_tables(self,
	comp_stats, table_id,
	sample_sources, output_dir, output_formats,
	caption=table_title, caption.placement="top")

comp_stats <- comparison_statistics(
	sample_sources, f, c("don_chem_type_name"), "AHdist", tests)
table_id <- paste("AHdist_chem_type_comparison", "by_don_chem_type", sep="_")
table_title <- "H-Bond Distance by Donor Chemical Type\nDistribution Comparison, B-Factor < 30"
save_tables(self,
	comp_stats, table_id, sample_sources, output_dir, output_formats,
	caption=table_title, caption.placement="top")

comp_stats <- comparison_statistics(
	sample_sources, f, c("acc_chem_type_name"), "AHdist", tests)
table_id <- paste("AHdist_chem_type_comparison", "by_acc_chem_type", sep="_")
table_title <- "H-Bond Distance by Acceptor Chemical Type\nDistribution Comparison, B-Factor < 30"
save_tables(self,
	comp_stats, table_id, sample_sources, output_dir, output_formats,
	caption=table_title, caption.placement="top")

comp_stats <- comparison_statistics(
	sample_sources, f, c("acc_hybrid"), "AHdist", tests)
table_id <- paste("AHdist_chem_type_comparison", "by_acc_hybrid", sep="_")
table_title <- "H-Bond Distance by Acceptor Hybrid\nDistribution Comparison, B-Factor < 30"
save_tables(self,
	comp_stats, table_id, sample_sources, output_dir, output_formats,
	caption=table_title, caption.placement="top")

comp_stats <- comparison_statistics(
	sample_sources, f, c("don_chem_type_name", "acc_chem_type_name"), "AHdist", tests)
table_id <- paste("AHdist_chem_type_comparison", "by_don_chem_type_acc_chem_type", sep="_")
table_title <- "H-Bond Distance by Donor and Acceptor Chemical Types\nDistribution Comparison, B-Factor < 30"
save_tables(self,
	comp_stats, table_id, sample_sources, output_dir, output_formats,
	caption=table_title, caption.placement="top")

})) # end FeaturesAnalysis

