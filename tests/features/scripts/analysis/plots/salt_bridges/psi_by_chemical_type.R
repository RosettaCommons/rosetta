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
id = "salt_bridge_geo_dim_1d",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("StructureFeatures", "ResidueFeatures", "SaltBridgeFeatures", "PdbDataFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


sele <-"
SELECT
	struct.tag, don.resNum AS don_resNum, acc.resNum AS acc_resNum,
	sb.psi,
	don.name3 AS don_res_type, acc.HBChemType AS acc_chem_type,
	CASE don.resNum - acc.resNum
		WHEN -1 THEN '-1' WHEN -2 THEN '-2' WHEN -3 THEN '-3' WHEN -4 THEN '-4'
		WHEN 1 THEN '1' WHEN 2 THEN '2' WHEN 3 THEN '3' WHEN 4 THEN '4'
		ELSE 'long' END AS seq_sep
FROM
	structures AS struct,
	salt_bridges AS sb,
	hbond_sites AS acc, residues as don
--	residue_pdb_confidence AS don_con, hbond_sites_pdb AS acc_con
WHERE
	sb.struct_id = struct.struct_id AND
	don.struct_id = sb.struct_id AND don.resNum = sb.don_resNum AND
	acc.struct_id = sb.struct_id AND acc.site_id = sb.acc_id AND
--	don_con.struct_id = sb.struct_id AND don_con.residue_number = sb.don_resNum AND
--	acc_con.struct_id = sb.struct_id AND acc_con.site_id = sb.acc_id AND
--	don_con.max_sc_temperature < 30 AND acc_con.heavy_atom_temperature < 30 AND
	sb.rho > 2;"

f <- query_sample_sources(sample_sources, sele)
f <- na.omit(f, method="r")

# give more descriptive plot labels
f$acc_chem_type <- factor(f$acc_chem_type,
	levels = c("hbacc_CXA", "hbacc_CXL"),
	labels = c("aCXA: n,q", "aCXL: d,e"))

f$acc_type <- with(f, interaction(acc_chem_type, orbital, sep="/"))

f$seq_sep <- factor(f$seq_sep,
	levels = c("-4", "-3", "-2", "-1", "1", "2", "3", "4", "long"),
	labels = c("-4", "-3", "-2", "-1", "1", "2", "3", "4", "long"))

plot_parts <- list(
	theme_bw(),
	geom_line(aes(x=x, y=y, colour=sample_source)),
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)),
	scale_y_continuous("FeatureDensity"))


plot_id <- "salt_bridge_psi"
dens <- estimate_density_1d_wrap(f, c("sample_source"), "psi", xlim=c(-180, 180), adjust=.2)
p <- ggplot(data=dens) + plot_parts +
	ggtitle("Salt Bridge PSI, B-Factor < 30") +
	scale_x_continuous("Angle Around Donor")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


plot_id <- "salt_bridge_psi_by_acc_type_don_res_type"
dens <- estimate_density_1d_wrap(f, c("sample_source", "don_res_type", "acc_type"), "psi", xlim=c(-180, 180), adjust=.2)
p <- ggplot(data=dens) + plot_parts +
	ggtitle("Salt Bridge PSI, B-Factor < 30") +
	facet_grid( acc_type ~ don_res_type ) +
	scale_x_continuous("Angle Around Donor")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


plot_id <- "salt_bridge_psi_by_orbital_don_res_type"
dens <- estimate_density_1d_wrap(f, c("sample_source", "don_res_type", "orbital"), "psi", xlim=c(-180, 180), adjust=.2)
p <- ggplot(data=dens) + plot_parts +
	ggtitle("Salt Bridge PSI, B-Factor < 30") +
	facet_grid( orbital ~ don_res_type ) +
	scale_x_continuous("Angle Around Donor")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


plot_id <- "salt_bridge_psi_by_acc_chem_type_don_res_type"
dens <- estimate_density_1d_wrap(f, c("sample_source", "don_res_type", "acc_chem_type"), "psi", xlim=c(-180, 180), adjust=.2)
p <- ggplot(data=dens) + plot_parts +
	ggtitle("Salt Bridge PSI, B-Factor < 30") +
	facet_grid( acc_chem_type ~ don_res_type ) +
	scale_x_continuous("Angle Around Donor")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

plot_id <- "salt_bridge_psi_by_orbital_acc_chem_type"
dens <- estimate_density_1d_wrap(f, c("sample_source", "orbital", "acc_chem_type"), "psi", xlim=c(-180, 180), adjust=.2)
p <- ggplot(data=dens) + plot_parts +
	ggtitle("Salt Bridge PSI, B-Factor < 30") +
	facet_grid( orbital ~ acc_chem_type ) +
	scale_x_continuous("Angle Around Donor")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

plot_id <- "salt_bridge_psi_by_seq_sep_don_res_type"
dens <- estimate_density_1d(f, c("sample_source", "seq_sep", "don_res_type"), "psi", xlim=c(-180, 180), adjust=.2)
p <- ggplot(data=dens) + plot_parts +
	ggtitle("Salt Bridge PSI, B-Factor < 30") +
	facet_grid( seq_sep ~ don_res_type ) +
	scale_x_continuous("Angle Around Donor")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)



}))
