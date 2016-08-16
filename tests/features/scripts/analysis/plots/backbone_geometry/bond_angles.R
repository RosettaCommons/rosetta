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
id = "bond_angles",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("ResidueFeatures", "PdbDataFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


sele <-"
SELECT
	res.name3 AS res_type,
	dssp_code.label as dssp_label,
	b_ang.observed
FROM
	residues AS res,
--	residue_pdb_confidence AS res_conf,
	residue_secondary_structure AS ss,
	dssp_codes AS dssp_code,
	bond_intrares_angles AS b_ang
WHERE
--	res_conf.struct_id = res.struct_id AND res_conf.residue_number = res.resNum AND
--	res_conf.max_temperature < 30 AND
	ss.struct_id = res.struct_id AND ss.resNum == res.resNum AND
	dssp_code.code = ss.dssp AND
	b_ang.struct_id = res.struct_id AND b_ang.resNum = res.resNum AND
	b_ang.outAtm1Num = 1 AND b_ang.cenAtmNum = 2 AND b_ang.outAtm2Num = 3;"


f <- query_sample_sources(sample_sources, sele)
f$observed <- f$observed*180/pi

f <- f[as.character(f$res_type) != "MLY",]

f$is_helix <- as.character(f$dssp_label) == "H: a-Helix"

plot_parts <- list(
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)),
	geom_vline(x=111.2, size=1.5),
	scale_x_continuous("Bond Angle (degrees)", limit=c(95,125), breaks=c(95, 100, 105, 110, 115, 120, 125)),
	scale_y_continuous("Feature Density"),
	theme_bw())

plot_id <- "backbone_geometry_bond_angle_NCaC"
dens <- estimate_density_1d(f, c("sample_source"), "observed")
p <- ggplot(data=dens) + plot_parts +
	geom_line(aes(x, y, colour=sample_source), size=1.2) +
	ggtitle("Backbone N-Ca-C Bond Angle; B-Factor < 30")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

plot_id <- "backbone_geometry_bond_angle_NCaC_is_helix"
dens <- estimate_density_1d(f, c("sample_source", "is_helix"), "observed")
p <- ggplot(data=dens) + plot_parts +
	geom_line(aes(x, y, colour=sample_source), size=1.2) +
	facet_wrap(~is_helix, ncol=1) +
	ggtitle("Backbone N-Ca-C Bond Angle a-Helix vs Other; B-Factor < 30")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


plot_id <- "backbone_geometry_bond_angle_NCaC_by_res_type"
dens <- estimate_density_1d(f, c("sample_source", "res_type"), "observed")
p <- ggplot(data=dens) + plot_parts +
	geom_line(aes(x, y, colour=sample_source), size=1.2) +
	facet_wrap(~res_type) +
	ggtitle("Backbone N-Ca-C Bond Angle by Residue Type; B-Factor < 30")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

plot_id <- "backbone_geometry_bond_angle_NCaC_by_res_type_together"
dens <- estimate_density_1d(f, c("sample_source", "res_type"), "observed")
p <- ggplot(data=dens) + plot_parts +
	geom_line(aes(x, y, colour=res_type), size=1.2) +
	facet_wrap(~sample_source, ncol=1) +
	ggtitle("Backbone N-Ca-C Bond Angle by Residue Type; B-Factor < 30")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

plot_id <- "backbone_geometry_bond_angle_NCaC_by_ss"
dens <- estimate_density_1d(f, c("sample_source", "dssp_label"), "observed")
p <- ggplot(data=dens) + plot_parts +
	geom_line(aes(x, y, colour=sample_source), size=1.2) +
	facet_wrap(~dssp_label) +
	ggtitle("Backbone N-Ca-C Bond Angle by DSSP; B-Factor < 30")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

plot_id <- "backbone_geometry_bond_angle_NCaC_by_ss_together"
dens <- estimate_density_1d(f, c("sample_source", "dssp_label"), "observed")
p <- ggplot(data=dens) + plot_parts +
	geom_line(aes(x, y, colour=dssp_label), size=1.2) +
	facet_wrap(~sample_source, ncol=1) +
	ggtitle("Backbone N-Ca-C Bond Angle by DSSP; B-Factor < 30")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


plot_id <- "backbone_geometry_bond_angle_NCaC_by_res_type_ss"
dens <- estimate_density_1d(f[f$res_type != "PRO",], c("sample_source", "res_type", "dssp_label"), "observed")
p <- ggplot(data=dens) + plot_parts +
	geom_line(aes(x, y, colour=sample_source), size=1) +
	facet_grid(res_type ~ dssp_label) +
	ggtitle("Backbone N-Ca-C Bond Angle by ResType and DSSP; B-Factor < 30")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)





})) # end FeaturesAnalysis
