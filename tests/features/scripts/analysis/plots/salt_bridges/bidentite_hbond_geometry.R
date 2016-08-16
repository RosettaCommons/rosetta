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
id = "bidentite_hbond_geometry",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("StructureFeatures", "ResidueFeatures", "SaltBridgeFeatures", "PdbDataFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

problem_region <- data.frame(xmin=-121, xmax=-115)

sele <-paste("
CREATE TABLE IF NOT EXISTS CXL_ARG_salt_bridges AS SELECT
	struct.tag AS tag,
	sb.struct_id AS struct_id,
	sb.don_resNum AS don_residue_number,
	sb.acc_id AS acc_id,
	acc.resNum AS acc_residue_number,
	sb.psi * 180/3.14159 AS psi
FROM
	structures AS struct,
	salt_bridges AS sb,
	residues AS don,
	hbond_sites AS acc
WHERE
	sb.struct_id = struct.struct_id AND
	don.struct_id = sb.struct_id AND
	don.resNum = sb.don_resNum AND
	don.name3 = 'ARG' AND
	acc.struct_id = sb.struct_id AND
	acc.site_id = sb.acc_id AND
	acc.HBChemType = 'hbacc_CXL' AND
	ABS(acc.resNum - don.resNum) > 4;
	
CREATE INDEX IF NOT EXISTS hbonds_struct_id_acc_id ON
	hbonds (struct_id, acc_id);

SELECT
	sb.tag,
	sb.struct_id AS struct_id,
	sb.don_residue_number AS don_residue_number,
	sb.acc_residue_number AS acc_residue_number,
	sb.psi,
	hb1_geo.AHdist AS hb1_AHdist, hb2_geo.AHdist AS hb2_AHdist,
	hb1_geo.cosAHD AS hb1_cosAHD, hb2_geo.cosAHD AS hb2_cosAHD,
	hb1_geo.cosBAH AS hb1_cosBAH, hb2_geo.cosBAH AS hb2_cosBAH,
	hb1_geo.chi    AS hb1_chi,    hb2_geo.chi    AS hb2_chi,
	CASE WHEN
		sb.psi < ", problem_region$xmax[1], " AND
		", problem_region$xmin, " < sb.psi THEN 1
	ELSE 0 END AS in_bifurcated_region
FROM
	CXL_ARG_salt_bridges AS sb,
	hbonds AS hb1,
	hbonds AS hb2,
	hbond_sites AS don1,
	hbond_sites AS don2,
	hbond_geom_coords AS hb1_geo,
	hbond_geom_coords AS hb2_geo
WHERE
	-- hb1 and hb2 are hydrogen bonds that are bifurcated at the acceptor
	-- between an ASP or GLU hbond site to ARG with seqsep > 4
	-- the ASP -> ARG forms a salt bridge which is evaluated
	hb1.struct_id = sb.struct_id AND hb1.acc_id = sb.acc_id AND
	hb2.struct_id = sb.struct_id AND hb2.acc_id = sb.acc_id AND
	don1.struct_id = hb1.struct_id AND don1.site_id = hb1.don_id AND
	don2.struct_id = hb2.struct_id AND don2.site_id = hb2.don_id AND
	don1.site_id != don2.site_id AND
	hb1_geo.struct_id = hb1.struct_id AND
	hb1_geo.hbond_id = hb1.hbond_id AND
	hb2_geo.struct_id = hb2.struct_id AND
	hb2_geo.hbond_id = hb2.hbond_id", sep="")
f_bifurcated <- query_sample_sources(sample_sources, sele)
f_bifurcated <- na.omit(f_bifurcated, method="r")
f_bifurcated$salt_bridge_type <- factor("bifurcated")
f_bifurcated$in_bifurcated_region <- factor(f_bifurcated$in_bifurcated_region)
f_bifurcated$hb1_AHD <- acos(f_bifurcated$hb1_cosAHD)*180/pi
f_bifurcated$hb1_BAH <- acos(f_bifurcated$hb1_cosBAH)*180/pi



sele <- paste("
SELECT
	sb1.tag,
	sb1.struct_id AS struct_id,
	sb1.don_residue_number AS don_residue_number,
	sb1.acc_residue_number AS acc_residue_number,
	sb1.psi AS psi,
	sb2.psi AS sb2_psi,
	hb1_geo.AHdist AS hb1_AHdist, hb2_geo.AHdist AS hb2_AHdist,
	hb1_geo.cosAHD AS hb1_cosAHD, hb2_geo.cosAHD AS hb2_cosAHD,
	hb1_geo.cosBAH AS hb1_cosBAH, hb2_geo.cosBAH AS hb2_cosBAH,
	hb1_geo.chi    AS hb1_chi,    hb2_geo.chi    AS hb2_chi,
	CASE WHEN
		sb1.psi < ", problem_region$xmax[1], " AND
		", problem_region$xmin, " < sb1.psi THEN 1
	ELSE 0 END AS in_bifurcated_region,
	CASE WHEN
		sb2.psi < ", problem_region$xmax[1], " AND
		", problem_region$xmin, " < sb2.psi THEN 1
	ELSE 0 END AS sb2_in_bifurcated_region
FROM
	CXL_ARG_salt_bridges AS sb1, CXL_ARG_salt_bridges AS sb2,
	hbonds AS hb1, hbonds AS hb2,
	hbond_sites AS don1, hbond_sites AS don2,
	hbond_geom_coords AS hb1_geo, hbond_geom_coords AS hb2_geo
WHERE
	sb2.struct_id = sb1.struct_id AND
	sb1.acc_id != sb2.acc_id AND
	sb2.don_residue_number = sb1.don_residue_number AND
	sb2.acc_residue_number = sb1.acc_residue_number AND
	hb1.struct_id = sb1.struct_id AND hb1.acc_id = sb1.acc_id AND
	hb2.struct_id = sb2.struct_id AND hb2.acc_id = sb2.acc_id AND
	don1.struct_id = hb1.struct_id AND don1.site_id = hb1.don_id AND
	don2.struct_id = hb2.struct_id AND don2.site_id = hb2.don_id AND
	don1.resNum = sb1.don_residue_number AND
	don2.resNum = sb2.don_residue_number AND
	don1.site_id != don2.site_id AND
	hb1_geo.struct_id = hb1.struct_id AND
	hb1_geo.hbond_id = hb1.hbond_id AND
	hb2_geo.struct_id = hb2.struct_id AND
	hb2_geo.hbond_id = hb2.hbond_id;", sep="")
f_bidentite <- query_sample_sources(sample_sources, sele)
f_bidentite$salt_bridge_type <- factor("bidentite")
f_bidentite <- na.omit(f_bidentite, method="r")
f_bidentite$in_bifurcated_region <- factor(f_bidentite$in_bifurcated_region)
f_bidentite$hb1_AHD <- 180-acos(f_bidentite$hb1_cosAHD)*180/pi
f_bidentite$hb1_BAH <- 180-acos(f_bidentite$hb1_cosBAH)*180/pi

plot_id <- "salt_bridge_bifurcated_ASP_GLU_to_ARG_psi_vs_AHdist"
p <- ggplot(data=f_bifurcated) + theme_bw() +
	geom_point(aes(x=psi, y=hb1_AHdist), size=.8) +
	ggtitle("ASP/GLU -> ARG Bifurcated Salt Bridges; SeqSep > 4") +
	scale_x_continuous("Angle around Donor (degrees)") +
	scale_y_continuous("A-H Distance") +
	facet_wrap(~sample_source, ncol=1)
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

plot_id <- "salt_bridge_bidentite_ASP_GLU_to_ARG_psi_vs_AHdist"
p <- ggplot(data=f_bidentite) + theme_bw() +
	geom_point(aes(x=psi, y=hb1_AHdist), size=.8) +
	ggtitle("ASP/GLU -> ARG Bidentite Salt Bridges; SeqSep > 4") +
	scale_x_continuous("Angle around Donor (degrees)") +
	scale_y_continuous("A-H Distance") +
	facet_wrap(~sample_source, ncol=1)
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


plot_id <- "salt_bridge_bifurcated_ASP_GLU_to_ARG_psi_vs_AHD"
p <- ggplot(data=f_bifurcated) + theme_bw() +
	geom_point(aes(x=psi, y=acos(hb1_cosAHD)*180/pi), size=.8) +
	ggtitle("ASP/GLU -> ARG Bifurcated Salt Bridges; SeqSep > 4") +
	scale_x_continuous("Angle around Donor (degrees)") +
	scale_y_continuous("A-H-D Angle (degrees)") +
	facet_wrap(~sample_source, ncol=1)
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

plot_id <- "salt_bridge_bidentite_ASP_GLU_to_ARG_psi_vs_AHD"
p <- ggplot(data=f_bidentite) + theme_bw() +
	geom_point(aes(x=psi, y=acos(hb1_cosAHD)*180/pi), size=.8) +
	ggtitle("ASP/GLU -> ARG Bidentite Salt Bridges; SeqSep > 4") +
	scale_x_continuous("Angle around Donor (degrees)") +
	scale_y_continuous("A-H-D Angle (degrees)") +
	facet_wrap(~sample_source, ncol=1)
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


plot_id <- "salt_bridge_bifurcated_ASP_GLU_to_ARG_psi_vs_BAH"
p <- ggplot(data=f_bifurcated) + theme_bw() +
	geom_point(aes(x=psi, y=acos(hb1_cosBAH)*180/pi), size=.8) +
	ggtitle("ASP/GLU -> ARG Bifurcated Salt Bridges; SeqSep > 4") +
	scale_x_continuous("Angle around Donor (degrees)") +
	scale_y_continuous("B-A-H Angle (degrees)") +
	facet_wrap(~sample_source, ncol=1)
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

plot_id <- "salt_bridge_bidentite_ASP_GLU_to_ARG_psi_vs_BAH"
p <- ggplot(data=f_bidentite) + theme_bw() +
	geom_point(aes(x=psi, y=acos(hb1_cosBAH)*180/pi), size=.8) +
	ggtitle("ASP/GLU -> ARG Bidentite Salt Bridges; SeqSep > 4") +
	scale_x_continuous("Angle around Donor (degrees)") +
	scale_y_continuous("B-A-H Angle (degrees)") +
	facet_wrap(~sample_source, ncol=1)
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


plot_id <- "salt_bridge_bifurcated_ASP_GLU_to_ARG_psi_in_bifurcated_region_AHD_AHdist"
p <- ggplot(data=f_bifurcated) + theme_bw() +
	geom_point(
		aes(x=acos(hb1_cosAHD)*180/pi, y=hb1_AHdist, colour=in_bifurcated_region),
		size=.8) +
	ggtitle("ASP/GLU -> ARG Bifurcated Salt Bridges; SeqSep > 4") +
	scale_x_continuous("A-H-D Angle (degrees)") +
	scale_y_continuous("A-H Distance") +
	facet_wrap(~sample_source, ncol=1)
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

plot_id <- "salt_bridge_bidentite_ASP_GLU_to_ARG_psi_in_bifurcated_region_AHD_AHdist"
p <- ggplot(data=f_bidentite) + theme_bw() +
	geom_point(
		aes(x=acos(hb1_cosAHD)*180/pi, y=hb1_AHdist, colour=in_bifurcated_region),
		size=.8) +
	ggtitle("ASP/GLU -> ARG Bidentite Salt Bridges; SeqSep > 4") +
	scale_x_continuous("A-H-D Angle (degrees)") +
	scale_y_continuous("A-H Distance") +
	facet_wrap(~sample_source, ncol=1)
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


#equal area projection
f_bifurcated <- transform(f_bifurcated,
	capx = 2*sin(acos(hb1_cosBAH)/2)*cos(hb1_chi),
	capy = 2*sin(acos(hb1_cosBAH)/2)*sin(hb1_chi))

capx_limits <- c(-1.5,1.5)
capy_limits <- c(-1.5,1.5)

max_BAH_angle = 110 # in degrees
abs_cap_limit <- 2*sin((max_BAH_angle*pi/180)/2)
capx_limits <- c(-abs_cap_limit, abs_cap_limit)
capy_limits <- c(-abs_cap_limit, abs_cap_limit)

plot_parts <- list(
	theme_bw(),
	polar_equal_area_grids_bw(),
	geom_point(aes(x=capx, y=capy)),
	coord_equal(ratio=1),
	facet_grid(sample_source ~ in_bifurcated_region),
	scale_x_continuous('', limits=capx_limits, breaks=c()),
	scale_y_continuous('', limits=capy_limits, breaks=c()),
	theme(
		axis.text.x=theme_blank(),
		axis.text.y=theme_blank(),
		axis.title.x=theme_blank(),
		axis.title.y=theme_blank(),
		axis.ticks.x = theme_blank(),
		axis.ticks.y = theme_blank()))

plot_id <- "salt_bridge_bifurcated_ASP_GLU_to_ARG__BAH_vs_chi_by_in_bifurcated_sample_source"
ggplot(data=f_bifurcated) + plot_parts +
	ggtitle("Acceptor Hydrogen Bonds seq_sep > 4: chi vs BAH bifurcated salt bridges") +
save_plots(self, plot_id, sample_sources, output_dir, output_formats)



plot_id <- "salt_bridge_bifurcated_ASP_GLU_to_ARG__AHD"
dens <- estimate_density_1d_reflect_boundary(
	f_bifurcated,
	c("sample_source", "in_bifurcated_region"),
	"hb1_AHD",
	reflect_left=TRUE,
	right_boundary=0,
	conical_3d_normalization,
	adjust=.5)

dens <- estimate_density_1d(
	f_bifurcated,
	c("sample_source", "in_bifurcated_region"),
	"hb1_AHD",
	conical_3d_normalization,
	adjust=.8)

p <- ggplot(data=dens) + theme_bw() +
	ggtitle("Bifurcated ASP/GLU -> ARG Salt Bridge AHD; SeqSep > 4") +
	geom_line(aes(x=x, y=y, colour=sample_source)) +
	facet_wrap(~in_bifurcated_region) +
	scale_y_continuous("FeatureDensity") +          
	scale_x_continuous("A-H-D Angle (degrees)")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


plot_id <- "salt_bridge_bidentite_ASP_GLU_to_ARG__AHD"
dens <- estimate_density_1d_reflect_boundary(
	f_bidentite,
	c("sample_source", "in_bifurcated_region"),
	"hb1_AHD",
	reflect_left=TRUE,
	right_boundary=0,
	conical_3d_normalization,
	adjust=.5)

p <- ggplot(data=dens) + plot_parts +
	ggtitle("Bidentite ASP/GLU -> ARG Salt Bridge AHD; SeqSep > 4") +
	geom_line(aes(x=x, y=y, colour=sample_source)) +
	facet_wrap(~in_bifurcated_region) +
	scale_y_continuous("FeatureDensity") +          
	scale_x_continuous("A-H-D Angle (degrees)")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)



})) # end FeaturesAnalysis
