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
id = "alpha_helices",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){



sele <-"
CREATE TEMPORARY TABLE seq_sep_4_hbs AS SELECT
	hb.struct_id, don.resNum AS don_resNum, acc.resNum AS acc_resNum,
	geom.AHdist, geom.cosBAH, geom.cosAHD, geom.chi
FROM
	hbonds AS hb,
	hbond_geom_coords AS geom,
	hbond_sites AS don, hbond_sites AS acc,
	hbond_sites_pdb AS don_pdb, hbond_sites_pdb AS acc_pdb
WHERE
	geom.struct_id = hb.struct_id AND geom.hbond_id = hb.hbond_id AND
	don.struct_id = hb.struct_id AND don.site_id = hb.don_id AND
	don.HBChemType = 'hbdon_PBA' AND
	acc.struct_id = hb.struct_id AND acc.site_id = hb.acc_id AND
	acc.HBChemType = 'hbacc_PBA' AND
	don_pdb.struct_id = hb.struct_id AND don_pdb.site_id = hb.don_id AND
	don_pdb.heavy_atom_temperature < 30 AND
	acc_pdb.struct_id = hb.struct_id AND acc_pdb.site_id = hb.acc_id AND
	acc_pdb.heavy_atom_temperature < 30 AND
	don.resNum - acc.resNum = 4;

SELECT
	hb.AHdist, hb.cosBAH, hb.cosAHD, hb.chi,
	hb.struct_id, hb.don_resNum as don_resNum, hb.acc_resNum as acc_resNum,
	hb_m1.don_resNum as m1_don_resNum, hb_m1.acc_resNum as m1_acc_resNum,
	hb_p1.don_resNum as p1_don_resNum, hb_p1.acc_resNum as p1_acc_resNum,
	CASE WHEN hb_m1.struct_id IS NULL THEN 0 ELSE 1 END AS hb_m1,
	CASE WHEN hb_p1.struct_id IS NULL THEN 0 ELSE 1 END AS hb_p1
FROM
	seq_sep_4_hbs AS hb
	LEFT JOIN seq_sep_4_hbs AS hb_m1 ON
		hb.struct_id = hb_m1.struct_id AND
		hb.don_resNum - 1 = hb_m1.don_resNum AND
		hb.acc_resNum - 1 = hb_m1.acc_resNum
	LEFT JOIN seq_sep_4_hbs AS hb_p1 ON
		hb.struct_id = hb_p1.struct_id AND
		hb.don_resNum + 1 = hb_p1.don_resNum AND
		hb.acc_resNum + 1 = hb_p1.acc_resNum;"

f <- query_sample_sources(sample_sources, sele)

f$alpha_type <- factor(1 + f$hb_m1 + 2*f$hb_p1,
	levels = c(1, 2, 3, 4),
	labels = c("alpha_turn", "alpha_helix_begin", "alpha_helix_end", "alpha_helix_middle"))


dens <- estimate_density_1d(
	f, c("sample_source", "alpha_type"),
	"AHdist", weight_fun = radial_3d_normalization, adjust=.5)

plot_id <- "hbond_backbone_backbone_alpha_type_AHdist"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x, y=y, colour=alpha_type)) +
	geom_indicator(aes(indicator=counts, colour=alpha_type, group=alpha_type)) +
	facet_wrap( ~ sample_source) +
	ggtitle("Backbone-Backbone H-Bonds A-H Distance by Alpha Type\nB-Factor < 30 normalized for equal weight per unit distance") +
	scale_y_continuous("FeatureDensity", limits=c(0,6), breaks=c(1,3,5)) +
	scale_x_continuous(expression(paste('Acceptor -- Proton Distance (', ring(A), ')')), limits=c(1.4,2.7), breaks=c(1.6, 1.9, 2.2, 2.6))

if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)




dens <- estimate_density_1d_reflect_boundary(
	f, c("sample_source", "alpha_type"),
	"cosAHD", reflect_right=T, right_boundary=1, adjust=.5)

plot_id <- "hbond_backbone_backbone_alpha_type_cosAHD"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=180-acos(x)*180/pi, y=y, colour=alpha_type)) +
	geom_indicator(aes(indicator=counts, colour=alpha_type, group=alpha_type)) +
	facet_wrap( ~ sample_source) +
	ggtitle("Backbone-Backbone H-Bonds AHD Angle by Alpha Type\nB-Factor < 30 normalized for equal weight per unit distance") +
	scale_y_continuous("FeatureDensity", limits=c(0,20), breaks=c(0,5,10,15)) +
	scale_x_continuous("Acceptor -- Hydrogen -- Donor (degrees)", trans="reverse")

if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)



dens <- estimate_density_1d(
	f, c("sample_source", "alpha_type"), "cosBAH", adjust=.5)

plot_id <- "hbond_backbone_backbone_alpha_type_cosBAH"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=acos(x)*180/pi, y=y, colour=alpha_type)) +
	geom_indicator(aes(indicator=counts, colour=alpha_type, group=alpha_type)) +
	facet_wrap( ~ sample_source) +
	ggtitle("Backbone-Backbone H-Bonds BAH Angle by Alpha Type\nB-Factor < 30 normalized for equal weight per unit distance") +
	scale_x_continuous(paste('Base -- Acceptor -- Hydrogen (degrees)')) +
	scale_y_continuous("FeatureDensity")

if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)


f$chi_deg <- f$chi*180/pi
dens <- estimate_density_1d_wrap(
	f, c("sample_source", "alpha_type"), "chi_deg", adjust=.5)

plot_id <- "hbond_backbone_backbone_alpha_type_chi"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x, y=y, colour=alpha_type)) +
	geom_indicator(aes(indicator=counts, colour=alpha_type, group=alpha_type)) +
	facet_wrap( ~ sample_source) +
	ggtitle("Backbone-Backbone H-Bonds CHI Angle by Alpha Type\nB-Factor < 30 normalized for equal weight per unit distance") +
	scale_x_continuous('Acceptor Base -- Acceptor Torsion (degrees)', breaks=c(90,270)) +
	scale_y_continuous("FeatureDensity")

if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)



f <- transform(f,
	capx = 2*sin(acos(cosBAH)/2)*cos(chi),
	capy = 2*sin(acos(cosBAH)/2)*sin(chi))

capx_limits <- range(f$capx)
capy_limits <- range(f$capy)
f <- ddply(f, c("sample_source", "alpha_type"),
	transform, counts = length(sample_source))

plot_id <- "hbond_backbone_backbone_alpha_type_chi_BAH"
d_ply(f, .(sample_source), function(sub_f){
	ss_id <- sub_f[1,"sample_source"]

	p <- ggplot(data=sub_f) + theme_bw() +
		theme(panel.background=element_rect(fill="#00007F", colour="#00007F")) +
		stat_bin2d(aes(x=capx, y=capy, fill=log(..density..)), binwidth=c(.01, .01)) +
		polar_equal_area_grids_bw() +
		geom_indicator(aes(indicator=counts, colour="white")) +
		coord_equal(ratio=1) +
		facet_wrap( ~ alpha_type) +
		ggtitle(paste("Backbone-Backbone H-Bonds CHI Angle by Alpha Type\nB-Factor < 30 equal coordinate projection Sample Source: ", ss_id, sep="")) +
		scale_x_continuous('Longitude: Acceptor Base -- Acceptor Torsion (degrees) ', limits=capx_limits, breaks=c()) +
		scale_y_continuous("Latitude: Acceptor Base -- Acceptor -- Hydrogen (degrees)", limits=capy_limits, breaks=c())

	if(nrow(sample_sources) <= 3){
		p <- p + theme(legend.position="bottom", legend.direction="horizontal")
	}

	save_plots(self, plot_id, sample_sources[sample_sources$sample_source == ss_id,], output_dir, output_formats)
})

})) # end FeaturesAnalysis
