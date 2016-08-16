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
id = "salt_bridge_geo_dim_2d",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("ResidueFeatures", "SaltBridgeFeatures", "PdbDataFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

source("scripts/analysis/plots/salt_bridges/salt_bridge_scales.R")

sele <-"
SELECT
	sb.psi, sb.theta, sb.rho, sb.orbital,
	don.name3 AS don_res_type, acc.HBChemType AS acc_chem_type,
	CASE don.resNum - acc.resNum
		WHEN -1 THEN '-1' WHEN -2 THEN '-2' WHEN -3 THEN '-3' WHEN -4 THEN '-4'
		WHEN 1 THEN '1' WHEN 2 THEN '2' WHEN 3 THEN '3' WHEN 4 THEN '4'
		ELSE 'long' END AS seq_sep
FROM
	salt_bridges AS sb,
	hbond_sites AS acc,
	residues as don,
	residue_pdb_confidence AS don_con,
	hbond_sites_pdb AS acc_con
WHERE
	don.struct_id = sb.struct_id AND don.resNum = sb.don_resNum AND
	acc.struct_id = sb.struct_id AND acc.site_id = sb.acc_id AND
	don_con.struct_id = sb.struct_id AND don_con.residue_number = don.resNum AND
	acc_con.struct_id = sb.struct_id AND acc_con.site_id = sb.acc_id AND
	don_con.max_sc_temperature < 30 AND acc_con.heavy_atom_temperature < 30 AND
	sb.rho > 2;"

f <- query_sample_sources(sample_sources, sele)

# give more descriptive plot labels
f$acc_chem_type <- factor(f$acc_chem_type,
	levels = c("hbacc_CXA", "hbacc_CXL"),
	labels = c("aCXA: n,q", "aCXL: d,e"))

f$acc_type <- with(f, interaction(acc_chem_type, orbital, sep="/"))


plot_parts <- list(
	theme_bw(),
	geom_point(size=.7),
	stat_density2d(size=.2),
	geom_indicator(aes(indicator=counts), group=1))


d_ply(f, .(sample_source), function(sub_f){
	ss_id <- sub_f$sample_source[1]
	ss <- sample_sources[sample_sources$sample_source == ss_id, ]

	plot_id <- paste("salt_bridge_psi_rho", ss_id, sep="_")
	sub_f$counts <- nrow(sub_f)
	p <- ggplot(data=sub_f, aes(x=psi, y=rho)) + plot_parts +
		ggtitle(paste("Salt Bridge PSI vs RHO, B-Factor < 30\nss_id: ", ss_id,sep="")) +
		scale_x_psi + scale_y_rho
	save_plots(self, plot_id, ss, output_dir, output_formats)


	plot_id <- paste("salt_bridge_psi_rho_by_acc_type_don_res_type", ss_id, sep="_")
	sub_f <- ddply(sub_f, .(acc_type, don_res_type), transform, counts=length(sample_source))
	p <- ggplot(data=sub_f, aes(x=psi, y=rho)) + plot_parts +
		ggtitle(paste("Salt Bridge PSI vs RHO, B-Factor < 30\nss_id: ", ss_id, sep="")) +
		facet_grid( acc_type ~ don_res_type ) +
		scale_x_psi + scale_y_rho
	save_plots(self, plot_id, ss, output_dir, output_formats)


	plot_id <- paste("salt_bridge_psi_rho_by_orbital_don_res_type", ss_id, sep="_")
	sub_f <- ddply(sub_f, .(orbital, don_res_type), transform, counts=length(sample_source))
	p <- ggplot(data=sub_f, aes(x=psi, y=rho)) + plot_parts +
		ggtitle(paste("Salt Bridge PSI vs RHO, B-Factor < 30\nss_id: ", ss_id, sep="")) +
		facet_grid( orbital ~ don_res_type ) +
		scale_x_psi + scale_y_rho
	save_plots(self, plot_id, ss, output_dir, output_formats)


	plot_id <- paste("salt_bridge_psi_rho_by_acc_chem_type_don_res_type", ss_id, sep="_")
	sub_f <- ddply(sub_f, .(acc_chem_type, don_res_type), transform, counts=length(sample_source))
	p <- ggplot(data=sub_f, aes(x=psi, y=rho)) + plot_parts +
		ggtitle(paste("Salt Bridge PSI vs RHO, B-Factor < 30\nss_id: ", ss_id, sep="")) +
		facet_grid( acc_chem_type ~ don_res_type ) +
		scale_x_psi + scale_y_rho
	save_plots(self, plot_id, ss, output_dir, output_formats)

	plot_id <- paste("salt_bridge_psi_rho_by_orbital_acc_chem_type", ss_id, sep="_")
	sub_f <- ddply(sub_f, .(orbital, acc_chem_type), transform, counts=length(sample_source))
	p <- ggplot(data=sub_f, aes(x=psi, y=rho)) + plot_parts +
		ggtitle(paste("Salt Bridge PSI vs RHO, B-Factor < 30\nss_id: ", ss_id, sep="")) +
		facet_grid( orbital ~ acc_chem_type ) +
		scale_x_psi + scale_y_rho
	save_plots(self, plot_id, ss, output_dir, output_formats)

	plot_id <- paste("salt_bridge_psi_rho_by_seq_sep_don_res_type", ss_id, sep="_")
	sub_f <- ddply(sub_f, .(seq_sep, don_res_type), transform, counts=length(sample_source))
	p <- ggplot(data=sub_f, aes(x=psi, y=rho)) + plot_parts +
		ggtitle(paste("Salt Bridge PSI vs RHO, B-Factor < 30\nss_id: ", ss_id, sep="")) +
		facet_grid( seq_sep ~ don_res_type ) +
		scale_x_psi + scale_y_rho
	save_plots(self, plot_id, ss, output_dir, output_formats)








	plot_id <- paste("salt_bridge_theta_rho", ss_id, sep="_")
	sub_f$counts <- nrow(sub_f)
	p <- ggplot(data=sub_f, aes(x=theta, y=rho)) + plot_parts +
		ggtitle(paste("Salt Bridge THETA vs RHO, B-Factor < 30\nss_id: ", ss_id, sep="")) +
		scale_x_theta + scale_y_rho

	save_plots(self, plot_id, ss, output_dir, output_formats)


	plot_id <- paste("salt_bridge_theta_rho_by_acc_type_don_res_type", ss_id, sep="_")
	sub_f <- ddply(sub_f, .(orbital, don_res_type), transform, counts=length(sample_source))
	p <- ggplot(data=sub_f, aes(x=theta, y=rho)) + plot_parts +
		ggtitle(paste("Salt Bridge THETA vs RHO, B-Factor < 30\nss_id: ", ss_id, sep="")) +
		facet_grid( acc_type ~ don_res_type ) +
		scale_x_theta + scale_y_rho
	save_plots(self, plot_id, ss, output_dir, output_formats)

	plot_id <- paste("salt_bridge_theta_rho_by_orbital_don_res_type", ss_id, sep="_")
	sub_f <- ddply(sub_f, .(orbital, don_res_type), transform, counts=length(sample_source))
	p <- ggplot(data=sub_f, aes(x=theta, y=rho)) + plot_parts +
		ggtitle(paste("Salt Bridge THETA vs RHO, B-Factor < 30\nss_id: ", ss_id, sep="")) +
		facet_grid( orbital ~ don_res_type ) +
		scale_x_theta + scale_y_rho
	save_plots(self, plot_id, ss, output_dir, output_formats)


	plot_id <- paste("salt_bridge_theta_rho_by_acc_chem_type_don_res_type", ss_id, sep="_")
	sub_f <- ddply(sub_f, .(acc_chem_type, don_res_type), transform, length(sample_source))
	p <- ggplot(data=sub_f, aes(x=theta, y=rho)) + plot_parts +
		ggtitle(paste("Salt Bridge THETA vs RHO, B-Factor < 30\nss_id: ", ss_id, sep="")) +
		facet_grid( acc_chem_type ~ don_res_type ) +
		scale_x_theta + scale_y_rho
	save_plots(self, plot_id, ss, output_dir, output_formats)

	plot_id <- paste("salt_bridge_theta_rho_by_orbital_acc_chem_type", ss_id, sep="_")
	sub_f <- ddply(sub_f, .(orbital, acc_chem_type), transform, length(sample_source))
	p <- ggplot(data=sub_f, aes(x=theta, y=rho)) + plot_parts +
		ggtitle(paste("Salt Bridge THETA vs RHO, B-Factor < 30\nss_id: ", ss_id, sep=""))+
		facet_grid( orbital ~ acc_chem_type ) +
		scale_x_theta + scale_y_rho
	save_plots(self, plot_id, ss, output_dir, output_formats)

	plot_id <- paste("salt_bridge_theta_rho_by_seq_sep_don_res_type", ss_id, sep="_")
	sub_f <- ddply(sub_f, .(seq_sep, don_res_type), transform, length(sample_source))
	p <- ggplot(data=sub_f, aes(x=theta, y=rho)) + plot_parts +
		ggtitle(paste("Salt Bridge THETA vs RHO, B-Factor < 30\nss_id: ", ss_id, sep=""))+
		facet_grid( seq_sep ~ don_res_type ) +
		scale_x_theta + scale_y_rho
	save_plots(self, plot_id, ss, output_dir, output_formats)





	plot_parts <- list(
		theme_bw(),
		geom_line(aes(x=x, y=log(y+1), colour=sample_source)),
		geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source), group=1),
		scale_y_continuous("log(FeatureDensity + 1)"),
		scale_x_rho)

	plot_id <- paste("salt_bridge_rho", ss_id, sep="_")
	dens <- estimate_density_1d(f, c("sample_source"), "rho", radial_3d_normalization)
	dens$y <- log(dens$y + 1)
	p <- ggplot(data=dens) + plot_parts +
		ggtitle("Salt Bridge RHO, B-Factor < 30\nnormalized for equal weight per unit distance")
	if(nrow(sample_sources) <= 3){
		p <- p + theme(legend.position="bottom", legend.direction="horizontal")
	}
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)


	plot_id <- paste("salt_bridge_rho_by_orbital_don_res_type", ss_id, sep="_")
	dens <- estimate_density_1d(f, c("sample_source", "don_res_type", "orbital"), "rho", radial_3d_normalization)
	p <- ggplot(data=dens) + plot_parts +
		ggtitle("Salt Bridge RHO, B-Factor < 30\nnormalized for equal weight per unit distance") +
		facet_grid( orbital ~ don_res_type )
	if(nrow(sample_sources) <= 3){
		p <- p + theme(legend.position="bottom", legend.direction="horizontal")
	}
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)


	plot_id <- paste("salt_bridge_rho_by_acc_chem_type_don_res_type", ss_id, sep="_")
	dens <- estimate_density_1d(f, c("sample_source", "don_res_type", "acc_chem_type"), "rho", radial_3d_normalization)
	p <- ggplot(data=dens) + plot_parts +
		ggtitle("Salt Bridge RHO, B-Factor < 30\nnormalized for equal weight per unit distance") +
		facet_grid( acc_chem_type ~ don_res_type )
	if(nrow(sample_sources) <= 3){
		p <- p + theme(legend.position="bottom", legend.direction="horizontal")
	}
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)

	plot_id <- paste("salt_bridge_rho_by_orbital_acc_chem_type", ss_id, sep="_")
	dens <- estimate_density_1d(f, c("sample_source", "orbital", "acc_chem_type"), "rho", radial_3d_normalization)
	p <- ggplot(data=dens) + plot_parts +
		ggtitle("Salt Bridge RHO, B-Factor < 30\nnormalized for equal weight per unit distance") +
		facet_grid( orbital ~ acc_chem_type )
	if(nrow(sample_sources) <= 3){
		p <- p + theme(legend.position="bottom", legend.direction="horizontal")
	}
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)

	plot_id <- paste("salt_bridge_rho_by_seq_sep_don_res_type", ss_id, sep="_")
	dens <- estimate_density_1d(f, c("sample_source", "seq_sep", "don_res_type"), "rho", radial_3d_normalization)
	p <- ggplot(data=dens) + plot_parts +
		ggtitle("Salt Bridge RHO, B-Factor < 30\nnormalized for equal weight per unit distance") +
		facet_grid( seq_sep ~ don_res_type )
	if(nrow(sample_sources) <= 3){
		p <- p + theme(legend.position="bottom", legend.direction="horizontal")
	}
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
})


})) # end FeaturesAnalysis
