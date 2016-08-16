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
id = "rotamer_recovery_differences_salt_bridge",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("ResidueFeatures", "RotamerRecoveryFeatures", "SaltBridgeFeatures", "PdbDataFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
	don_rr.divergence AS don_rr, acc_rr.divergence AS acc_rr,
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
	hbond_sites_pdb AS acc_con,
	rotamer_recovery AS don_rr, rotamer_recovery AS acc_rr
WHERE
	don.struct_id = sb.struct_id AND don.resNum = sb.don_resNum AND
	acc.struct_id = sb.struct_id AND acc.site_id = sb.acc_id AND
	don_con.struct_id = sb.struct_id AND don_con.residue_number = don.resNum AND
	acc_con.struct_id = sb.struct_id AND acc_con.site_id = sb.acc_id AND
	don_con.max_sc_temperature < 30 AND acc_con.heavy_atom_temperature < 30 AND
	sb.rho > 2 AND
	don_rr.struct_id = don.struct_id AND don_rr.resNum = don.resNum AND
	acc_rr.struct_id = acc.struct_id AND acc_rr.resNum = acc.resNum;"

f <- query_sample_sources(sample_sources, sele)

f[f$don_res_type == 'LYS',] <- transform(f[f$don_res_type == 'LYS',],
	theta = theta * 180/pi,
	psi = psi * 180/pi)

f[f$don_res_type == 'HIS' | f$don_res_type == 'ARG',] <- transform(f[f$don_res_type == 'HIS' | f$don_res_type == 'ARG',],
	theta = theta - 90)

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
	geom_point(size=.7),
	stat_density2d(size=.2),
	geom_indicator(aes(indicator=counts)),
	scale_y_continuous(expression(
		paste('Cental Carbon -- Oxygen Distance (', ring(A), ')')), limit=c(2,6)))


d_ply(f, .(sample_source), function(sub_f){
	ss_id <- sub_f$sample_source[1]
	ss <- sample_sources[sample_sources$sample_source == ss_id, ]

	plot_id <- "rotamer_recovery_donor_salt_bridge_psi_rho"
	sub_f$counts <- nrow(sub_f)
	p <- ggplot(data=sub_f, aes(x=psi, y=rho, color=don_rr, group=don_rr)) + plot_parts +
		ggtitle(paste("Salt Bridge PSI vs RHO, B-Factor < 30\nss_id: ", ss_id,sep="")) +
		scale_x_continuous("Angle Around Donor (Degrees)")
	save_plots(self, plot_id, ss, output_dir, output_formats)

	plot_id <- "rotamer_recovery_donor_salt_bridge_psi_rho_by_acc_type_don_res_type"
	sub_f <- ddply(sub_f, .(acc_type, don_res_type), transform, length(sample_source))
	p <- ggplot(data=sub_f, aes(x=psi, y=rho, color=don_rr, group=don_rr)) + plot_parts +
		ggtitle(paste("Salt Bridge PSI vs RHO, B-Factor < 30\nss_id: ", ss_id, sep="")) +
		facet_grid( acc_type ~ don_res_type ) +
		scale_x_continuous("Angle Around Donor (Degrees)")
	save_plots(self, plot_id, ss, output_dir, output_formats)

	plot_id <- "rotamer_recovery_donor_salt_bridge_psi_rho_by_seq_sep_don_res_type"
	sub_f <- ddply(sub_f, .(seq_sep, don_res_type), transform, length(sample_source))
	p <- ggplot(data=sub_f, aes(x=psi, y=rho, color=don_rr, group=don_rr)) + plot_parts +
		ggtitle(paste("Salt Bridge PSI vs RHO, B-Factor < 30\nss_id: ", ss_id, sep="")) +
		facet_grid( seq_sep ~ don_res_type ) +
		scale_x_continuous("Angle Around Donor")
	save_plots(self, plot_id, ss, output_dir, output_formats)








	plot_id <- "rotamer_recovery_donor_salt_bridge_theta_rho"
	sub_f$counts <- nrow(sub_f)
	p <- ggplot(data=sub_f, aes(x=theta, y=rho, color=don_rr, group=don_rr)) + plot_parts +
		ggtitle(paste("Salt Bridge THETA vs RHO, B-Factor < 30\nss_id: ", ss_id, sep="")) +
		scale_x_continuous("Angle Out of Donor Plane")
	save_plots(self, plot_id, ss, output_dir, output_formats)


	plot_id <- "rotamer_recovery_donor_salt_bridge_theta_rho_by_acc_type_don_res_type"
	sub_f <- ddply(sub_f, .(orbital, don_res_type), transform, length(sample_source))
	p <- ggplot(data=sub_f, aes(x=theta, y=rho, color=don_rr, group=don_rr)) + plot_parts +
		ggtitle(paste("Salt Bridge THETA vs RHO, B-Factor < 30\nss_id: ", ss_id, sep="")) +
		facet_grid( acc_type ~ don_res_type ) +
		scale_x_continuous("Angle Out of Donor Plane (Degrees)")
	save_plots(self, plot_id, ss, output_dir, output_formats)


	plot_id <- "rotamer_recovery_donor_salt_bridge_theta_rho_by_seq_sep_don_res_type"
	sub_f <- ddply(sub_f, .(seq_sep, don_res_type), transform, length(sample_source))
	p <- ggplot(data=sub_f, aes(x=theta, y=rho, color=don_rr, group=don_rr)) + plot_parts +
		ggtitle(paste("Salt Bridge THETA vs RHO, B-Factor < 30\nss_id: ", ss_id, sep=""))+
		facet_grid( seq_sep ~ don_res_type ) +
		scale_x_continuous("Angle Out of Donor Plane")
	save_plots(self, plot_id, ss, output_dir, output_formats)
})

})) # end FeaturesAnalysis
