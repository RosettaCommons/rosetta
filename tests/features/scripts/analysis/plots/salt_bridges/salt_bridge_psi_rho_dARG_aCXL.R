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
id = "salt_bridge_psi_rho_dARG_aCXL",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("ResidueFeatures", "SaltBridgeFeatures", "PdbDataFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
	sb.psi, sb.rho,
	CASE don.resNum - acc.resNum
		WHEN -1 THEN '-1' WHEN -2 THEN '-2' WHEN -3 THEN '-3' WHEN -4 THEN '-4'
		WHEN 1 THEN '1' WHEN 2 THEN '2' WHEN 3 THEN '3' WHEN 4 THEN '4'
		ELSE 'long' END AS seq_sep
FROM
	salt_bridges AS sb,
	hbond_sites AS acc,
	residues AS don
--	residue_pdb_confidence AS don_con,
--	hbond_sites_pdb AS acc_con
WHERE
	don.struct_id = sb.struct_id AND don.resNum = sb.don_resNum AND
	acc.struct_id = sb.struct_id AND acc.site_id = sb.acc_id AND
	don.name3 == 'ARG' AND acc.HBChemType == 'hbacc_CXL';"
#	don_con.struct_id = sb.struct_id AND don_con.residue_number = don.resNum AND
#	acc_con.struct_id = sb.struct_id AND acc_con.site_id = sb.acc_id AND
#	don_con.max_sc_temperature < 30 AND acc_con.heavy_atom_temperature < 30;"

f <- query_sample_sources(sample_sources, sele)
f$psi_degrees <- f$psi*180/pi
f$cluster <- add_salt_bridge_clusters(f)
f <- ddply(f, .(sample_source), transform, counts = length(sample_source))


plot_id <- "salt_bridge_psi_rho"

p <- ggplot() +
#	geom_polygon(
#		data=clusters,
#		aes(x=x, y=y, fill=cluster_label, group=cluster_label),
#		show_guide=F) +
	geom_point(
		data=f,
		aes(x=psi_degrees, y=rho, colour=cluster),
		size=.7) +
	stat_density2d(
		data=f,
		aes(x=psi_degrees, y=rho, colour=cluster),
		size=.2) +
	geom_indicator(
		data=f,
		aes(indicator=counts),
		group=1) +
	ggtitle("Salt Bridge PSI vs RHO, B-Factor < 30") +
	facet_wrap(~sample_source, ncol=2) +
	theme_bw() +
	scale_y_continuous(expression(
		paste('Cental Carbon -- Oxygen Distance (', ring(A), ')')),
		limit=c(2.8,4.5), breaks=seq(2.8, 4.5, .1)) +
	scale_x_continuous("Angle Around Donor (Degrees)", breaks=seq(-180, 180, 20))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)



})) # end FeaturesAnalysis
