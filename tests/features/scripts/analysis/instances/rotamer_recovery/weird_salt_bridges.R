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
id = "weird_salt_bridges",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("StructureFeatures", "ResidueFeatures", "RotamerRecoveryFeatures", "PdbDataFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
	sb.psi, sb.rho,
	struct.tag, '' AS chain, don.resNum,
	'CA' AS CA, 'C' AS C, 'CB' AS CB
FROM
	salt_bridges AS sb,
	residues as don, hbond_sites AS acc,
	structures AS struct
WHERE
	struct.struct_id = sb.struct_id AND
	don.struct_id = sb.struct_id AND don.resNum = sb.don_resNum AND
	acc.struct_id = sb.struct_id AND acc.site_id = sb.acc_id AND
	don.name3 = 'LYS' AND acc.HBChemType = 'hbacc_CXL' AND sb.orbital = 'syn' AND
	sb.rho > 2;"

f <-  query_sample_sources(sample_sources, sele)

psi_range <- c(-10, 10)
rho_range <- c(3.3, 3.6)

f$psi <- f$psi * 180/pi

f$instance <-
	f$psi > psi_range[1] & f$psi < psi_range[2] &
	f$rho > rho_range[1] & f$rho < rho_range[2]


plot_parts <- list(
	theme_bw(),
	geom_point(aes(colour=instance), size=.7),
	stat_density2d(size=.2),
	geom_indicator(aes(indicator=counts)),
	scale_y_continuous(expression(
		paste('Cental Carbon -- Oxygen Distance (', ring(A), ')')), limit=c(2,6)))

d_ply(f, .(sample_source), function(sub_f){
	ss_id <- sub_f$sample_source[1]
	ss <- sample_sources[sample_sources$sample_source == ss_id, ]

	sub_f$counts <- sum(sub_f$instance)
	plot_id <- "salt_bridge_instances_psi_rho_LYS_CXL"
	p <- ggplot(data=sub_f, aes(x=psi, y=rho)) + plot_parts +
		ggtitle(paste("Salt Bridge LYS donor D/E acceptor, PSI vs RHO by Acceptor\nss_id: ", ss_id,sep="")) +
		scale_x_continuous("Angle Around Donor (Degrees)")
	save_plots(self, plot_id, ss, output_dir, output_formats)
})

n_examples <- 15

g <- f[f$instance,]
g$id <- 1:nrow(g)
g <- melt(g[g$id <= n_examples,],
	id.vars=c("sample_source", "tag", "id", "chain", "resNum"),
	measure.vars=c("CA", "C", "CB"),
	variable_name = "atom")

instances_id <- "salt_bridge_instances_psi_rho_LYS_CXL"
prepare_feature_instances(instances_id, sample_sources, g, output_dir)

})) # end FeaturesAnalysis
