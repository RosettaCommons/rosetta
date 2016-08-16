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
id = "rama_zoom_helix_LEU_chi1-60",
author = "Matthew O'Meara",
brief_description = "Ramachandran plots conditional on the first sidechain torsional angle for LEU when chi1 is in the -60 degree bin",
feature_reporter_dependencies = c("ResidueFeatures", "ProteinBackboneTorsionAngleFeatures", "ResidueSecondaryStructureFeatures", "PdbDataFeatures", "ProteinResidueConformationFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


sele <-"
SELECT
	res.name3 AS res_type,
	CASE
-- figure out which is gauche+, gauche- and trans
		WHEN 0 < res_dofs.chi1 AND res_dofs.chi1 < 120 THEN 'chi1 60'
		WHEN -120 < res_dofs.chi1 AND res_dofs.chi1 < 0 THEN 'chi1 -60'
	  ELSE 'chi1 180' END AS chi1_bin,
	bb.phi, bb.psi
FROM
	residues AS res,
	residue_pdb_confidence AS res_conf,
	protein_residue_conformation AS res_dofs,
	protein_backbone_torsion_angles AS bb
WHERE
	res_conf.struct_id = res.struct_id AND res_conf.residue_number = res.resNum AND
	res_conf.max_temperature < 30 AND
	res_dofs.struct_id = res.struct_id AND res_dofs.seqpos == res.resNum AND
	bb.struct_id = res.struct_id AND bb.resNum == res.resNum AND
	res.name3 == 'LEU' AND	-120 < res_dofs.chi1 AND res_dofs.chi1 < 0 AND
	-75 < bb.phi AND bb.phi < -55 AND
	-50 < bb.psi AND bb.psi < -25;"

f <- query_sample_sources(sample_sources, sele)
narrow_output_formats <- transform(output_formats, width=height)


d_ply(f, .(sample_source, res_type, chi1_bin), function(sub_f){
	if(nrow(sub_f) < 5) {
		print("skipping this group:")
		print(summary(sub_f))
		return()
	}

	res_type <- as.character(sub_f[1, c("res_type")])
	chi1_bin_numeric <- as.numeric(sub_f[1, c("chi1_bin")])
	chi1_bin_character <- as.character(sub_f[1, c("chi1_bin")])
	ss_id <- as.character(sub_f[1, c("sample_source")])
	sub_plot_id <- paste(
		"ramachandran_zoom_helix_region_res_truncated_peaks-", res_type, "_chi1-", chi1_bin_numeric, "_", ss_id,
		sep="")

	dens <- MASS::kde2d(sub_f$phi, sub_f$psi, n=300, h=1)
	densdf <- data.frame(expand.grid(x=dens$x, y=dens$y), z=as.vector(dens$z))
	densdf$counts <- nrow(sub_f)
	ggplot(data=densdf) +
		theme_bw() +
		theme(panel.background=element_rect(fill="#00007F", colour="#00007F")) +
		geom_tile(aes(x=x, y=y, fill=z)) +
		geom_indicator(aes(indicator=counts), color="white") +
		coord_equal(ratio=1) +
		scale_fill_gradientn('Density', colours=jet.colors(15), limits=c(0, .01)) +
		theme(legend.position="bottom", legend.direction="horizontal") +
		ggtitle(paste(
			"Backbone Torsion Angles Res (truncated peaks): ", res_type, " ",
			chi1_bin_character, " B-Factor < 30\nSample Source: ",
			ss_id, sep="")) +
		scale_x_continuous(expression(paste("phi Angle (Degrees)", sep="")), limits=c(-75, -55)) +
		scale_y_continuous(expression(paste("psi Angle (Degrees)", sep="")), limits=c(-50, -25))

	save_plots(self, sub_plot_id, sample_sources, output_dir, narrow_output_formats)


	sub_plot_id <- paste(
		"ramachandran_zoom_helix_region_res-", res_type, "_chi1-", chi1_bin_numeric, "_", ss_id,
		sep="")

	ggplot(data=densdf) +
		theme_bw() +
		theme(panel.background=element_rect(fill="#00007F", colour="#00007F")) +
		geom_tile(aes(x=x, y=y, fill=z)) +
		geom_indicator(aes(indicator=counts), color="white") +
		coord_equal(ratio=1) +
		scale_fill_gradientn('Density', colours=jet.colors(15)) +
		theme(legend.position="bottom", legend.direction="horizontal") +
		ggtitle(paste(
			"Backbone Torsion Angles Res: ", res_type, " ",
			chi1_bin_character, " B-Factor < 30\nSample Source: ",
			ss_id, sep="")) +
		scale_x_continuous(expression(paste("phi Angle (Degrees)", sep="")), limits=c(-75, -55)) +
		scale_y_continuous(expression(paste("psi Angle (Degrees)", sep="")), limits=c(-50, -25))

	save_plots(self, sub_plot_id, sample_sources, output_dir, narrow_output_formats)

})

})) # end FeaturesAnalysis
