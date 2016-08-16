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
id = "rama_zoom_all_beta",
author = "Matthew O'Meara",
brief_description = "Ramachandran plots conditional on the first sidechain torsional angle for THR when chi1 is in the -60 degree bin",
feature_reporter_dependencies = c("ResidueFeatures", "ProteinBackboneTorsionAngleFeatures", "ResidueSecondaryStructureFeatures", "PdbDataFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
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
	-150 < bb.phi AND bb.phi < -75 AND
	100 < bb.psi AND bb.psi < 150;"

phi_range <- c(-150, -75)
psi_range <- c(100, 150)

f <- query_sample_sources(sample_sources, sele)

dens <- ddply(f, .(sample_source), function(sub_f) {
	dens <- MASS::kde2d(sub_f$phi, sub_f$psi, n=500, h=1)
	densdf <- data.frame(expand.grid(x=dens$x, y=dens$y), z=as.vector(dens$z)/max(dens$z))
	densdf$counts <- nrow(sub_f)
	densdf
})

if (nrow(sample_sources) > 24) {
	stop("Unable to label facet panel because there arn't enough letters in the alphebet!")
}

facet_labels <- data.frame(
	sample_source = sample_sources$sample_source,
	label=toupper(letters[1:nrow(sample_sources)]))


plot_id <- "ramachandran_zoom_beta"
ggplot(data=dens) +
	theme_bw() +
	geom_raster(aes(x=x, y=y, fill=z)) +
	geom_indicator(aes(indicator=counts), group=1, color="black", ypos=.99) +
	geom_indicator(
		data=facet_labels,
		aes(indicator=label),
		group=1,
		color="black",
		xpos=.03,
		ypos=.03,
		size=20) +
	coord_equal(ratio=1) +
	scale_fill_gradient('Scaled Density', low="white", high="black") +
	theme(legend.position="bottom", legend.direction="horizontal") +
#	ggtitle(paste(
#		"Backbone Torsion Angles in the Beta-Region; B-Factor < 30", sep="")) +
	scale_x_continuous(expression(paste("phi Angle (Degrees)", sep="")), limits=phi_range) +
	scale_y_continuous(expression(paste("psi Angle (Degrees)", sep="")), limits=psi_range) +
	facet_wrap(~sample_source)
save_plots(self, plot_id, sample_sources, output_dir, output_formats)



})) # end FeaturesAnalysis
