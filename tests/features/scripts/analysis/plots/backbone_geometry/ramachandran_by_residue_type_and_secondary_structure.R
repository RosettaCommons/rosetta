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
id = "ramachandran_by_residue_type_and_secondary_structure",
filename = "scripts/analysis/plots/backbone_geometry/ramachandran_by_residue_type_and_secondary_structure.R",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("ResidueFeatures", "ProteinBackboneTorsionAngleFeatures", "ResidueSecondaryStructureFeatures", "PdbDataFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


sele <-"
SELECT
	res.name3 AS res_type,
	ss.dssp,
	bb.phi, bb.psi
FROM
	residues AS res,
	residue_pdb_confidence AS res_conf,
	residue_secondary_structure AS ss,
	protein_backbone_torsion_angles AS bb
WHERE
	res_conf.struct_id = res.struct_id AND res_conf.residue_number = res.resNum AND
	res_conf.max_temperature < 30 AND
	ss.struct_id = res.struct_id AND ss.resNum == res.resNum AND
	bb.struct_id = res.struct_id AND bb.resNum == res.resNum;"

f <- query_sample_sources(sample_sources, sele)



f <- ddply(
	f, .(sample_source, res_type, dssp),
	transform, counts = length(sample_source))

plot_parts <- list(
	theme_bw(),
	theme(panel.background=element_rect(fill="#00007F", colour="#00007F")),
	stat_density2d(
		aes(x=phi, y=psi, fill=..density..), geom="tile", contour=FALSE),
	geom_indicator(aes(indicator=counts), color="white"),
	coord_equal(ratio=1),
	scale_x_continuous(expression(paste("phi Angle (Degrees)", sep="")), limits=c(-180, 180)),
	scale_y_continuous(expression(paste("psi Angle (Degrees)", sep="")), limits=c(-180, 180)),
	scale_fill_gradientn('Density', colours=jet.colors(15)),
	theme(legend.position="bottom", legend.direction="horizontal"))

narrow_output_formats <- transform(output_formats, width=height)

d_ply(f, .(sample_source, res_type, dssp), function(sub_f){
	if(nrow(sub_f) < 5) return()

	res_type <- as.character(sub_f[1, c("res_type")])
	dssp <- as.character(sub_f[1, c("dssp")])
	ss_id <- as.character(sub_f[1, c("sample_source")])
	sub_plot_id <- paste("ramachandran_res-", res_type, "_dssp-", dssp,"_", ss_id, sep="")
	ggplot(data=sub_f) + plot_parts +
		ggtitle(paste("Backbone Torsion Angles Res: ", res_type, " DSSP: ", dssp, " B-Factor < 30\nSample Source: ", ss_id, sep=""))
	save_plots(self, sub_plot_id, sample_sources, output_dir, narrow_output_formats)
})

})) # end FeaturesAnalysis
