# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

check_setup()
feature_analyses <- c(feature_analyses, new("FeaturesAnalysis",
id = "omega_angle",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("ProteinBackboneTorsionAngleFeatures", "ResidueSecondaryStructureFeatures"),
run=function(self){


sele <-"
SELECT
	res.name3 AS res_type,
	ss.dssp,
	bb.omega,
	CASE WHEN abs(bb.omega) < 90 THEN 'cis' ELSE 'trans' END AS isomerization
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

f$omega <- ifelse(f$omega < 0, f$omega + 360, f$omega)

dens <- estimate_density_1d(f[f$isomerization == 'trans',],
	c("sample_source"), "omega")

plot_id <- "backbone_geometry_omega_angle_trans"
p <- ggplot(dens) + theme_bw() +
	geom_line(aes(x, y, colour=sample_source), size=2) +
	opts(title = "Trans Omega Angle For B-Factor < 30") +
	geom_vline(x=180) +
        scale_x_continuous("Omega Angle (degrees)", limit=c(150,210), breaks=c(150, 160, 170, 180, 190, 200, 210)) +
	scale_y_continuous("Feature Density")
if(nrow(sample_sources) <= 3){
	p <- p + opts(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
