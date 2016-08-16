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
id = "rama_not_helix_tight_kenrel",
filename = "scripts/analysis/plots/backbone_geometry/rama_not_helix_tight_kernel.R",
author = "Matthew O'Meara",
brief_description = "Overview of ramachandran distribution for all non Alpha-helix residues. The Alpha-helix residues are not included because their density concentration obscures the other features of the ramachandran distribution.",
feature_reporter_dependencies = c("ResidueFeatures", "ProteinBackboneTorsionAngleFeatures", "PdbDataFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


sele <-"
SELECT
	bb.phi, bb.psi,
	CASE ss.dssp WHEN 'H' THEN 1 ELSE 0 END AS is_helix
FROM
	residues AS res,
	residue_pdb_confidence AS res_conf,
	residue_secondary_structure AS ss,
	protein_backbone_torsion_angles AS bb
WHERE
	res_conf.struct_id = res.struct_id AND res_conf.residue_number = res.resNum AND
	res_conf.max_temperature < 30 AND
	ss.struct_id = res.struct_id AND ss.resNum == res.resNum AND
	ss.dssp != 'H' AND
	bb.struct_id = res.struct_id AND bb.resNum == res.resNum;"

f <- query_sample_sources(sample_sources, sele)

f <- ddply(
	f, .(sample_source),
	transform, counts = length(sample_source))

plot_parts <- list(
	theme_bw(),
	geom_raster(aes(x=x, y=y, fill=log(log(z+1)+1))),
	geom_indicator(aes(indicator=counts), color="white"),
	coord_equal(ratio=1),
	scale_x_continuous(expression(paste("phi Angle (Degrees)", sep=""))),
	scale_y_continuous(expression(paste("psi Angle (Degrees)", sep=""))),
	scale_fill_gradientn('Log(Density)', colours=jet.colors(15)),
	theme(
		legend.position="bottom",
		legend.direction="horizontal",
		panel.grid.major=theme_blank(),
		panel.grid.minor=theme_blank(),
		panel.background=theme_rect(fill="#00007F")))

narrow_output_formats <- transform(output_formats, width=height)

d_ply(f, .(sample_source), function(sub_f){
	if(nrow(sub_f) < 5) return()

	ss_id <- as.character(sub_f[1, c("sample_source")])
	sub_plot_id <- paste("ramachandran_", ss_id, sep="")

	t <- system.time({
		dens <- MASS::kde2d(sub_f$phi, sub_f$psi, n=500, h=1)
		densdf <- data.frame(expand.grid(x=dens$x, y=dens$y), z=as.vector(dens$z)/max(dens$z))
		densdf$counts <- nrow(sub_f)
	})
	cat("estimate density ... ", as.character(round(t[3], 2)), "s\n", sep="")

	ggplot(data=densdf) + plot_parts +
		ggtitle(paste("Backbone Torsion Angles, non Alpha-Helix, B-Factor < 30\nSample Source: ", ss_id, sep=""))

	save_plots(self, sub_plot_id, sample_sources, output_dir, narrow_output_formats)
})

})) # end FeaturesAnalysis
