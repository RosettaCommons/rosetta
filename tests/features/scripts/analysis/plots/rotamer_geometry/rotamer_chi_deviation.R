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
id = "rotamer_chi_deviation",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("ResidueFeatures", "ProteinResidueConformationFeatures", "PdbDataFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


sele <- "
SELECT
	res.name3 AS res_type,
	rot.rotamer_bin, rot.rotamer_bin_probability,
	'chi1' AS variable,
	rot.chi1_deviation / rot.chi1_standard_deviation AS value
FROM
	residues AS res, residue_rotamers AS rot
WHERE
	rot.struct_id = res.struct_id AND rot.residue_number = res.resNum AND
	rot.rotamer_bin_probability > .001 AND
	rot.chi1_deviation IS NOT NULL
UNION
SELECT
	res.name3 AS res_type,
	rot.rotamer_bin, rot.rotamer_bin_probability,
	'chi2' AS variable,
	rot.chi2_deviation / rot.chi2_standard_deviation AS value
FROM
	residues AS res, residue_rotamers AS rot
WHERE
	rot.struct_id = res.struct_id AND rot.residue_number = res.resNum AND
	rot.rotamer_bin_probability > .001 AND
	rot.chi2_deviation IS NOT NULL
UNION
SELECT
	res.name3 AS res_type,
	rot.rotamer_bin, rot.rotamer_bin_probability,
	'chi3' AS variable,
	rot.chi3_deviation / rot.chi3_standard_deviation AS value
FROM
	residues AS res, residue_rotamers AS rot
WHERE
	rot.struct_id = res.struct_id AND rot.residue_number = res.resNum AND
	rot.rotamer_bin_probability > .001 AND
	rot.chi3_deviation IS NOT NULL
UNION
SELECT
	res.name3 AS res_type,
	rot.rotamer_bin, rot.rotamer_bin_probability,
	'chi4' AS variable,
	rot.chi4_deviation / rot.chi4_standard_deviation AS value
FROM
	residues AS res, residue_rotamers AS rot
WHERE
	rot.struct_id = res.struct_id AND rot.residue_number = res.resNum AND
	rot.rotamer_bin_probability > .001 AND
	rot.chi4_deviation IS NOT NULL;"
f <- query_sample_sources(sample_sources, sele)


f <- ddply(f, .(res_type), function(df){
	if(nrow(df) < 50){
		return(data.frame())
	} else {
		return(df)
	}
})

d_ply(f, .(variable), function(df){
	chi <- as.character(df$variable[1])

	plot_id <- paste("rotamer_normalized_", chi, "_deviation", sep="")
	p <- ggplot(df) +
		theme_bw() +
		geom_boxplot(
			aes(x=sample_source, y=value, color=sample_source),
			outlier.size=.7) +
		facet_wrap( ~ res_type, ncol=2) +
		theme(
#			title=paste("Rotamer ", chi, " Deviation", sep=""),
			legend.position="none") +
		scale_y_continuous("Normalized Angle Deviation", limits=c(-10, 10)) +
		scale_x_discrete("Energy Function") +
		coord_flip()
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
})
})) # end FeaturesAnalysis
