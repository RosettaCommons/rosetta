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
id = "dssp_content",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("ResidueSecondaryStructureFeatures"),
run=function(self){

sele <-"
SELECT
  ss.dssp,
  COUNT( ss.dssp ) AS counts,
  num_res.total AS num_res
FROM
  residue_secondary_structure AS ss,
  (SELECT COUNT(*) AS total FROM residue_secondary_structure) AS num_res
GROUP BY
  ss.dssp;"

f <- query_sample_sources(sample_sources, sele)

# filter out(?) bad data
f <- f[f$counts > 10,]
f$fraction <- f$counts / f$num_res

f$dssp_description <- NA
f[f$dssp == "H", "dssp_description"] <- "H: alpha helix"
f[f$dssp == "B", "dssp_description"] <- "B: isolated beta-bridge"
f[f$dssp == "E", "dssp_description"] <- "E: extended strand, participates in beta ladder"
f[f$dssp == "G", "dssp_description"] <- "G: 3-helix (3/10 helix)"
f[f$dssp == "I", "dssp_description"] <- "I: 5 helix (pi helix)"
f[f$dssp == "T", "dssp_description"] <- "T: hydrogen bonded turn"
f[f$dssp == "S", "dssp_description"] <- "S: bend"
f[f$dssp == " ", "dssp_description"] <- "loop or irregular"
f$dssp_description <- factor(f$dssp_description)

plot_id <- "dssp_content"
ggplot(data=f) + theme_bw() +
	geom_bar(aes(x=sample_source, y=fraction, fill=sample_source)) +
	facet_wrap(~dssp_description, ncol=1) +
	opts(title = "Secondary Structure Content") +
	coord_flip() +
	opts(legend.position="none") +
        labs(y="Fraction of All Residues in Each Sample Source")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

})) # end FeaturesAnalysis
