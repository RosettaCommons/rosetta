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
id = "component_scores_one_body_by_radius_of_gyration",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("ResidueScoresFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


sele_1b <-"
SELECT
	rg.radius_of_gyration,
	rs.score_type,
	rs.score_value
FROM
	radius_of_gyration AS rg,
	residue_scores_1b as rs
WHERE
	rs.struct_id = rg.struct_id;"

scores <- query_sample_sources(sample_sources, sele_1b)
scores$score_type <- factor(scores$score_type)

d_ply(scores, .(score_type), function(sub_scores){
	score_type <- sub_scores$score_type[1]
	plot_id <- paste("component_scores_", score_type, "_by_radius_of_gyration", sep="")
	p <- ggplot(data=sub_scores) + theme_bw() +
#          	geom_point(aes(x=radius_of_gyration, y=score_value, colour=sample_source)) +
		stat_smooth(aes(x=radius_of_gyration, y=score_value, colour=sample_source), method=lm) +
#		geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
		ggtitle(paste("Rosetta ", score_type, " Scores", sep="")) +
		labs(x="Radius of Gyration (Angstroms)") +
		scale_y_continuous(paste(score_type, " Rosetta Energy", sep=""))
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)

})


})) # end FeaturesAnalysis
