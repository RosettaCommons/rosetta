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
id = "rotamer_recovery_vs_score",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("ResidueFeatures", "ResidueScoresFeatures", "RotamerRecoveryFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


res_types <- c(
					"LYS", "VAL", "ILE", "ASN", "GLU", "GLN", "ARG", "SER", "ASP",
					"LEU", "HIS", "THR", "PRO", "TYR", "TRP", "CYS", "PHE", "MET")

ldply(res_types, function(res_type){

  subplot_id <- paste(plot_id, res_type, sep="_")
	print( paste( "Begin processing res type '", res_type, "'", sep="" ))

	sele <-paste("
	SELECT
	  rr.divergence AS divergence,
	  sc.score_type AS score_type,
	  sc.score_value AS score_value
	FROM
	  rotamer_recovery AS rr,
	  residue_scores_1b AS sc,
	  residues AS res
	WHERE
	  rr.resNum = res.resNum AND
	  rr.struct_id = res.struct_id AND
	  sc.resNum = res.resNum AND
	  sc.struct_id = res.struct_id AND
		res.name3 = '",res_type,"' AND
		score_type != 'ref' AND
	  score_value > -5 AND
	  score_value < 5;",sep="")

	all_1b <-  query_sample_sources(sample_sources, sele)

  sele <-paste("
  SELECT
    rr.divergence AS divergence,
		CASE WHEN sc.score_type == 'fa_rep' THEN 'fa_atr/rep'
		     WHEN sc.score_type == 'fa_atr' THEN 'fa_atr/rep'
		     ELSE sc.score_type END AS score_type,
    sc.score_value AS score_value
  FROM
    rotamer_recovery AS rr,
    residue_scores_2b AS sc,
    residues AS res
  WHERE
    rr.resNum = res.resNum AND
    rr.struct_id = res.struct_id AND
    sc.resNum1 = res.resNum AND
    sc.struct_id = res.struct_id AND
		res.name3 = '",res_type,"' AND
    score_value > -5 AND
    score_value < 5;", sep="")

  all_2b <-  query_sample_sources(sample_sources, sele)

  all <- rbind( all_1b, all_2b )
  all$score_type <- factor(all$score_type)

	some <- ddply(all, c("score_type"), function(df) sample_rows( df, 4000 ) )

	plot_id <- "rotamer_recovery_vs_score"
  p <- ggplot(
    data=some,
    aes(x=score_value, y=log(divergence + 1)))
  p <- p + geom_point(size=.5)
  p <- p + facet_wrap( ~ score_type, scales="free_x" )
  p <- p + ggtitle(paste("Rotamer Recovery for", res_type, "by Score Type"))
  p <- p + labs(x="Score",
                y="log(AutomorphicRMSD + 1)")
  p <- p + theme_bw()

  save_plots(self, subplot_id, sample_sources, output_dir, output_formats)
})

})) # end FeaturesAnalysis
