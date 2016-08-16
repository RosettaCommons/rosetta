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
id = "pi_pi_hist_scOrbCpi_scOrbCpi_by_angle",
author = "Matthew O'Meara, Steven Combs",
brief_description = "",
feature_reporter_dependencies = c("OrbitalFeatures"),
run=function(self, sample_sources, output_dir, output_formats){
	
	sele <- "
			SELECT
			resNum1,
			resName1,
			resNum2,
			resName2,
			OrbName1,
			OrbOrbdist,
			AOO_angle
			FROM
			orbital_orbital
			WHERE 
			OrbOrbdist < 2.5 AND 
			(
			(OrbName1 = 'C.pi.sp2' AND OrbName2 = 'C.pi.sp2')  
			) AND 
			ABS(resNum1 - resNum2) > 5;"
	
	
	all_geom <- query_sample_sources(sample_sources, sele)
#ss_id <- all_geom$sample_sources[1]
#plot_id = paste("hist_AOH_scOrbPi_scHaro_all_residues", ss_id, sep="_")
	plot_id = "pi_pi_hist_AOO_scOrbCpi_scOrbCpi_all_residues"
	ggplot(data=all_geom) +
			geom_freqpoly(aes(x=AOO_angle, fill=sample_source, color = sample_source), binwidth=0.1) +
			ggtitle("AOO Pi-Pi measured via Orb - Orb \n",
					"scOrbCpi to scOrbCpi Combined: TYR, PHE, TRP") +
			scale_x_continuous(breaks=c(-1,-.9,-.8,-.7,-.6,-.5,-.4, -.3, -.2, -.1, 0, .1, .2, .3, .4, .5, .6, .7, .8), paste('Acceptor -- Orbital -- Orbital (degrees)'))
	scale_y_continuous("scOrbCpi ScOrbCpi Counts < 2.5 A from Orbital -- Orbital")
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
	
	
	
	
})) # end FeaturesAnalysis
