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
id = "cation_pi_hist_scOrbPi_scOrbNpi_by_angle",
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
			(OrbName1 = 'C.pi.sp2' AND OrbName2 = 'N.pi.sp2') OR 
			(OrbName2 = 'C.pi.sp2' AND OrbName1 = 'N.pi.sp2') 
			) AND
	(((resName1 = 'TYR' OR resName1 = 'PHE' OR resName1 = 'TRP') AND (resName2 = 'ARG')) OR 
	((resName1 = 'ARG' ) AND (resName2 = 'TYR' OR resName2 = 'PHE' OR resName2 = 'TRP')));" 
	
	
	all_geom <- query_sample_sources(sample_sources, sele)
	plot_id = "cation_pi_hist_AOO_scOrbPi_scOrbNpi_all_residues"
	ggplot(data=all_geom) +
			geom_freqpoly(aes(x=AOO_angle, fill=sample_source, color = sample_source ), binwidth=0.1) +
			ggtitle("AOO Cation - Pi Measured Via Orb - Orb \n",
						 "scOrbCpi to scOrbNpi TYR/PHE/TRP to LYS/ARG") +
			scale_x_continuous(breaks=c(-1,-.9,-.8,-.7,-.6,-.5,-.4, -.3, -.2, -.1, 0, .1, .2, .3, .4, .5, .6, .7, .8), paste('Acceptor -- Orbital -- Orbital'))
	scale_y_continuous("scOrbNpi scOrbCpi Counts < 2.5 A from Orbital -- Orbital")
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
		
})) # end FeaturesAnalysis
