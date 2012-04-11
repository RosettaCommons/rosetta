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
id = "cation_pi_heatmap_scOrbCpi_Hpol_by_angle_dist",
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
			OrbHdist,
			AOH_angle
			FROM
			HPOL_orbital
			WHERE
			(((resName1 = 'TYR' OR resName1 = 'PHE' OR resName1 = 'TRP') AND (resName2 = 'ARG' OR resName2 = 'LYS')) OR 
			((resName1 = 'LYS' OR resName1 = 'ARG' ) AND (resName2 = 'TYR' OR resName2 = 'PHE' OR resName2 = 'TRP'))) 
			AND OrbHdist < 3;"
	
	
	f <- query_sample_sources(sample_sources, sele)
	
	plot_id = "stat_bin2d_AOH_scOrbPi_Hpol_all_residues"
	ggplot(f, aes(x=AOH_angle, y=OrbHdist)) +
			stat_bin2d(aes(x=AOH_angle, y=OrbHdist, fill=..count..)) + 
			scale_fill_gradient(low = "white", high = "steelblue")
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
	
	
	plot_id = "stat_density2d_AOH_scOrbPi_Hpol_all_residues"
	ggplot(f, aes(x=AOH_angle, y=OrbHdist)) + stat_density2d(aes(x=AOH_angle, y=OrbHdist, fill = ..density..), geom="raster", contour=F) + 
			scale_fill_gradient(low = "white", high = "steelblue")	
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
	
	
#	(p <- ggplot(nba.m, aes(variable, Name)) + 
#				geom_tile(aes(fill = rescale), + colour = "white") + scale_fill_gradient(low = "white",
#						+     high = "steelblue"))
#	plot_id = "hist_DHO_scOrbPi_Hpol_all_residues"
#	ggplot(data=all_geom) +
#			geom_bar(aes(x=DHO_angle, fill=sample_source), binwidth=0.1, position="dodge") +
#			opts(title = "DHO Cation - Pi Measured Via Orb - Hpol \n",
#						"ScOrbCpi to Hpol TYR/PHE/TRP to LYS/ARG") +
#			scale_x_continuous(breaks=c(-1,-.9,-.8,-.7,-.6,-.5,-.4, -.3, -.2, -.1, 0, .1, .2, .3, .4, .5, .6, .7, .8), paste('Donor -- Hydrogen -- Orbital'))
#	scale_y_continuous("scOrbCpi scHpol Counts < 3.0 A from Orbital -- Hydrogen")
#	save_plots(self, plot_id, sample_sources, output_dir, output_formats)

})) # end FeaturesAnalysis
