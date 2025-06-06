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
id = "pi_pi_heatmap_scOrbCpi_Haro_by_angle_dist",
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
	AOH_angle,
	DHO_angle
FROM
	HARO_orbital
WHERE
	(((resName1 = 'TYR' OR resName1 = 'PHE' OR resName1 = 'TRP') AND (resName2 = 'TYR' OR resName2 = 'PHE' OR resName2 = 'TRP')) OR 
	((resName1 = 'TYR' OR resName1 = 'PHE' OR resName1= 'TRP') AND (resName2 = 'TYR' OR resName2 = 'PHE' OR resName2 = 'TRP'))) 
	AND OrbHdist < 4;"
	
	
	f <- query_sample_sources(sample_sources, sele)
	
	plot_id = "pi_pi_stat_bin2d_heatmap_AOH_scOrbPi_Haro_all_residues"
	ggplot(f, aes(x=AOH_angle, y=OrbHdist)) +
			stat_bin2d(aes(x=AOH_angle, y=OrbHdist, fill=..count..), binwidth = c(0.1, 0.5), xlim =c(-1, 0), ylim = c(0, 3)) + 
			scale_fill_gradient(low = "white", high = "steelblue") +
			scale_x_discrete(paste('Acceptor -- Orbital -- Hydrogen')) +  
			scale_y_discrete(paste('scOrbCpi scHaro Counts < 3.0 A from Orbital -- Hydrogen'))
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
	
	
	plot_id = "pi_pi_stat_density2d_heatmap_AOH_scOrbPi_Haro_all_residues"
	ggplot(f, aes(x=AOH_angle, y=OrbHdist)) + stat_density2d(aes(x=AOH_angle, y=OrbHdist, fill = ..density..), geom="raster", contour=F, binwidth = c(0.1, 0.5), xlim =c(-1, 0), ylim = c(0, 3)) + 
			scale_fill_gradient(low = "white", high = "steelblue") +
			scale_x_discrete(paste('Acceptor -- Orbital -- Hydrogen')) +  
			scale_y_discrete(paste('scOrbCpi scHaro Counts < 3.0 A from Orbital -- Hydrogen'))
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
	
})) # end FeaturesAnalysis
