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
id = "pi_pi_hist_ScOrbCpi_Haro_by_angle",
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


all_geom <- query_sample_sources(sample_sources, sele)
#ss_id <- all_geom$sample_sources[1]
#plot_id = paste("hist_AOH_scOrbPi_scHaro_all_residues", ss_id, sep="_")
plot_id = "hist_AOH_scOrbPi_scHaro_all_residues"
ggplot(data=all_geom) +
		geom_freqpoly(aes(x=AOH_angle, fill=sample_source, color = sample_source), binwidth=0.1) +
  ggtitle("AOH Pi-Pi measured via Orb - Haro \n", 
		  "ScOrbCpi to Haro Combined: TYR, PHE, TRP") +
  scale_x_continuous(breaks=c(-1,-.9,-.8,-.7,-.6,-.5,-.4, -.3, -.2, -.1, 0, .1, .2, .3, .4, .5, .6, .7, .8), paste('Acceptor -- Orbital -- Hydrogen (degrees)'))
  scale_y_continuous("scOrbCpi scHaro Counts < 4.0 A from Orbital -- Hydrogen")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

plot_id = "hist_DHO_scOrbPi_scHaro_all_residues"
ggplot(data=all_geom) +
		geom_freqpoly(aes(x=DHO_angle, fill=sample_source, color = sample_source), binwidth=0.1) +
  ggtitle("DHO Pi-Pi measured via Orb - Haro \n", 
		  "ScOrbCpi to Haro Combined: TYR, PHE, TRP") +
  scale_x_continuous(breaks=c(-1,-.9,-.8,-.7,-.6,-.5,-.4, -.3, -.2, -.1, 0, .1, .2, .3, .4, .5, .6, .7, .8), paste('Donor -- Hydrogen -- Orbital (degrees)'))
  scale_y_continuous("scOrbCpi scHaro Counts < 4.0 A from Orbital -- Hydrogen")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

#####Phe - Phe

#sele <- "
#SELECT
#  resNum1,
#  resName1,
#  resNum2,
#  resName2,
#  OrbHdist,
#  AOH_angle,
#  DHO_angle
#FROM
#  HARO_orbital
#WHERE
#  (((resName1 = 'PHE') AND (resName2 = 'PHE'))) 
#  AND OrbHdist < 4;"
#
#
#all_geom <- query_sample_sources(sample_sources, sele)
#
##plot_id = paste("hist_AOH_scOrbPi_scHaro_Phe_Phe", ss_id, sep="_")
#
#plot_id = "hist_AOH_scOrbPi_scHaro_Phe_Phe"
#ggplot(data=all_geom) +
#  geom_bar(aes(x=AOH_angle, fill=sample_source), binwidth=0.1, position="dodge") +
#  ggtitle("AOH ScOrbCpi to Haro Phe to Phe") +
#  scale_x_continuous(breaks=c(-1,-.9,-.8,-.7,-.6,-.5,-.4, -.3, -.2, -.1, 0, .1, .2, .3, .4, .5, .6, .7, .8), paste('Acceptor -- Orbital -- Hydrogen (degrees)'))
#  scale_y_continuous("scOrbCpi scHaro Counts < 4.0 A from Orbital -- Hydrogen")
#save_plots(self, plot_id, sample_sources, output_dir, output_formats)
#
#plot_id = "hist_DHO_scOrbPi_scHaro_Phe_Phe"
##plot_id = paste("hist_DHO_scOrbPi_scHaro_Phe_Phe", ss_id, sep="_")
#ggplot(data=all_geom) +
#  geom_bar(aes(x=DHO_angle, fill=sample_source), binwidth=0.1, position="dodge") +
#  ggtitle("DHO ScOrbCpi to Haro Phe to Phe") +
#  scale_x_continuous(breaks=c(-1,-.9,-.8,-.7,-.6,-.5,-.4, -.3, -.2, -.1, 0, .1, .2, .3, .4, .5, .6, .7, .8), paste('Donor -- Hydrogen -- Orbital (degrees)'))
#  scale_y_continuous("scOrbCpi scHaro Counts < 4.0 A from Orbital -- Hydrogen")
#save_plots(self, plot_id, sample_sources, output_dir, output_formats)




})) # end FeaturesAnalysis
