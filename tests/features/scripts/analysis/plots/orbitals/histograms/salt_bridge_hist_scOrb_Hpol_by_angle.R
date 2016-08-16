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
id = "hist_salt_bridges_by_angle",
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
	DHO_angle,
	htype2
FROM
	HPOL_orbital
WHERE
	(((resName1 = 'ARG' OR resName1 = 'LYS' OR resName1 = 'HIS') AND (resName2 = 'ASP' OR resName2 = 'GLU')) OR 
	((resName1 = 'GLU' OR resName1 = 'ASP') AND (resName2 = 'ARG' OR resName2 = 'LYS' OR resName2 = 'HIS'))) 
	AND OrbHdist < 4 AND orbName1 = 'O.p.sp2' AND htype2 = 'Hpol';"


all_geom <- query_sample_sources(sample_sources, sele)



plot_id = "hist_AOH_scOrb_scH_salt_bridge_all_residues"
ggplot(data=all_geom) +
 #geom_freqpoly(aes(x=(acos(AOH_angle))*180/pi, fill=sample_source, color = sample_source), binwidth=5, position="dodge") +
	geom_freqpoly(aes(x=AOH_angle, fill=sample_source, color = sample_source), binwidth=0.1) +
  ggtitle("AOH SaltBridges Combined: His, Arg, Lys, Glu, Asp") +
 # scale_x_continuous(breaks=c(0,30,60,90,120,150, 180), paste('Acceptor -- Orbital -- Hydrogen (degrees)'))
	scale_x_continuous(breaks=c(-1,-.9,-.8,-.7,-.6,-.5,-.4, -.3, -.2, -.1, 0, .1, .2, .3, .4, .5, .6, .7, .8), paste('Acceptor -- Orbital -- Hydrogen (degrees)'))
  scale_y_continuous("Salt Bridge Counts < 4.0 A from Orbital -- Hydrogen")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

plot_id = "hist_DHO_scOrb_scH_salt_bridge_all_residues"
ggplot(data=all_geom) +
  geom_freqpoly(aes(x=DHO_angle, fill=sample_source, color = sample_source), binwidth=0.1) +
  ggtitle("DHO SaltBridges Combined: His, Arg, Lys, Glu, Asp") +
  scale_x_continuous(breaks=c(-1,-.9,-.8,-.7,-.6,-.5,-.4, -.3, -.2, -.1, 0, .1, .2, .3, .4, .5, .6, .7, .8), paste('Donor -- Hydrogen -- Orbital (degrees)'))
  scale_y_continuous("Salt Bridge Counts < 4.0 A from Orbital -- Hydrogen")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)
####Only Histidine
sele <- "
SELECT
  resNum1,
  resName1,
  resNum2,
  resName2,
  OrbHdist,
  AOH_angle,
  DHO_angle,
	htype2
FROM
  HPOL_orbital
WHERE
  (((resName1 = 'HIS') AND (resName2 = 'ASP' OR resName2 = 'GLU')) OR 
  ((resName1 = 'GLU' OR resName1 = 'ASP') AND (resName2 = 'HIS'))) 
  AND OrbHdist < 4 AND orbName1 = 'O.p.sp2' AND htype2 = 'Hpol';"

all_geom <- query_sample_sources(sample_sources, sele)


plot_id = "hist_AOH_scOrb_scH_salt_bridge_His"
ggplot(data=all_geom) +
  geom_freqpoly(aes(x=AOH_angle, fill=sample_source, color = sample_source), binwidth=.1) +
  ggtitle("AOH SaltBridges Combined: His, Glu, Asp") +
  scale_x_continuous(breaks=c(-1,-.9,-.8,-.7,-.6,-.5,-.4, -.3, -.2, -.1, 0, .1, .2, .3, .4, .5, .6, .7, .8), paste('Acceptor -- Orbital -- Hydrogen (degrees)'))
  scale_y_continuous("Salt Bridge Counts < 4.0 A from Orbital -- Hydrogen")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

plot_id = "hist_DHO_scOrb_scH_salt_bridge_His"
ggplot(data=all_geom) +
  geom_freqpoly(aes(x=DHO_angle, fill=sample_source, color = sample_source), binwidth=0.1) +
  ggtitle("DHO SaltBridges Combined: His, Glu, Asp") +
  scale_x_continuous(breaks=c(-1,-.9,-.8,-.7,-.6,-.5,-.4, -.3, -.2, -.1, 0, .1, .2, .3, .4, .5, .6, .7, .8), paste('Donor -- Hydrogen -- Orbital (degrees)'))
  scale_y_continuous("Salt Bridge Counts < 4.0 A from Orbital -- Hydrogen")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

####Only Lys
sele <- "
SELECT
  resNum1,
  resName1,
  resNum2,
  resName2,
  OrbHdist,
  AOH_angle,
  DHO_angle,
	htype2
FROM
  HPOL_orbital
WHERE
  (((resName1 = 'LYS') AND (resName2 = 'ASP' OR resName2 = 'GLU')) OR 
  ((resName1 = 'GLU' OR resName1 = 'ASP') AND (resName2 = 'LYS'))) 
  AND OrbHdist < 4 AND orbName1 = 'O.p.sp2' AND htype2 = 'Hpol';"

all_geom <- query_sample_sources(sample_sources, sele)


plot_id = "hist_AOH_scOrb_scH_salt_bridge_Lys"
ggplot(data=all_geom) +
  geom_freqpoly(aes(x=AOH_angle, fill=sample_source, color = sample_source), binwidth=0.1) +
  ggtitle("AOH SaltBridges Combined: Lys, Glu, Asp") +
  scale_x_continuous(breaks=c(-1,-.9,-.8,-.7,-.6,-.5,-.4, -.3, -.2, -.1, 0, .1, .2, .3, .4, .5, .6, .7, .8), paste('Acceptor -- Orbital -- Hydrogen (degrees)'))
  scale_y_continuous("Salt Bridge Counts < 4.0 A from Orbital -- Hydrogen")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

plot_id = "hist_DHO_scOrb_scH_salt_bridge_Lys"
ggplot(data=all_geom) +
  geom_freqpoly(aes(x=DHO_angle, fill=sample_source, color = sample_source), binwidth=0.1) +
  ggtitle("DHO SaltBridges Combined: Lys, Glu, Asp") +
  scale_x_continuous(breaks=c(-1,-.9,-.8,-.7,-.6,-.5,-.4, -.3, -.2, -.1, 0, .1, .2, .3, .4, .5, .6, .7, .8), paste('Donor -- Hydrogen -- Orbital (degrees)'))
  scale_y_continuous("Salt Bridge Counts < 4.0 A from Orbital -- Hydrogen")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

####Only Lys
sele <- "
SELECT
  resNum1,
  resName1,
  resNum2,
  resName2,
  OrbHdist,
  AOH_angle,
  DHO_angle,
	htype2
FROM
  HPOL_orbital
WHERE
  (((resName1 = 'ARG') AND (resName2 = 'ASP' OR resName2 = 'GLU')) OR 
  ((resName1 = 'GLU' OR resName1 = 'ASP') AND (resName2 = 'ARG'))) 
  AND OrbHdist < 4 AND orbName1 = 'O.p.sp2' AND htype2 = 'Hpol';"

all_geom <- query_sample_sources(sample_sources, sele)


plot_id = "hist_AOH_scOrb_scH_salt_bridge_Arg"
ggplot(data=all_geom) +
  geom_freqpoly(aes(x=AOH_angle, fill=sample_source, color = sample_source), binwidth=0.1) +
  ggtitle("AOH SaltBridges Combined: Arg, Glu, Asp") +
  scale_x_continuous(breaks=c(-1,-.9,-.8,-.7,-.6,-.5,-.4, -.3, -.2, -.1, 0, .1, .2, .3, .4, .5, .6, .7, .8), paste('Acceptor -- Orbital -- Hydrogen (degrees)'))
  scale_y_continuous("Salt Bridge Counts < 4.0 A from Orbital -- Hydrogen")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

plot_id = "hist_DHO_scOrb_scH_salt_bridge_Arg"
ggplot(data=all_geom) +
  geom_freqpoly(aes(x=DHO_angle, fill=sample_source, color = sample_source), binwidth=0.1) +
  ggtitle("DHO SaltBridges Combined: Arg, Glu, Asp") +
  scale_x_continuous(breaks=c(-1,-.9,-.8,-.7,-.6,-.5,-.4, -.3, -.2, -.1, 0, .1, .2, .3, .4, .5, .6, .7, .8), paste('Donor -- Hydrogen -- Orbital (degrees)'))
  scale_y_continuous("Salt Bridge Counts < 4.0 A from Orbital -- Hydrogen")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)




})) # end FeaturesAnalysis
