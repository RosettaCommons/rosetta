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
id = "pi_pi_eqpoldens_scOrbCpi_Haro_chiBAHD_cosAHD",
author = "Matthew O'Meara, Steven Combs",
brief_description = "",
feature_reporter_dependencies = c("OrbitalFeatures"),
run=function(self, sample_sources, output_dir, output_formats){




#####Pi Stacking interactions 2.0

sele <- "
SELECT
  resNum1,
  resName1,
  resNum2,
  resName2,
	OrbName1,
	cosAHD,
  OrbHdist,
	chiBAHD,
	htype2
FROM
  Haro_orbital
WHERE
  (((resName1 = 'PHE' OR resName1 = 'TYR' OR resName1 = 'TRP' ) AND (resName2 = 'PHE' OR resName2 = 'TYR' OR resName2 = 'TRP')))  AND 
	OrbHdist < 2 AND 
	OrbName1 = 'C.pi.sp2' AND
	ABS(resNum1 - resNum2) > 5;"

f <- query_sample_sources(sample_sources, sele)

f <- transform(f,
  capx = 2*sin(acos(cosAHD)/2)*cos(chiBAHD),
  capy = 2*sin(acos(cosAHD)/2)*sin(chiBAHD))

capx_limits <- c(-1.5,1.5)
capy_limits <- capx_limits


plot_id = "pi_pi_eqpoldens_scOrbCpi_Haro_chiBAHD_cosAHD_all_2.0A"

f_first <- f[ f$sample_source == levels(sample_sources$sample_source), ]

ggplot(data=f_first) + theme_bw() +
		theme(panel.background=element_rect(fill="#00007F", colour="#00007F")) +
		stat_density2d(
				aes(x=capx, y=capy, fill=..density..), geom="tile", contour=FALSE ) +
		polar_equal_area_grids_bw() +
		ggtitle(
						paste("Pi Stacking chiBAHD vs AHD Angles with Sequence Separation > 5\n",
    "Sidechain Donors to Sidechain sp2 Acceptors, Equal Coordinate Projection\n",
		"PHE+TYR+TRP to PHE+TYR+TRP at 2.0A", sep="")) +
		scale_x_continuous(
				'2*sin(AHD/2) * cos(chiBAHD)', limits=capx_limits, breaks=c(-1, 0, 1)) +
		scale_y_continuous(
				'2*sin(AHD/2) * sin(chiBAHD)', limits=capy_limits, breaks=c(-1, 0, 1)) +
		coord_fixed(ratio = 1) +
		scale_fill_gradientn('Density', colours=jet.colors(10))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


#####Pi Stacking interactions 3.5

sele <- "
SELECT
  resNum1,
  resName1,
  resNum2,
  resName2,
  OrbName1,
  cosAHD,
  OrbHdist,
  chiBAHD,
  htype2
FROM
  Haro_orbital
WHERE
  (((resName1 = 'PHE' OR resName1 = 'TYR' OR resName1 = 'TRP' ) AND (resName2 = 'PHE' OR resName2 = 'TYR' OR resName2 = 'TRP')))  AND 
  OrbHdist < 3.5 AND 
  OrbName1 = 'C.pi.sp2' AND
  ABS(resNum1 - resNum2) > 5;"

f <- query_sample_sources(sample_sources, sele)

f <- transform(f,
  capx = 2*sin(acos(cosAHD)/2)*cos(chiBAHD),
  capy = 2*sin(acos(cosAHD)/2)*sin(chiBAHD))

capx_limits <- c(-1.5,1.5)
capy_limits <- capx_limits


plot_id = "pi_pi_eqpoldens_scOrbCpi_Haro_chiBAHD_cosAHD_all_3.5A"

f_first <- f[ f$sample_source == levels(sample_sources$sample_source), ]


ggplot(data=f_first) + theme_bw() +
		theme(panel.background=element_rect(fill="#00007F", colour="#00007F")) +
		stat_density2d(
				aes(x=capx, y=capy, fill=..density..), geom="tile", contour=FALSE ) +
		polar_equal_area_grids_bw() +
		ggtitle(
						paste("Pi Stacking chiBAHD vs AHD Angles with Sequence Separation > 5\n",
    "Sidechain Donors to Sidechain sp2 Acceptors, Equal Coordinate Projection\n",
    "PHE+TYR+TRP to PHE+TYR+TRP at 3.5A", sep="")) +
		scale_x_continuous(
				'2*sin(AHD/2) * cos(chiBAHD)', limits=capx_limits, breaks=c(-1, 0, 1)) +
		scale_y_continuous(
				'2*sin(AHD/2) * sin(chiBAHD)', limits=capy_limits, breaks=c(-1, 0, 1)) +
		coord_fixed(ratio = 1) +
		scale_fill_gradientn('Density', colours=jet.colors(10))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

###Phe-Phe 2.0
sele <- "
SELECT
  resNum1,
  resName1,
  resNum2,
  resName2,
  OrbName1,
  cosAHD,
  OrbHdist,
  chiBAHD,
  htype2
FROM
  Haro_orbital
WHERE
  (((resName1 = 'PHE' ) AND (resName2 = 'PHE' )))  AND 
  OrbHdist < 2 AND 
  OrbName1 = 'C.pi.sp2' AND
  ABS(resNum1 - resNum2) > 5;"

f <- query_sample_sources(sample_sources, sele)

f <- transform(f,
  capx = 2*sin(acos(cosAHD)/2)*cos(chiBAHD),
  capy = 2*sin(acos(cosAHD)/2)*sin(chiBAHD))

capx_limits <- c(-1.5,1.5)
capy_limits <- capx_limits


plot_id = "pi_pi_eqpoldens_scOrbCpi_Haro_chiBAHD_cosAHD_phe_phe_2.0A"

f_first <- f[ f$sample_source == levels(sample_sources$sample_source), ]

ggplot(data=f_first) + theme_bw() +
		theme(panel.background=element_rect(fill="#00007F", colour="#00007F")) +
		stat_density2d(
				aes(x=capx, y=capy, fill=..density..), geom="tile", contour=FALSE ) +
		polar_equal_area_grids_bw() +
		ggtitle(
						paste("Pi Stacking chiBAHD vs AHD Angles with Sequence Separation > 5\n",
    "Sidechain Donors to Sidechain sp2 Acceptors, Equal Coordinate Projection\n",
    "PHE to PHE at 2.0A", sep="")) +
		scale_x_continuous(
				'2*sin(AHD/2) * cos(chiBAHD)', limits=capx_limits, breaks=c(-1, 0, 1)) +
		scale_y_continuous(
				'2*sin(AHD/2) * sin(chiBAHD)', limits=capy_limits, breaks=c(-1, 0, 1)) +
		coord_fixed(ratio = 1) +
		scale_fill_gradientn('Density', colours=jet.colors(10))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)




###Phe -Phe 3.5
sele <- "
SELECT
  resNum1,
  resName1,
  resNum2,
  resName2,
  OrbName1,
  cosAHD,
  OrbHdist,
  chiBAHD,
  htype2
FROM
  Haro_orbital
WHERE
  (((resName1 = 'PHE' ) AND (resName2 = 'PHE' )))  AND 
  OrbHdist < 3.5 AND 
  OrbName1 = 'C.pi.sp2' AND
  ABS(resNum1 - resNum2) > 5;"

f <- query_sample_sources(sample_sources, sele)

f <- transform(f,
  capx = 2*sin(acos(cosAHD)/2)*cos(chiBAHD),
  capy = 2*sin(acos(cosAHD)/2)*sin(chiBAHD))

capx_limits <- c(-1.5,1.5)
capy_limits <- capx_limits


plot_id = "pi_pi_eqpoldens_scOrbCpi_Haro_chiBAHD_cosAHD_phe_phe_3.5A"

f_first <- f[ f$sample_source == levels(sample_sources$sample_source), ]


ggplot(data=f_first) + theme_bw() +
		theme(panel.background=element_rect(fill="#00007F", colour="#00007F")) +
		stat_density2d(
				aes(x=capx, y=capy, fill=..density..), geom="tile", contour=FALSE ) +
		polar_equal_area_grids_bw() +
		ggtitle(
						paste("Pi Stacking chiBAHD vs AHD Angles with Sequence Separation > 5\n",
    "Sidechain Donors to Sidechain sp2 Acceptors, Equal Coordinate Projection\n",
    "PHE to PHE at 3.5A", sep="")) +
		scale_x_continuous(
				'2*sin(AHD/2) * cos(chiBAHD)', limits=capx_limits, breaks=c(-1, 0, 1)) +
		scale_y_continuous(
				'2*sin(AHD/2) * sin(chiBAHD)', limits=capy_limits, breaks=c(-1, 0, 1)) +
		coord_fixed(ratio = 1) +
		scale_fill_gradientn('Density', colours=jet.colors(10))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

})) # end FeaturesAnalysis
