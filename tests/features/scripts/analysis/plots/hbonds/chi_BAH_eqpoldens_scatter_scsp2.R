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
id = "chi_BAH_eqpoldens_scatter_scsp2",
author = "Matthew O'Meara",
brief_description = "The rank of a hydrogen bond at donor site or acceptor site is rank of the relative Rosetta HBond energy of the hydrogen bond at the site.",

long_description = "
The plot is of hydrogen bonds with sidechain donor groups and Sp2 sidechain acceptors (eg. in proteins, Asp, Glu, Asn or Gln), filtered for so the sequence separation between the residues of the donor and acceptor groups is greater than 5.

Dimensions of the plot is the direction of the Hydrogen atom from the Acceptor atom. This direction is projected via the Lambert Azmuthal projection with the center point where the Acceptor Base, the Acceptor and the Hydrogen are colinear, the positive X axis is the Anti orbital and the negative X axis is the Syn orbital.T he Syn orbital is in the Acceptor Sp2 plane the Acceptor Base 2, Acceptor Base, Acceptor, and Orbital positions form a trans dihedral angle.

The plot organization is a smooth 2D density for the first sample source and then put points on to it from the second sample source. This is useful, for example, when one wants to determine the likelihood that a few samples (the first sample source) came from the a larger sample source (the second sample source).",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
	geom.cosBAH, geom.chi,
	acc.HBChemType AS acc_chem_type, don.HBChemType AS don_chem_type
FROM
	hbond_geom_coords AS geom,
	hbonds AS hbond,
	hbond_sites AS don, hbond_sites AS acc
WHERE
  don.HBChemType != 'hbdon_PBA' AND
	(acc.HBChemType == 'hbacc_CXL' OR acc.HBChemType == 'hbacc_CXA') AND
	hbond.struct_id = geom.struct_id AND hbond.hbond_id = geom.hbond_id AND
	hbond.struct_id = don.struct_id AND hbond.don_id = don.site_id AND
	hbond.struct_id = acc.struct_id AND hbond.acc_id = acc.site_id AND
	ABS(don.resNum - acc.resNum) > 5;";
f <- query_sample_sources(sample_sources, sele)

#equal area projection
f <- transform(f,
	capx = 2*sin(acos(cosBAH)/2)*cos(chi),
	capy = 2*sin(acos(cosBAH)/2)*sin(chi))

capx_limits <- c(-1.5,1.5)
capy_limits <- capx_limits

plot_id = "chi_BAH_eqpoldens_and_scatter_scsc_to_sp2"

f_first <- f[ f$sample_source == levels(sample_sources$sample_source)[1], ]
f_second <- f[ f$sample_source == levels(sample_sources$sample_source)[2], ]

ggplot(data=f_first) + theme_bw() +
	theme(panel.background=element_rect(fill="#00007F", colour="#00007F")) +
  stat_density2d(
		aes(x=capx, y=capy, fill=..density..), geom="tile", contour=FALSE ) +
	polar_equal_area_grids_bw() +
	geom_point(data=f_second,aes(x=capx,y=capy),colour="white",size=2) +
	theme(title =
		paste("Hydrogen Bonds chi vs BAH Angles with Sequence Separation > 5\n",
		"Sidechain Donors to Sidechain sp2 Acceptors, Equal Coordinate Projection\n",
		"Reference (density) vs Test (white circles)", sep="")) +
	scale_x_continuous(
		'2*sin(BAH/2) * cos(CHI)', limits=capx_limits, breaks=c(-1, 0, 1)) +
	scale_y_continuous(
		'2*sin(BAH/2) * sin(CHI)', limits=capy_limits, breaks=c(-1, 0, 1)) +
	coord_fixed(ratio = 1) +
	scale_fill_gradientn('Density', colours=jet.colors(10))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
