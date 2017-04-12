# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

# check_setup()


library(ggplot2)
library(plyr)
library("scales")
library(grid)
library(viridis)

feature_analyses <- c(feature_analyses, new("FeaturesAnalysis",
id = "rama_zoom_helix_stbb_hbonds",
author = "Andrew Leaver-Fay",
brief_description = "Scatter the second sample source on the phi/psi distribution of the first sample source",
feature_reporter_dependencies = c("ResidueFeatures", "ProteinBackboneTorsionAngleFeatures", "ResidueSecondaryStructureFeatures", "PdbDataFeatures", "ProteinResidueConformationFeatures","HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

print("test1")
print(typeof(sample_sources))
print(class(sample_sources))
print(typeof(sample_sources[1]))
print(class(sample_sources[1]))
print("test2")

sele <- "
SELECT
  bb.phi, bb.psi
FROM
  protein_backbone_torsion_angles AS bb
WHERE
	-75 < bb.phi AND bb.phi < -55 AND
	-50 < bb.psi AND bb.psi < -25;"

ss1 = data.frame(sample_sources[1,])
f <- query_sample_sources( ss1, sele)


sele <-"
CREATE TEMPORARY TABLE bif_hbonds_to_bb AS SELECT
  hbond.struct_id AS struct_id,
  hbond.don_id AS don_id,
  hbond.acc_id AS acc_id,
  hbond.energy AS energy
FROM
  hbonds as hbond,
  hbond_sites AS acc_site,
  hbond_site_environment as acc_env
WHERE
  hbond.struct_id = acc_site.struct_id AND hbond.acc_id = acc_site.site_id AND
  acc_site.HBChemType = 'hbacc_PBA' AND
  acc_site.struct_id = acc_env.struct_id AND acc_site.site_id = acc_env.site_id AND
  acc_env.num_hbonds >= 2;

SELECT
	bb.phi,
  bb.psi,
  don_site.resType as don_res_type,
  ( acc_site.resNum - don_site.resNum ) as seq_sep
FROM
  bif_hbonds_to_bb AS hbond,
  bif_hbonds_to_bb AS hbond2,
  hbond_sites AS acc_site,
  hbond_sites AS don_site,
  hbond_sites AS don_site2,
	protein_backbone_torsion_angles AS bb
WHERE
  hbond.struct_id = don_site.struct_id AND hbond.don_id = don_site.site_id AND
  hbond.struct_id = acc_site.struct_id AND hbond.acc_id = acc_site.site_id AND
  hbond2.struct_id = hbond.struct_id AND hbond2.acc_id = hbond.acc_id AND
  hbond2.struct_id = don_site2.struct_id AND hbond2.don_id = don_site2.site_id AND
  don_site.HBChemType = 'hbdon_HXL' AND
  don_site2.HBChemType = 'hbdon_PBA' AND
	bb.struct_id = acc_site.struct_id AND bb.resNum == acc_site.resNum AND
	-75 < bb.phi AND bb.phi < -55 AND
	-50 < bb.psi AND bb.psi < -25
;"
ss2 = data.frame(sample_sources[2,])
g <- query_sample_sources( ss2, sele)


narrow_output_formats <- transform(output_formats, width=height)

# #equal area projection
# f <- transform(f,
# 	capx = 2*sin(acos(cosBAH)/2)*cos(chi),
# 	capy = 2*sin(acos(cosBAH)/2)*sin(chi))
#
# capx_limits <- c(-1.5,1.5)
# capy_limits <- capx_limits


# f_first <- f[ f$sample_source == levels(sample_sources$sample_source)[1], ]
# f_second <- f[ f$sample_source == levels(sample_sources$sample_source)[2], ]
#
#f_first <- ddply(
#  f_first, .(sample_source), transform, counts = length(sample_source))

# plot_id = "chi_BAH_eqpoldens_and_scatter_lrbb"
# ggplot(data=f_first) + theme_bw() +
# 	theme(panel.background=element_rect(fill="#00007F", colour="#00007F")) +
# 	stat_density2d(
# 		aes(x=capx,y=capy, fill=..density..), geom="tile", contour=FALSE ) +
# 	polar_equal_area_grids_bw() +
# 	geom_indicator(aes(indicator=counts), color="white") +
# 	geom_point( data=f_second,aes(x=capx,y=capy),colour="white",size=2) +
# 	ggtitle(
# 		paste("Hydrogen Bonds chi vs BAH Angles with Sequence Separation > 5\n",
# 		"Backbone/Backbone Hydrogen Bonds, Equal Coordinate Projection\n",
# 		"Reference (density) vs Test (white circles)", sep="")) +
# 	scale_x_continuous(
# 		'2*sin(BAH/2) * cos(CHI)', limits=capx_limits, breaks=c(-1, 0, 1)) +
# 	scale_y_continuous(
# 		'2*sin(BAH/2) * sin(CHI)', limits=capy_limits, breaks=c(-1, 0, 1)) +
# 	coord_fixed(ratio = 1) +
# 	scale_fill_gradientn('Density', colours=jet.colors(10))
# save_plots(self, plot_id, sample_sources, output_dir, output_formats)

dens <- MASS::kde2d(f$phi, f$psi, n=300, h=1)
densdf <- data.frame(expand.grid(x=dens$x, y=dens$y), z=as.vector(dens$z))
densdf$counts <- nrow(f)

# gdens <- MASS::kde2d(g$phi, g$psi, n=300, h=1)
# gdensdf <- data.frame(expand.grid(g=gdens$x, g=gdens$y), z=as.vector(gdens$z))
# gdensdf$counts <- nrow(g)

#jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
#  "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

ggplot(data=densdf) +
	theme_bw() +
	theme(panel.background=element_rect(fill="#00007F", colour="#00007F")) +
	geom_tile(aes(x=x, y=y, fill=z)) +
	geom_indicator(aes(indicator=counts), color="white") +
	geom_point( data=g,aes(x=phi,y=psi),colour="white",size=0.5 ) +
	geom_density_2d( data=g,aes(x=phi,y=psi),colour="white") +
	coord_equal(ratio=1) +
	scale_fill_viridis('Density', limits=c(0, .01)) +
	theme(legend.position="bottom", legend.direction="horizontal") +
	ggtitle(paste(
		"Ser/Thr to BB hbonds with competing BB/BB hbonds.\n Torsion Angles: ", sample_sources[2,2], sep="")) +
	scale_x_continuous(expression(paste("phi Angle (Degrees)", sep="")), limits=c(-75, -55)) +
	scale_y_continuous(expression(paste("psi Angle (Degrees)", sep="")), limits=c(-50, -25))

plot_id <- "rama_for_helix_w_ST_bb_hbonds"
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
