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
id = "cosAHD_chem_type",
filename = "scripts/analysis/plots/hbonds/cosAHD_chem_type.R",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self){

sele <-"
SELECT
  geom.cosAHD,
  acc_site.HBChemType AS acc_chem_type,
  don_site.HBChemType AS don_chem_type
FROM
  hbond_geom_coords AS geom,
  hbonds AS hbond,
  hbond_sites AS don_site,
  hbond_sites AS acc_site
WHERE
  hbond.struct_id = geom.struct_id AND
  hbond.hbond_id =  geom.hbond_id AND
  hbond.struct_id = don_site.struct_id AND
  hbond.don_id = don_site.site_id AND
  hbond.struct_id = acc_site.struct_id AND
  hbond.acc_id = acc_site.site_id AND
  abs( don_site.resNum - acc_site.resNum ) > 5;"

f <- query_sample_sources(sample_sources, sele)

# This is deprecated please use the hbond_chem_types table for the lables instead
# Order the plots better and give more descriptive labels
f$don_chem_type <- factor(f$don_chem_type,
	levels = c("hbdon_IMD", "hbdon_IME", "hbdon_GDE", "hbdon_GDH",
		"hbdon_AHX", "hbdon_HXL", "hbdon_IND", "hbdon_AMO", "hbdon_CXA", "hbdon_PBA"),
	labels = c("dIMD: h", "dIME: h", "dGDE: r", "dGDH: r",
		"dAHX: y", "dHXL: s,t", "dIND: w", "dAMO: k", "dCXA: n,q", "dPBA: bb"))

# This is deprecated please use the hbond_chem_types table for the lables instead
# Order the plots better and give more descriptive labels
f$acc_chem_type <- factor(f$acc_chem_type,
	levels = c("hbacc_IMD", "hbacc_IME", "hbacc_AHX", "hbacc_HXL",
		"hbacc_CXA", "hbacc_CXL", "hbacc_PBA"),
	labels = c("aIMD: h", "aIME: h", "aAHX: y", "aHXL: s,t",
		"aCXA: n,q", "aCXL: d,e", "aPBA: bb"))

#coAHD goes from 0 to 1, where 1 is linear
#since there is significant density at 1,
#to accurately model a discontinuity, reflect
#the data across the right boundary, in computing the density esitmation
dens <- estimate_density_1d_reflect_boundary(
 data=f,
 ids = c("sample_source", "acc_chem_type", "don_chem_type"),
 variable = "cosAHD",
 reflect_right=TRUE,
 right_boundary=1)

plot_id = "cosAHD_chem_type"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=180-acos(x)*180/pi, y=y, colour=sample_source)) +
	geom_indicator(aes(colour=sample_source, indicator=counts)) +
	facet_grid(don_chem_type ~ acc_chem_type) +
	opts(title = "Hydrogen Bonds AHD Angle by Chemical Type; SeqSep > 5\n(normalized for equal volume per unit distance)") +
	scale_y_continuous("FeatureDensity", limits=c(0,20), breaks=c(0,5,10,15)) +
	scale_x_continuous("Acceptor -- Hydrogen -- Donor (degrees)", trans="reverse")

if(nrow(sample_sources) <= 3){
	p <- p + opts(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
