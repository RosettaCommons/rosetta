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
id = "kmeans_cluster",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <- "
SELECT
  geo.AHdist, geo.cosBAH, geo.cosAHD, geo.chi,
  don.HBChemType AS don_chem_type,
  acc.HBChemType AS acc_chem_type
FROM
  hbonds AS hb,
  hbond_sites AS don,
  hbond_sites AS acc,
  hbond_geom_coords AS geo
WHERE
  don.struct_id = hb.struct_id AND don.site_id = hb.don_id AND
  acc.struct_id = hb.struct_id AND acc.site_id = hb.acc_id AND
  geo.struct_id = hb.struct_id AND geo.hbond_id = hb.hbond_id;"

# Execute the SQL query on each sample source.
f <- query_sample_sources(sample_sources, sele)

f$don_chem_type <- factor(f$don_chem_type,
	levels = c("hbdon_IMD", "hbdon_IME", "hbdon_GDE", "hbdon_GDH",
		"hbdon_AHX", "hbdon_HXL", "hbdon_IND", "hbdon_AMO", "hbdon_CXA", "hbdon_PBA"),
	labels = c("dIMD: h", "dIME: h", "dGDE: r", "dGDH: r",
		"dAHX: y", "dHXL: s,t", "dIND: w", "dAMO: k", "dCXA: n,q", "dPBA: bb"))

f$acc_chem_type <- factor(f$acc_chem_type,
	levels = c("hbacc_IMD", "hbacc_IME", "hbacc_AHX", "hbacc_HXL",
		"hbacc_CXA", "hbacc_CXL", "hbacc_PBA"),
	labels = c("aIMD: h", "aIME: h", "aAHX: y", "aHXL: s,t",
		"aCXA: n,q", "aCXL: d,e", "aPBA: bb"))


#k <- kmeans(f[,c("AHdist", "cosBAH", "cosAHD", "chi")], 2, iter.max=10, nstarts=10)



})) # end FeaturesAnalysis
