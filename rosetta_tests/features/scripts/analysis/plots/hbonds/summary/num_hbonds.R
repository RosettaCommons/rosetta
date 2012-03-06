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
id = "num_hbonds",
author = "Matthew O'Meara",
brief_description = "Count the number of hydrogen bonds formed conditional on the donor and acceptor chemical types.",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self){

source("scripts/analysis/plots/hbonds/hbond_geo_dim_scales.R")

sele <-"
SELECT
  acc_site.HBChemType AS acc_chem_type,
  count(acc_site.HBChemType) AS acc_chem_type_count
FROM
  hbonds AS hbond,
  hbond_sites AS acc_site
WHERE
  hbond.struct_id == acc_site.struct_id AND hbond.acc_id == acc_site.site_id
GROUP BY
  acc_site.HBChemType;"

f <- query_sample_sources(sample_sources, sele)

f$acc_chem_type_name <- acc_chem_type_name_linear(f$acc_chem_type)
f <- na.omit(f, method="r")

plot_id <- "hbonds_num_by_acc_chem_type"
ggplot(f, aes(acc_chem_type_name, log(acc_chem_type_count))) + theme_bw() +
  geom_line(aes(colour=sample_source, group=sample_source), size=.2) +
	geom_point(aes(colour=sample_source), size=3) +
  opts(title = "Number of HBonds by Acceptor Type") +
  labs(x="Acceptor Type", y="log(Counts)")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

sele <-"
SELECT
  don_site.HBChemType AS don_chem_type,
  count(don_site.HBChemType) AS don_chem_type_count
FROM
  hbonds AS hbond,
  hbond_sites AS don_site
WHERE
  hbond.struct_id == don_site.struct_id AND hbond.don_id == don_site.site_id
GROUP BY
  don_site.HBChemType;"

f <- query_sample_sources(sample_sources, sele)

f$don_chem_type_name <- don_chem_type_name_linear(f$don_chem_type)
f <- na.omit(f, method="r")

plot_id <- "hbonds_num_by_don_chem_type"
ggplot(f, aes(don_chem_type_name, log(don_chem_type_count))) + theme_bw() +
  geom_line(aes(colour=sample_source, group=sample_source), size=.2) +
	geom_point(aes(colour=sample_source), size=3) +
  opts(title = "Number of HBonds by Donor Type") +
  labs(x="Donor Type", y="log(Counts)")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)



})) # end FeaturesAnalysis
