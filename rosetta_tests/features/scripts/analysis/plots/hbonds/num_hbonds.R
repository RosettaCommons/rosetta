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
filename = "scripts/analysis/plots/hbonds/num_hbonds.R",
author = "Matthew O'Meara",
brief_description = "Count the number of hydrogen bonds formed conditional on the donor and acceptor chemical types.",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self){

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

# This is deprecated please use the hbond_chem_types table for the lables instead
# Order the plots better and give more descriptive labels
f$acc_chem_type <- factor(f$acc_chem_type,
	levels = c("hbacc_IMD", "hbacc_IME", "hbacc_AHX", "hbacc_HXL",
		"hbacc_CXA", "hbacc_CXL", "hbacc_PBA"),
	labels = c("aIMD: h", "aIME: h", "aAHX: y", "aHXL: s,t",
		"aCXA: n,q", "aCXL: d,e", "aPBA: bb"))


plot_id <- "hbonds_num_by_acc_chem_type"
ggplot(f, aes(acc_chem_type, log(acc_chem_type_count))) + theme_bw() +
  geom_bar(aes(fill=sample_source), position="dodge") +
  opts(title = "HBond Counts for each Acceptor Type") +
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

# This is deprecated please use the hbond_chem_types table for the lables instead
# Order the plots better and give more descriptive labels
f$don_chem_type <- factor(f$don_chem_type,
	levels = c("hbdon_IMD", "hbdon_IME", "hbdon_GDE", "hbdon_GDH",
		"hbdon_AHX", "hbdon_HXL", "hbdon_IND", "hbdon_AMO", "hbdon_CXA", "hbdon_PBA"),
	labels = c("dIMD: h", "dIME: h", "dGDE: r", "dGDH: r",
		"dAHX: y", "dHXL: s,t", "dIND: w", "dAMO: k", "dCXA: n,q", "dPBA: bb"))

plot_id <- "hbonds_num_by_don_chem_type"
ggplot(f, aes(don_chem_type, log(don_chem_type_count))) + theme_bw() +
  geom_bar(aes(fill=sample_source), position="dodge") +
  opts(title = "HBond Counts for Each Donor Type") +
  labs(x="Donor Type", y="log(Counts)")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)



})) # end FeaturesAnalysis
