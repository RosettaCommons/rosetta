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
id = "num_hbonds",
author = "Matthew O'Meara",
brief_description = "Count the number of hydrogen bonds formed conditional on the donor and acceptor chemical types.",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

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
  ggtitle("Number of HBonds by Acceptor Type") +
  labs(x="Acceptor Type", y="log(Counts)")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

t <- cast(f, sample_source ~ acc_chem_type_name, value="acc_chem_type_count")
save_tables(self, t, plot_id, sample_sources, output_dir, output_formats, caption="Number of HBonds by Acceptor Type", caption.placement="top")

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
  ggtitle("Number of HBonds by Donor Type") +
  labs(x="Donor Type", y="log(Counts)")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

t <- cast(f, sample_source ~ don_chem_type_name, value="don_chem_type_count")
save_tables(self, t, plot_id, sample_sources, output_dir, output_formats, caption="Number of HBonds by Donor Type", caption.placement="top")


sele <-"
SELECT
	don_counts.num_don_sites,
	acc_counts.num_acc_sites,
  bond_counts.don_chem_type,
  bond_counts.acc_chem_type,
	bond_counts.bond_counts
FROM
	(SELECT
	  don_site.HBChemType AS don_chem_type,
	  acc_site.HBChemType AS acc_chem_type,
	  count(*) AS bond_counts
	FROM
	  hbonds AS hbond,
	  hbond_sites AS don_site,
	  hbond_sites AS acc_site
	WHERE
	  hbond.struct_id == don_site.struct_id AND hbond.don_id == don_site.site_id AND
	  hbond.struct_id == acc_site.struct_id AND hbond.acc_id == acc_site.site_id
	GROUP BY
	  don_site.HBChemType,
	  acc_site.HBChemType) AS bond_counts,
	(SELECT
		don.HBChemType AS don_chem_type,
		count(*) AS num_don_sites
	FROM
		hbond_sites AS don
	GROUP BY
		don_chem_type) AS don_counts,
	(SELECT
		acc.HBChemType AS acc_chem_type,
		count(*) AS num_acc_sites
	FROM
		hbond_sites AS acc
	GROUP BY
		acc_chem_type) AS acc_counts
WHERE
	bond_counts.don_chem_type = don_counts.don_chem_type AND
	bond_counts.acc_chem_type = acc_counts.acc_chem_type;"

f <- query_sample_sources(sample_sources, sele)

f$don_chem_type_name <- don_chem_type_name_linear(f$don_chem_type)
f$acc_chem_type_name <- acc_chem_type_name_linear(f$acc_chem_type)
f <- na.omit(f, method="r")

plot_id <- "hbonds_num_don_acc_by_chem_type"
ggplot(f, aes(don_chem_type_name, log(bond_counts))) + theme_bw() +
  geom_line(aes(colour=sample_source, group=sample_source), size=.2) +
	geom_point(aes(colour=sample_source), size=2) +
	facet_wrap( ~ acc_chem_type_name) +
  ggtitle("Number of HBonds by Donor and Acceptor Type") +
  labs(x="Donor Type", y="log(Counts)") +
	coord_flip()

save_plots(self, plot_id, sample_sources, output_dir, output_formats)

plot_id <- "hbonds_num_acc_don_by_chem_type"
ggplot(f, aes(acc_chem_type_name, log(bond_counts))) + theme_bw() +
  geom_line(aes(colour=sample_source, group=sample_source), size=.2) +
	geom_point(aes(colour=sample_source), size=2) +
	facet_wrap( ~ don_chem_type_name) +
  ggtitle("Number of HBonds by Acceptor and Donor Type") +
  labs(x="Acceptor Type", y="log(Counts)") +
	coord_flip() +

save_plots(self, plot_id, sample_sources, output_dir, output_formats)

table_id <- "hbond_num_bonds_by_don_acc_chem_type"
t <- cast(f,
	sample_source +
	acc_chem_type_name ~ don_chem_type_name,
	value="bond_counts")
save_tables(self, t, plot_id, sample_sources, output_dir, output_formats, caption="Number of HBonds by Acceptor and Donor Type", caption.placement="top")


})) # end FeaturesAnalysis
