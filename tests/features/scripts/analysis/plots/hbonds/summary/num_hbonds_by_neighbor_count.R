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
id = "num_hbonds_by_neighbor_count",
author = "Matthew O'Meara",
brief_description = "Count the number of hydrogen bonds formed conditional on the donor and acceptor chemical types.",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

source("scripts/analysis/plots/hbonds/hbond_geo_dim_scales.R")

sele <-"
SELECT
  acc_site.HBChemType AS chem_type,
  count(acc_site.HBChemType) AS chem_type_count,
	CASE
		WHEN burial.ten_a_neighbors < 17 THEN 'exposed'
		WHEN burial.ten_a_neighbors BETWEEN 17 AND 23 THEN 'partial'
		WHEN 24 <= burial.ten_a_neighbors THEN 'buried' END AS solvent_exposure,
	'Acceptor' AS site_type
FROM
  hbonds AS hbond,
  hbond_sites AS acc_site,
	residue_burial AS burial
WHERE
  hbond.struct_id == acc_site.struct_id AND hbond.acc_id == acc_site.site_id AND
	burial.struct_id == acc_site.struct_id AND burial.resNum == acc_site.resNum
GROUP BY
  acc_site.HBChemType,
	CASE
		WHEN burial.ten_a_neighbors < 17 THEN 'exposed'
		WHEN burial.ten_a_neighbors BETWEEN 17 AND 23 THEN 'partial'
		WHEN 24 <= burial.ten_a_neighbors THEN 'buried' END
UNION
SELECT
  don_site.HBChemType AS chem_type,
  count(don_site.HBChemType) AS chem_type_count,
	CASE
		WHEN burial.ten_a_neighbors < 17 THEN 'exposed'
		WHEN burial.ten_a_neighbors BETWEEN 17 AND 23 THEN 'partial'
		WHEN 24 <= burial.ten_a_neighbors THEN 'buried' END AS solvent_exposure,
	'Donor' AS site_type
FROM
  hbonds AS hbond,
  hbond_sites AS don_site,
	residue_burial AS burial
WHERE
  hbond.struct_id == don_site.struct_id AND hbond.don_id == don_site.site_id AND
	burial.struct_id == don_site.struct_id AND burial.resNum == don_site.resNum
GROUP BY
  don_site.HBChemType,
	CASE
		WHEN burial.ten_a_neighbors < 17 THEN 'exposed'
		WHEN burial.ten_a_neighbors BETWEEN 17 AND 23 THEN 'partial'
		WHEN 24 <= burial.ten_a_neighbors THEN 'buried' END;"

f <- query_sample_sources(sample_sources, sele)

f$solvent_exposure <- factor(f$solvent_exposure,
	levels = c("exposed", "partial", "buried"),
	labels = c("Exposed (#Nbrs < 17)", "Partial (17 <= #Nbrs < 24)", "Buried (24 <= #Nbrs)"))

f$chem_type_name <- chem_type_name_linear(f$chem_type)
f <- na.omit(f, method="r")
f <- f[order(f$chem_type_name),]

plot_id <- "hbonds_num_hbonds_by_neighbor_count"
p <- ggplot(f, aes(log(chem_type_count), chem_type_name)) + theme_bw() +
  geom_path(aes(colour=sample_source, group=sample_source), size=.2) +
	geom_point(aes(colour=sample_source), size=3) +
	facet_grid(site_type ~ solvent_exposure, scales="free_y", space="free") +
  ggtitle("Number of HBonds by Neighbor Count") +
  labs(y="Chemical Type of Hydrogen Bond Site", x="log(Counts)")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

})) # end FeaturesAnalysis
