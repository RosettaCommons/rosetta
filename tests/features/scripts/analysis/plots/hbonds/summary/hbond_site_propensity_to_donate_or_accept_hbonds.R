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
id = "hbond_site_propensity_to_donate_or_accept_hbonds",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

source("scripts/analysis/plots/hbonds/hbond_geo_dim_scales.R")


sele <- "
SELECT
	CASE WHEN site.satisfied == 1 THEN 'Sat' ELSE 'UnSat' END AS satisfied,
	site.buried_label AS buried,
	site.HBChemType AS chem_type,
	count(*) AS count
FROM
	hbond_heavy_sites AS site
WHERE
	site.HBChemType != 'hbdon_PBA' AND
	site.HBChemType != 'hbacc_PBA'
GROUP BY
	satisfied, buried, chem_type;"
f <- query_sample_sources(sample_sources, sele)
f <- f[f$count > 100,]


plot_parts <- list(
	theme_bw(),
	geom_bar(aes(y=count, fill=sample_source), position="dodge", stat="identity"),
	scale_x_discrete("Residue Type"),
	scale_y_continuous("Residue Counts"))

if(nrow(sample_sources) <= 3){
 plot_parts <- c(plot_parts,
	list(theme(legend.position="bottom", legend.direction="horizontal")))
}

f$chem_type_name <- chem_type_name_linear(f$chem_type)

plot_id <- "hbond_site_propensity_to_hydrogen_bond"
p <- ggplot(data=f, aes(x=satisfied)) + plot_parts +
	ggtitle("H-Bond Site Propensity to Hydrogen Bond") +
	theme(axis.text.x=theme_text(angle=-90, hjust=0)) +
	facet_grid(buried ~ chem_type_name)
save_plots(self, plot_id, sample_sources, output_dir, output_formats)



#####################################
table_id <- "hb_chem_type_fraction_satisfied"
p_sat <- ddply(f, .(chem_type_name, buried, sample_source), function(df){
	data.frame(
		p_sat = sum(df[df$satisfied == "Sat", "count"])/sum(df$count))
})
p_sat_wide <- cast(p_sat, buried + chem_type_name ~ sample_source, value.var=p_sat, value="p_sat")
save_tables(self, p_sat_wide, table_id, sample_sources, output_dir, output_formats)




})) # end FeaturesAnalysis
