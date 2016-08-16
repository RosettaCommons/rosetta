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
id = "sidechain_propensity_to_donate_or_accept_hbonds",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

source("scripts/analysis/plots/hbonds/hbond_geo_dim_scales.R")


sele <- "
SELECT
	don.res_type,
	don.satisfied_label AS don_satisfied,
	acc.satisfied_label AS acc_satisfied,
	CASE WHEN don.buried AND acc.buried THEN 'Buried' ELSE 'Exposed' END AS buried,
	don.HBChemType AS don_chem_type,
	acc.HBChemType AS acc_chem_type,
	count(*) AS count
FROM
	hbond_heavy_sites AS don,
	hbond_heavy_sites AS acc,
	residue_pdb_confidence AS res_b
WHERE
	don.struct_id = acc.struct_id AND
	don.residue_number = acc.residue_number AND
	don.is_donor = 1 AND
	acc.is_donor = 0 AND
	don.struct_id = res_b.struct_id AND don.residue_number = res_b.residue_number AND
	res_b.max_temperature < 30
GROUP BY
	don.res_type, don_satisfied, acc_satisfied, buried, don_chem_type, acc_chem_type;"
f <- query_sample_sources(sample_sources, sele)
f <- f[f$count > 100,]

print(summary(f))
table_id <- "hb_chem_type_don_acc_satisfaction_all"
save_tables(
	self,
	f,
	table_id, sample_sources, output_dir, output_formats)

table_id <- "hb_chem_type_don_acc_satisfaction_sc"
save_tables(
	self,
	f[f$don_chem_type != "hbdon_PBA" & f$acc_chem_type != "hbacc_PBA",],
	table_id, sample_sources, output_dir, output_formats)

table_id <- "hb_chem_type_don_acc_satisfaction_buried_sc"
save_tables(
	self,
	f[f$don_chem_type != "hbdon_PBA" & f$acc_chem_type != "hbacc_PBA" & f$buried == "Buried",],
	table_id, sample_sources, output_dir, output_formats)


plot_parts <- list(
	theme_bw(),
	geom_bar(aes(y=count, fill=sample_source), position="dodge", stat="identity"),
	scale_x_discrete("Residue Type"),
	scale_y_continuous("Residue Counts"))

if(nrow(sample_sources) <= 3){
 plot_parts <- c(plot_parts,
	list(theme(legend.position="bottom", legend.direction="horizontal")))
}


f_sc <- f[f$don_chem_type != "hbdon_PBA" & f$acc_chem_type != "hbacc_PBA",]
f_sc_burial_agnostic <- ddply(f_sc,
	.(res_type, don_satisfied, acc_satisfied, sample_source),
	transform, count=sum(count))

plot_id <- "sidechain_propensity_to_donate_or_accept_hbonds"
p <- ggplot(data=f_sc_burial_agnostic, aes(x=don_satisfied)) + plot_parts +
	ggtitle("Sidechain Propensity to Donate or Accept Hydrogen Bonds (unsat/SAT)") +
	facet_grid(acc_satisfied ~ res_type)
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


#####################################
plot_id <- "sidechain_propensity_to_donate_or_accept_hbonds_buried"
p <- ggplot(data=f_sc[f_sc$buried=="Buried",], aes(x=don_satisfied)) + plot_parts +
	ggtitle("Buried Sidechain Propensity to Donate or Accept Hydrogen Bonds (unsat/SAT)") +
	facet_grid(acc_satisfied ~ res_type)
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

###################################
plot_id <- "sidechain_propensity_to_donate_or_accept_hbonds_exposed"
p <- ggplot(data=f_sc[f_sc$buried=="Exposed",], aes(x=don_satisfied)) + plot_parts +
	ggtitle("Exposed Sidechain Propensity to Donate or Accept Hydrogen Bonds (unsat/SAT)") +
	facet_grid(acc_satisfied ~ res_type)
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

##################################
##f_bbsc <- f[f$don_chem_type == "hbdon_PBA" & f$acc_chem_type == "hbacc_PBA",]
##f_bb$don_acc_satisfied <- paste(f$don_satisfied, f$acc_satisfied, sep=" ")
##
##plot_id <- "backbone_propensity_to_donate_or_accept_hbonds"
##p <- ggplot(data=f_bb, aes(x=don_acc_satisfied)) + theme_bw() +
##	ggtitle("Backbone Propensity to Donate or Accept Hydrogen Bonds (SAT / unsat)") +
##	theme(axis.text.x=theme_text(angle=-90, hjust=0)) +
##	facet_grid(buried ~ res_type)
##save_plots(self, plot_id, sample_sources, output_dir, output_formats)
#
##############
#
#################################
sele <- "
SELECT
	CASE WHEN site.satisfied == 1 THEN 'Sat' ELSE 'UnSat' END AS satisfied,
	site.buried_label AS buried,
	site.HBChemType AS chem_type,
	count(*) AS count
FROM
	residue_pdb_confidence AS res_b
	hbond_heavy_sites AS site,
WHERE
	site.HBChemType != 'hbdon_PBA' AND
	site.HBChemType != 'hbacc_PBA' AND
	site.struct_id = res_b.struct_id AND site.residue_number = res_b.residue_number AND
	res_b.max_temperature < 30
GROUP BY
	satisfied, buried, chem_type;"
f <- query_sample_sources(sample_sources, sele)
f <- f[f$count > 100,]

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
		n_sat = sum(df[df$satisfied == "Sat", "count"]),
		p_sat = sum(df[df$satisfied == "Sat", "count"])/sum(df$count))
})


n_sat_wide <- cast(p_sat, buried + chem_type_name ~ sample_source, value.var=n_sat, value="n_sat")
z <- names(n_sat_wide)
names(n_sat_wide) <- c(
	"buried", "chem_type_name",
	paste("#(Sat) ", z[!(z %in% c("buried", "chem_type_name"))], collapse=""))

p_sat_wide <- cast(p_sat, buried + chem_type_name ~ sample_source, value.var=p_sat, value="p_sat")
z <- names(p_sat_wide)
names(p_sat_wide) <- c(
	"buried", "chem_type_name",
	paste("P(Sat) ", z[!(z %in% c("buried", "chem_type_name"))], collapse=""))


sat_wide <- merge(n_sat_wide, p_sat_wide)
save_tables(self, sat_wide, table_id, sample_sources, output_dir, output_formats)




})) # end FeaturesAnalysis
