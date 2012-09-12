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
id = "sidechain_propensity_to_donate_or_accept_hbonds",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

source("scripts/analysis/plots/hbonds/hbond_geo_dim_scales.R")

sele <- "
DROP TABLE IF EXISTS hbond_heavy_sites;
CREATE TABLE IF NOT EXISTS hbond_heavy_sites (
	site_id INTEGER PRIMARY KEY AUTOINCREMENT,
	struct_id BLOB,
	residue_number INTEGER,
	HBChemType TEXT,
	is_donor INTEGER,
	res_type TEXT,
	satisfied_label TEXT,
	satisfied INTEGER,
	buried_label TEXT,
	buried INTEGER);

INSERT INTO hbond_heavy_sites
SELECT
	NULL,
	site.struct_id,
	site.resNum,
	site.HBChemType,
	site.is_donor,
	res.name3 AS res_type,
	CASE WHEN env.num_hbonds == 0 THEN 'acc' ELSE 'ACC' END AS satisfied_label,
	env.num_hbonds != 0 AS satisfied,
	CASE WHEN env.sasa_r100 == 0 THEN 'Buried' ELSE 'Exposed' END AS buried_label,
	env.sasa_r100 == 0 AS buried
FROM
	residues AS res,
	hbond_sites AS site,
	hbond_site_environment AS env
WHERE
	site.struct_id = res.struct_id AND
	site.resNum = res.resNum AND
	env.struct_id = site.struct_id AND
	env.site_id = site.site_id AND
	site.is_donor = 0;

INSERT INTO hbond_heavy_sites
SELECT
	NULL,
	site.struct_id,
	site.resNum,
	site.HBChemType,
	site.is_donor,
	res.name3 AS res_type,
	CASE WHEN SUM(env.num_hbonds) == 0 THEN 'don' ELSE 'DON' END AS satisfied_label,
	SUM(env.num_hbonds) != 0 AS satisfied,
	CASE WHEN SUM(env.sasa_r100) == 0 THEN 'Buried' ELSE 'Exposed' END AS buried_label,
	SUM(env.sasa_r100) == 0 AS buried
FROM
	residues AS res,
	hbond_sites AS site,
	hbond_site_environment AS env,
	hbond_site_atoms AS atms
WHERE
	site.struct_id = res.struct_id AND
	site.resNum = res.resNum AND
	env.struct_id = site.struct_id AND
	env.site_id = site.site_id AND
	site.is_donor = 1 AND
	atms.struct_id = site.struct_id AND
	atms.site_id = site.site_id
GROUP BY
	atms.base_x, atms.base_y, atms.base_z;

DROP INDEX IF EXISTS hbond_heavy_sites_struct_id_residue_number;
CREATE INDEX hbond_heavy_sites_struct_id_residue_number ON
	hbond_heavy_sites ( struct_id, residue_number );

SELECT count(*) FROM hbond_heavy_sites;"
#heavy_sites <- query_sample_sources(sample_sources, sele);

#sele <- "
#SELECT
#	don.res_type,
#	don.satisfied_label AS don_satisfied,
#	acc.satisfied_label AS acc_satisfied,
#	CASE WHEN don.buried AND acc.buried THEN 'Buried' ELSE 'Exposed' END AS buried,
#	don.HBChemType AS don_chem_type,
#	acc.HBChemType AS acc_chem_type,
#	count(*) AS count
#FROM
#	hbond_heavy_sites AS don,
#	hbond_heavy_sites AS acc
#WHERE
#	don.struct_id = acc.struct_id AND
#	don.residue_number = acc.residue_number AND
#	don.is_donor = 1 AND
#	acc.is_donor = 0
#GROUP BY
#	don.res_type, don_satisfied, acc_satisfied, buried, don_chem_type, acc_chem_type;"
#f <- query_sample_sources(sample_sources, sele)
#
#print(summary(f))
#f <- f[f$count > 100,]
#
plot_parts <- list(
	theme_bw(),
	geom_bar(aes(y=count, fill=sample_source), position="dodge", stat="identity"),
	scale_x_discrete("Residue Type"),
	scale_y_continuous("Residue Counts"))
#
#if(nrow(sample_sources) <= 3){
# plot_parts <- c(plot_parts,
#	list(opts(legend.position="bottom", legend.direction="horizontal")))
#}
#
#
#f_sc <- f[f$don_chem_type != "hbdon_PBA" & f$acc_chem_type != "hbacc_PBA",]
#f_sc_burial_agnostic <- ddply(f_sc,
#	.(res_type, don_satisfied, acc_satisfied, sample_source),
#	transform, count=sum(count))
#
#plot_id <- "sidechain_propensity_to_donate_or_accept_hbonds"
#p <- ggplot(data=f_sc_burial_agnostic, aes(x=don_satisfied)) + plot_parts +
#	opts(title="Sidechain Propensity to Donate or Accept Hydrogen Bonds (unsat/SAT)") +
#	facet_grid(acc_satisfied ~ res_type)
#save_plots(self, plot_id, sample_sources, output_dir, output_formats)
#
#
######################################
#plot_id <- "sidechain_propensity_to_donate_or_accept_hbonds_buried"
#p <- ggplot(data=f_sc[f_sc$buried=="Buried",], aes(x=don_satisfied)) + plot_parts +
#	opts(title="Buried Sidechain Propensity to Donate or Accept Hydrogen Bonds (unsat/SAT)") +
#	facet_grid(acc_satisfied ~ res_type)
#save_plots(self, plot_id, sample_sources, output_dir, output_formats)
#
####################################
#plot_id <- "sidechain_propensity_to_donate_or_accept_hbonds_exposed"
#p <- ggplot(data=f_sc[f_sc$buried=="Exposed",], aes(x=don_satisfied)) + plot_parts +
#	opts(title="Exposed Sidechain Propensity to Donate or Accept Hydrogen Bonds (unsat/SAT)") +
#	facet_grid(acc_satisfied ~ res_type)
#save_plots(self, plot_id, sample_sources, output_dir, output_formats)
#
##################################
##f_bbsc <- f[f$don_chem_type == "hbdon_PBA" & f$acc_chem_type == "hbacc_PBA",]
##f_bb$don_acc_satisfied <- paste(f$don_satisfied, f$acc_satisfied, sep=" ")
##
##plot_id <- "backbone_propensity_to_donate_or_accept_hbonds"
##p <- ggplot(data=f_bb, aes(x=don_acc_satisfied)) + theme_bw() +
##	opts(title="Backbone Propensity to Donate or Accept Hydrogen Bonds (SAT / unsat)") +
##	opts(axis.text.x=theme_text(angle=-90, hjust=0)) +
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
	hbond_heavy_sites AS site
WHERE
	site.HBChemType != 'hbdon_PBA' AND
	site.HBChemType != 'hbacc_PBA'
GROUP BY
	satisfied, buried, chem_type;"
f <- query_sample_sources(sample_sources, sele)
f <- f[f$count > 100,]

f$chem_type_name <- chem_type_name_linear(f$chem_type)

plot_id <- "hbond_site_propensity_to_hydrogen_bond"
p <- ggplot(data=f, aes(x=satisfied)) + plot_parts +
	opts(title="H-Bond Site Propensity to Hydrogen Bond") +
	opts(axis.text.x=theme_text(angle=-90, hjust=0)) +
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


#f <- merge(f, sample_sources[,c("sample_source", "reference")])
#
#sat_test <- ddply(f, .(chem_type, buried), function(df){
#	adply(data.sample_sources[sample_sources$reference == 1, "sample_source"], 1, function(ref_ss){
#		adply(sample_sources[sample_sources$reference == 0, "sample_source"], 1, function(new_ss){
#			ref_df <- df[df$sample_source == ref_ss$sample_source,]
#			new_df <- df[df$sample_source == new_ss$sample_source,]
#			z <- matrix(c(
#				ref_df[ref_df$satisfied == "Sat", "count"],
#				ref_df[ref_df$satisfied == "UnSat", "count"],
#				new_df[new_df$satisfied == "Sat", "count"],
#				new_df[new_df$satisfied == "UnSat", "count"]), nrow=2)
#			tryCatch({
#				t <- fisher.test(z)
#			}, error=function(e){
##				cat(paste("ERROR: Computing Fisher's Exact Test:\n", e, sep=""))
##				cat("ref:\n")
##				print(ref_df)
##				cat("new:\n")
##				print(new_df)
#			})
#			data.frame(new_sample_source = new_ss$sample_source[1], p.value=t$p.value, odds_ratio=t$estimate)
#		})
#	})
#})
#
#
#sat_test <- ddply(f[f$reference == 1,], .(chem_type, buried, sample_source), function(ref_df){
#	ddply(f[f$reference == 0,], .(chem_type, buried, sample_source), function(new_df){
#		z <- matrix(c(
#			ref_df[ref_df$satisfied == "Sat", "count"],
#			ref_df[ref_df$satisfied == "UnSat", "count"],
#			new_df[new_df$satisfied == "Sat", "count"],
#			new_df[new_df$satisfied == "UnSat", "count"]), nrow=2)
#		tryCatch({
#			t <- fisher.test(z)
#		}, error=function(e){
##			cat(paste("ERROR: Computing Fisher's Exact Test:\n", e, sep=""))
##			cat("ref:\n")
##			print(ref_df)
##			cat("new:\n")
##			print(new_df)
#		})
#		data.frame(new_sample_source = new_df$sample_source[1], p.value=t$estimate)
#	})
#})
#
#	ddply(f[f$reference == 0,], .(chem_type, buried, sample_source), nrow)
#	})
#
#
#
#p_sat <- ddply(f, .(chem_type, buried, sample_source), function(df){
#	z <- data.frame(
#		p_sat = sum(df[df$satisfied == "Sat", "count"])/sum(df$count))
#	z$var <- z$p_sat * (1-z$p_sat)
#	z
#})
#


#p_sat <- dcast(f, buried + chem_type_name ~ sample_source, fun.aggregate=c(
#	function(df){
#		sum(df[df$satisfied == "Sat", "count"])/sum(df$count)
#	},
#	function(df){
#		p <- sum(df[df$satisfied == "Sat", "count"])/sum(df$count)
#		p*(1-p)
#	}))
#
#frac_sat <- function(df) sum(df[df$satisfied == "Sat", "count"])/sum(df$count)
#frac_sat_var <- function(df){
#	p <- sum(df[df$satisfied == "Sat", "count"])/sum(df$count)
#	p*(1-p)
#}
#p_sat <- dcast(f, buried + chem_type_name ~ sample_source, fun.aggregate=c("frac_sat"))


})) # end FeaturesAnalysis
