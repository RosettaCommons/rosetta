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
id = "AHdist_bbbb",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("StructureFeatures", "HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
	struct.tag,
	acc.site_id AS id,
	'' AS chain,
	acc.resNum,
	'CA' AS CA, 'C' AS C, 'O' AS O,
	CASE don.resNum - acc.resNum
		WHEN -1 THEN '-1' WHEN -2 THEN '-2' WHEN -3 THEN '-3' WHEN -4 THEN '-4'
		WHEN 1 THEN '1' WHEN 2 THEN '2' WHEN 3 THEN '3' WHEN 4 THEN '4'
		ELSE 'long' END AS seq_sep
FROM
	structures as struct,
	hbonds AS hbond,
	hbond_sites AS don, hbond_sites AS acc
WHERE
	hbond.struct_id = struct.struct_id AND
	don.struct_id = struct.struct_id AND don.site_id = hbond.don_id AND
	acc.struct_id = struct.struct_id AND acc.site_id = hbond.acc_id AND
	acc.HBChemType == 'hbacc_PBA' AND don.HBChemType == 'hbdon_PBA'
LIMIT 15;"
f <- query_sample_sources(sample_sources, sele)

#TODO fix bug:


#Error in `[.data.frame`(sort_df(data, variables), , c(variables, "value"),  : 
#  undefined columns selected
#Calls: <Anonymous> ... lapply -> is.vector -> condense -> [ -> [.data.frame
#Execution halted

#ss_ids <- as.character(unique(f$sample_source))
#g <- cast(f, tag + chain + resNum + CA + C + O ~ sample_source, value = "divergence")
#g$rel_div <- g[,ss_ids[1]] - g[,ss_ids[2]]
#g <- g[order(g$rel_div),]
#g$id <- 1:nrow(g)
#
#
#dens <- estimate_density_1d(
#	data = g, ids = c(), variable = "rel_div", n_pts=500)
#
#plot_id <- "rotamer_recovery_LYS_relative_divergence"
#p <- ggplot(data=dens) + theme_bw() +
#	geom_line(aes(x=x, y=log(y+1))) +
#	ggtitle(paste("Lysine Relative Rotamer Recovery\nbetween ", paste(ss_ids, collapse=" and "),", B-Factor < 20", sep="")) +
#	labs(x="Automorphic RMSD_1 - Automorphic RMSD_2", y="log(FeatureDensity + 1)")
#if(nrow(sample_sources) <= 3){
#	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
#}
#save_plots(self, plot_id, sample_sources[sample_sources$sample_source %in% ss_ids,], output_dir, output_formats)
#
#n_examples <- 50
#gm <- melt(g[g$id <= n_examples,],
#	id.vars=c("sample_source", "tag", "id", "chain", "resNum"),
#	measure.vars=c("CA", "C", "O"),
#	variable_name = "atom")
#
#instances_id <- "AHdist_bbbb"
#prepare_feature_instances(instances_id, sample_sources, gm, output_dir)

})) # end FeaturesAnalysis
