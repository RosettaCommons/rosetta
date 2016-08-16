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
id = "cosAHD_qq_chem_type",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

source("scripts/analysis/plots/hbonds/hbond_geo_dim_scales.R")

sele <-"
SELECT
	geom.cosAHD,
	don.HBChemType AS don_chem_type, acc.HBChemType AS acc_chem_type
FROM
	hbonds AS hb,
	hbond_geom_coords AS geom,
	hbond_sites AS don, hbond_sites AS acc,
	hbond_sites_pdb AS don_pdb, hbond_sites_pdb AS acc_pdb
WHERE
	geom.struct_id = hb.struct_id AND geom.hbond_id = hb.hbond_id AND
	don.struct_id = hb.struct_id AND don.site_id = hb.don_id AND
	acc.struct_id = hb.struct_id AND acc.site_id = hb.acc_id AND
	don_pdb.struct_id = hb.struct_id AND don_pdb.site_id = hb.don_id AND
	don_pdb.heavy_atom_temperature < 30 AND
	acc_pdb.struct_id = hb.struct_id AND acc_pdb.site_id = hb.acc_id AND
	acc_pdb.heavy_atom_temperature < 30 AND
	abs(don.resNum - acc.resNum) > 5;"

f <- query_sample_sources(sample_sources, sele)
f$AHD <- acos(f$cosAHD)
f <- na.omit(f, method="r")

compute_quantiles <- function(
	data,
	ids,
	variable,
	n_quantiles=300
) {
	ddply(data, ids, function(df){
		ps <- ppoints(n_quantiles)
		data.frame(
			probs=ps, quantiles=quantile(df[,variable], probs=ps))
	})
}


ss_id.ref <- f$sample_source[1]
qs <- compute_quantiles(
	f, c("sample_source", "don_chem_type", "acc_chem_type"), "AHD")

qs.ref <- qs[qs$sample_source == ss_id.ref,]
names(qs.ref) <- c(
	"sample_source.ref", "don_chem_type", "acc_chem_type", "probs", "quantiles.ref")
qs.new <- qs[qs$sample_source != ss_id.ref,]
names(qs.new) <- c(
	"sample_source.new", "don_chem_type", "acc_chem_type", "probs", "quantiles.new")
qs.merged <- merge(qs.ref, qs.new)
qs.merged$don_chem_type_name <- don_chem_type_name_linear(qs.merged$don_chem_type)
qs.merged$acc_chem_type_name <- acc_chem_type_name_linear(qs.merged$acc_chem_type)

plot_id = "hbond_AHD_qq_chem_type"
p <- ggplot(data=qs.merged) + theme_bw() +
	geom_abline() +
	geom_line(aes(x=quantiles.ref, y=quantiles.new, colour=sample_source.new)) +
	facet_grid(don_chem_type_name ~ acc_chem_type_name) +
	ggtitle(paste("HBond AHD Angle Quantiles Against ", ss_id.ref, " Quantiles SeqSep > 5, B-Fact < 30", sep="")) +
	scale_x_continuous(paste(ss_id.ref, " Quantiles", sep="")) +
	scale_y_continuous("New Sample Source Quantiles")

if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
