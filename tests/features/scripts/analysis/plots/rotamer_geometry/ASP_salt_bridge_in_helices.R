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
id = "ASP_salt_bridge_in_helices",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("ResidueFeatures", "ProteinResidueConformationFeatures", "PdbDataFeatures", "ResidueSecondaryStructureFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


sele <-"
SELECT
	res_dofs.chi1 AS chi1,
	res_dofs.chi2 AS chi2
FROM
	residues AS res,
	residue_pdb_confidence AS res_conf,
	residue_secondary_structure AS res_ss,
	protein_residue_conformation AS res_dofs
WHERE
	res_conf.struct_id = res.struct_id AND res_conf.residue_number = res.resNum AND
	res_conf.max_temperature < 30 AND
	res_ss.struct_id = res.struct_id AND res_ss.resNum == res.resNum AND
	res_ss.dssp == 'H' AND
	res_dofs.struct_id = res.struct_id AND res_dofs.seqpos == res.resNum AND
	res.name3 == 'ASP';"

f <- query_sample_sources(sample_sources, sele)

f <- ddply(f, .(sample_source), transform, counts=length(sample_source))

dens <- ddply(f, .(sample_source), function(sub_f) {
	dens <- MASS::kde2d(sub_f$chi1, sub_f$chi2, n=500, h=1)
	densdf <- data.frame(expand.grid(x=dens$x, y=dens$y), z=as.vector(dens$z)/max(dens$z))
	densdf$counts <- nrow(sub_f)
	densdf
})

plot_id <-"ASP_in_helices"
p <- ggplot(data=dens) +
	theme_bw() +
	geom_raster(aes(x=x, y=y, fill=z)) +
#          stat_density2d(
#		aes(x=chi1, y=chi2, fill=..density..),
#		geom="tile", contour=F, geom_params=list(),
#                stat_params=list(contour=F)) +
	geom_indicator(
		aes(indicator=counts), group=1, colour="black", ypos=.99) +
	facet_wrap(~sample_source) +
	ggtitle(("ASP chi1 vs chi2 For Helices;  BFact < 30", sep="")) +
	coord_equal(ratio=1) +
	scale_x_continuous("chi 1") +
	scale_y_continuous("chi 2") +
	scale_fill_gradient("Density", low="white", high="black") +
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
