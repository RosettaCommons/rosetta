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
id = "rotamer_recovery_by_res_type",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("ResidueFeatures", "RotamerRecoveryFeatures", "HBondFeatures", "ResidueSecondaryStructureFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
CREATE TEMPORARY TABLE max_residue_bfactors AS
SELECT
	hb_pdb_site.struct_id as struct_id,
  hb_pdb_site.resNum as resNum,
  MAX( hb_pdb_site.heavy_atom_temperature ) as max_temp
FROM
  hbond_sites_pdb as hb_pdb_site
GROUP BY
	hb_pdb_site.struct_id,
  hb_pdb_site.resNum;

SELECT
	rr.divergence AS divergence,
	res.name3 AS res_type
FROM
	rotamer_recovery AS rr,
	residues AS res,
	max_residue_bfactors AS b,
	residue_burial AS bur
WHERE
	res.resNum = rr.resNum AND
	res.struct_id = rr.struct_id AND
	b.struct_id = res.struct_id AND b.resNum = res.resNum AND
	b.max_temp < 20 AND
	bur.struct_id = res.struct_id AND bur.resNum = res.resNum AND
	bur.sasa_r140 = 0;"

#SELECT
#	rr.divergence AS divergence,
#	res.name3 AS res_type
#FROM
#	rotamer_recovery AS rr,
#	residues AS res,
#	residue_pdb_confidence AS res_confidence
#WHERE
#	res.resNum = rr.resNum AND
#	res.struct_id = rr.struct_id AND
#	res_confidence.struct_id = res.struct_id AND
#	res_confidence.residue_number = res.resNum AND
#	res_confidence.max_temperature < 20;"

f <- query_sample_sources(sample_sources, sele)

plot_id <- "rotamer_recovery_by_res_type"
p <- ggplot(data=f) + theme_bw() +
	geom_histogram(aes(x=divergence, y=log(..count..), fill=sample_source), breaks=c(0, 1, 2, 3, 4), position="dodge") +
#	geom_indicator(aes(indicator=counts, colour=sample_source)) +
#	geom_indicator(aes(indicator=mean, color=sample_source, xpos="left")) +
	facet_wrap( ~ res_type ) +
	ggtitle("Rotamer Recovery by Residue Type, 0 SASA, B-Factor < 20") +
	labs(x="Recovery Score", y="log(counts)")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)




exp <- ddply(f, .(sample_source, res_type),
	function(df) data.frame(mean=round(mean(df$divergence), 4)))

dens <- estimate_density_1d(
	data = f,
	ids = c("sample_source", "res_type"),
	variable = "divergence")

dens <- merge(dens, exp)

p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=log(x+1), y=log(y+1), colour=sample_source)) +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	geom_indicator(aes(indicator=mean, color=sample_source, group=sample_source), xpos="left") +
	facet_wrap( ~ res_type ) +
	ggtitle("Rotamer Recovery by Residue Type, B-Factor < 20") +
	labs(x="log(Automorphic RMSD + 1)", y="log(FeatureDensity + 1)")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

})) # end FeaturesAnalysis
