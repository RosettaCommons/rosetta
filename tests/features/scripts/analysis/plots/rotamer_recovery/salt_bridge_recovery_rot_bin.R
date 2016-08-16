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
id = "salt_bridge_recovery",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("StructureFeatures", "RotamerRecoveryFeatures", "HBondFeatures"),
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

CREATE TEMPORARY TABLE sc_hbond_card AS
SELECT
	don_site.struct_id AS struct_id,
	don_site.resNum AS don_resNum,
	acc_site.resNum AS acc_resNum,
	COUNT( hbond.hbond_id ) AS num_sc_hbs
FROM
	hbonds as hbond,
	hbond_sites as don_site,
	hbond_sites as acc_site
WHERE
	hbond.struct_id = don_site.struct_id and hbond.don_id = don_site.site_id and
	hbond.struct_id = acc_site.struct_id and hbond.acc_id = acc_site.site_id and
	don_site.HBChemType != 'hbdon_PBA' AND acc_site.HBChemType != 'hbacc_PBA'
GROUP BY
	hbond.struct_id,
	don_site.resNum,
	acc_site.resNum;

CREATE TEMPORARY TABLE arg_cxl_hbonds AS
SELECT
	hbond.hbond_id AS hbond_id,
	hbond.struct_id AS struct_id,
	hbond.don_id AS don_id,
	hbond.acc_id AS acc_id
FROM
	hbonds AS hbond,
	hbond_sites AS don_site,
	hbond_sites AS acc_site
WHERE
	hbond.struct_id = don_site.struct_id AND hbond.don_id = don_site.site_id AND
	hbond.struct_id = acc_site.struct_id AND hbond.acc_id = acc_site.site_id AND
	(don_site.HBChemType = 'hbdon_GDE' OR don_site.HBChemType = 'hbdon_GDH' ) AND
	acc_site.HBChemType = 'hbacc_CXL';

CREATE TEMPORARY TABLE arg_cxl_hbond_temps AS
SELECT
	hbond.hbond_id AS hbond_id,
	hbond.struct_id AS struct_id,
	hbond.don_id AS don_id,
	hbond.acc_id AS acc_id,
	don_site_pdb.resNum AS don_resNum,
	acc_site_pdb.resNum AS acc_resNum,
	don_max_temp.max_temp AS don_temp,
	acc_max_temp.max_temp AS acc_temp
FROM
	arg_cxl_hbonds AS hbond,
	hbond_sites_pdb AS don_site_pdb,
	hbond_sites_pdb AS acc_site_pdb,
	max_residue_bfactors AS don_max_temp,
	max_residue_bfactors AS acc_max_temp
WHERE
	hbond.struct_id = don_site_pdb.struct_id AND hbond.don_id = don_site_pdb.site_id AND
	hbond.struct_id = acc_site_pdb.struct_id AND hbond.acc_id = acc_site_pdb.site_id AND
	don_max_temp.struct_id = hbond.struct_id AND don_max_temp.resNum = don_site_pdb.resNum AND
	acc_max_temp.struct_id = hbond.struct_id AND acc_max_temp.resNum = acc_site_pdb.resNum;

SELECT
	structure.tag,
	don_site_pdb.chain, don_site_pdb.resNum, don_site_pdb.iCode,
	acc_site_pdb.chain, acc_site_pdb.resNum, acc_site_pdb.iCode,
	geom.AHdist,
	n_sc_hbonds.num_sc_hbs,
	don_res_recovery.divergence AS don_res_recovery,
	acc_res_recovery.divergence AS acc_res_recovery
FROM
	arg_cxl_hbond_temps AS hbond,
	hbond_geom_coords AS geom,
	sc_hbond_card as n_sc_hbonds,
	hbond_sites_pdb AS don_site_pdb,
	hbond_sites_pdb AS acc_site_pdb,
	structures AS structure,
	rotamer_recovery as don_res_recovery,
	rotamer_recovery as acc_res_recovery
WHERE
	hbond.don_temp < 20 AND hbond.acc_temp < 20 AND
	n_sc_hbonds.struct_id = hbond.struct_id AND
	n_sc_hbonds.don_resNum = hbond.don_resNum AND n_sc_hbonds.acc_resNum = hbond.acc_resNum AND
	geom.struct_id = hbond.struct_id AND geom.hbond_id = hbond.hbond_id AND
	don_site_pdb.struct_id = hbond.struct_id AND don_site_pdb.site_id = hbond.don_id AND
	acc_site_pdb.struct_id = hbond.struct_id AND acc_site_pdb.site_id = hbond.acc_id AND
	structure.struct_id = hbond.struct_id AND
	don_res_recovery.struct_id = hbond.struct_id AND don_res_recovery.resNum = don_site_pdb.resNum AND
	acc_res_recovery.struct_id = hbond.struct_id AND acc_res_recovery.resNum = acc_site_pdb.resNum;"

f <- query_sample_sources(sample_sources, sele)
f$don_res_recovery <- factor(f$don_res_recovery)


f  <- ddply(f, .(sample_source, num_sc_hbs), function(df){
  df$mean <- round(mean(as.numeric(df$don_res_recovery)), 4)
  df$counts <- nrow(df)
  df
})

plot_id <- "salt_bridge_ARG_with_ASP-GLU_recovery_rot_bin"
p <- ggplot(data=f) + theme_bw() +
	stat_bin(aes(x=don_res_recovery, colour=sample_source), geom="line") +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	geom_indicator(aes(indicator=mean, color=sample_source, group=sample_source), xpos="left") +
	facet_wrap( ~ num_sc_hbs, nrow=1 ) +
	ggtitle("ARG Rotamer Recovery when forming HBonds with ASP/GLU and ARG and Bfactors < 20\nBy Number of sc hbonds between the acceptor and donor residues") +
	labs(x="Chi Angles Not Recovered", y="FeatureDensity")

if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)




plot_id <- "salt_bridge_ASP-GLU_with_ARG_recovery"
exp <- ddply(f, .(sample_source, num_sc_hbs),
	function(df) data.frame(mean=round(mean(df$acc_res_recovery), 4)))
dens <- estimate_density_1d(
	f, c("sample_source", "num_sc_hbs" ),
	"acc_res_recovery", weight_fun = radial_3d_normalization)
dens <- merge(dens, exp)

p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x, y=y, colour=sample_source)) +
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	geom_indicator(aes(indicator=mean, color=sample_source, group=sample_source), xpos="left") +
	facet_wrap( ~ num_sc_hbs, nrow=1 ) +
	ggtitle("ASP/GLU Rotamer Recovery when forming HBonds with ARG and Bfactors < 20\nBy Number of sc hbonds between the acceptor and donor residues") +
  	labs(x="log(Automorphic RMSD + 1)", y="log(FeatureDensity + 1)")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)

})) # end FeaturesAnalysis
