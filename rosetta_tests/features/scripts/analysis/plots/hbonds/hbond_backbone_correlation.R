# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

check_setup()

# Organization of joins to form query
#rsd
#	acc_site
#		acc_hb
#			acc_geom
#		acc_env
#	don_site
#		don_hb
#			don_geom
#		don_env
sele <- "
SELECT
	acc_geom.AHdist AS acc_AHdist,
	don_geom.AHdist AS don_AHdist,
	acc_hb.energy AS acc_energy,
	don_hb.energy AS don_energy,
	acc_env.dssp AS dssp,
FROM
	residues AS rsd,
	hbond_sites AS acc_site,
	hbond_sites AS don_site,
	hbonds AS acc_hb,
	hbonds AS don_hb,
	hbond_geom_coords AS acc_geom,
	hbond_geom_coords AS don_geom,
	hbond_site_environment AS acc_env
WHERE
	acc_site.struct_id  = rsd.struct_id      AND acc_site.resNum   = rsd.resNum       AND
	acc_site.HBChemType = 'hbacc_PBA' AND
	acc_hb.struct_id    = acc_site.struct_id AND acc_hb.acc_id     = acc_site.site_id AND
	acc_geom.struct_id  = acc_hb.struct_id   AND acc_geom.hbond_id = acc_hb.hbond_id  AND
	acc_env.struct_id   = acc_site.struct_id AND acc_env.site_id   = acc_site.site_id AND
	don_site.struct_id  = rsd.struct_id      AND don_site.resNum   = rsd.resNum       AND
	don_site.HBChemType = 'hbdon_PBA' AND
	don_hb.struct_id    = don_site.struct_id AND don_hb.don_id     = don_site.site_id AND
	don_geom.struct_id  = don_hb.struct_id   AND don_geom.hbond_id = don_hb.hbond_id;"

all_geom <- query_sample_sources(sample_sources, sele)

plot_id <- "hbond_backbone_correlation"
p <- ggplot(data=all_geom) + theme_bw() +
	geom_point(aes(acc_AHdist, don_AHdist), size=.3) +
	facet_wrap(don_dssp) +
	opts(title = "HBond Backbone Amide Correlation by DSSP of residue") +
	labs(x=expression(paste('Backbone is Acceptor: Acceptor -- Proton Distance (', ring(A), ')')),
		y=expression(paste('Backbone is Donor: Acceptor -- Proton Distance (', ring(A), ')')))

save_plots(plot_id, sample_sources, output_dir, output_formats)
