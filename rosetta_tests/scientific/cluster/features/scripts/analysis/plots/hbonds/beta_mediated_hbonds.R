# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

check_setup()

sele <-"
SELECT
  structure.tag,
  don_site.chain,
  don_site.resNum,
  don_site.atmNum,
  acc_site.chain,
  acc_site.resNum,
  acc_site.atmNum,
  hbond.energy,
  hbond.HBEvalType,
  geom.AHdist,
  geom.cosBAH,
  geom.cosAHD,
  geom.chi,
  don_env.dssp AS don_dssp,
  acc_env.dssp AS acc_dssp
FROM
  hbond_geom_coords AS geom,
  hbonds AS hbond,
  structures AS structure,
  hbond_sites AS don_site,
  hbond_sites AS acc_site,
  hbond_site_environment AS don_env,
  hbond_site_environment AS acc_env
WHERE
  hbond.struct_id  = structure.struct_id AND
  hbond.struct_id  = geom.struct_id AND
  hbond.hbond_id   = geom.hbond_id AND
  hbond.struct_id  = don_site.struct_id AND
  hbond.don_id     = don_site.site_id AND
  hbond.struct_id  = acc_site.struct_id AND
  hbond.acc_id     = acc_site.site_id AND
  hbond.struct_id  = don_env.struct_id AND
  hbond.don_id     = don_env.site_id AND
  hbond.struct_id  = acc_env.struct_id AND
  hbond.acc_id     = acc_env.site_id;"
#--  don_env.dssp     = 'E' AND -- beta sheet
#--  acc_env.dssp     = 'E' AND -- beta sheet
#--  don_site.chain  != acc_site.chain AND -- to make the hbonds be across an interface
#--  don_site.HBChemType = 'hbdon_PBA' AND -- protein backbone amide
#--  acc_site.HBChemType = 'hbacc_PBA';"

f <- query_sample_sources(sample_sources, sele)

dens <- estimate_density_1d(
  data = f,
  ids = c("sample_source"),
  variable = "AHdist",
  weight_fun = radial_3d_normalization)

plot_id <- "beta_mediated_hbonds_AHdist"
ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x, y=log(y+1), colour=sample_source)) +
	geom_indicator(aes(indicator=counts, colour=sample_source)) +
	opts(title = "Hydrogen Bonds A-H Distance For beta-mediated HBonds\nnormalized for equal weight per unit distance") +
	labs(x=expression(paste('Acceptor -- Proton Distance (', ring(A), ')')),
	     y="log(FeatureDensity + 1)") +
	scale_y_continuous(limits=c(0,2.9), breaks=0:2) +
	scale_x_continuous(limits=c(1.4,2.7), breaks=c(1.6, 1.9, 2.2, 2.6)) +
save_plots(plot_id, sample_sources, output_dir, output_formats)
