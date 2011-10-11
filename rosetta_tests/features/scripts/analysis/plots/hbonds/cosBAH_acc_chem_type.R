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
  geom.cosBAH,
  acc_site.HBChemType AS acc_chem_type,
  don_site.HBChemType AS don_chem_type
FROM
  hbond_geom_coords AS geom,
  hbonds AS hbond,
  hbond_sites AS don_site,
  hbond_sites AS acc_site
WHERE
  hbond.struct_id = geom.struct_id AND
  hbond.hbond_id =  geom.hbond_id AND
  hbond.struct_id = don_site.struct_id AND
  hbond.don_id = don_site.site_id AND
  hbond.struct_id = acc_site.struct_id AND
  hbond.acc_id = acc_site.site_id;"

f <- query_sample_sources(sample_sources, sele)

f$acc_chem_type <- factor(f$acc_chem_type,
	levels =
		c("hbacc_HXL", "hbacc_CXL", "hbacc_IMD", "hbacc_PBA",
			"hbacc_AHX", "hbacc_CXA", "hbacc_IME"),
	labels =
		c("aHXL: s,t", "aCXL: d,e", "aIMD: h", "aPBA: bb",
			"aAHX: y",   "aCXA: n,q", "aIME: h"))

dens <- estimate_density_1d(
  f, c("sample_source", "acc_chem_type" ), "cosBAH")

plot_id = "cosBAH_acc_chem_type"
ggplot(data=dens) +
	geom_line(aes(x=acos(x)*180/pi, y=-log(y), colour=sample_source)) +
	geom_indicator(aes(colour=sample_source, indicator=counts)) +
	facet_grid( ~ acc_chem_type) +
	opts(title = "Hydrogen Bonds BAH Angle by Chemical Type\n(normalized for equal volume per unit distance)") +
	scale_x_continuous(paste('Base -- Acceptor -- Hydrogen (degrees)')) +
	scale_y_continuous("-log(FeatureDensity)", limits=c(-2.3,6)) +
	opts(legend.position=c(.58,.35)) +
	opts(legend.justification=c("left", "top"))
save_plots(plot_id, sample_sources, output_dir, output_formats)
