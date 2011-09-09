# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

check_setup()

plot_id = "hbond_AHD_chi_polar_scatter"

sele <-"
SELECT
  geom.cosBAH,
  geom.chi,
  acc_site.HBChemType AS acc_chem_type,
  don_site.HBChemType AS don_chem_type
FROM
  hbond_geom_coords AS geom,
  hbonds AS hbond,
  hbond_sites AS don_site,
  hbond_sites AS acc_site
WHERE
  hbond.struct_id = geom.struct_id AND hbond.hbond_id = geom.hbond_id AND
  hbond.struct_id = don_site.struct_id AND hbond.don_id = don_site.site_id AND
  hbond.struct_id = acc_site.struct_id AND hbond.acc_id = acc_site.site_id;";

f <- query_sample_sources(sample_sources, sele)

sub_f <- ddply(f, .variables=c("don_chem_type", "acc_chem_type"),
  function(df){sample_rows(df, 5000)})

p <- ggplot(data=some_geom, aes(x=chi, y=sin(acos(cosBAH)), colour=sample_source))
p <- p + geom_contour( size=.5)
p <- p + coord_polar("x")
p <- p + facet_grid(acc_chem_type ~  don_chem_type)
p <- p + opts(title = "Hydrogen Bonds chi vs sinBAH Angles by Chemical Type\n(normalized for equal volume per unit distance)")
p <- p + labs(x=expression(paste('Acceptor -- Hydrogen --   Donor (degrees)')),
              y="log(FeatureDensity + 1)")
p <- p + theme_bw()

save_plots(plot_id, sample_sources, output_dir, output_formats)
