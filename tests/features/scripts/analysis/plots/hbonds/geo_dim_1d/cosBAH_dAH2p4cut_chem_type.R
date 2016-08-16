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
id = "cosBAH_dAH2p4cut_chem_type",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

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
  geom.AHdist <= 2.4 AND
  hbond.struct_id = geom.struct_id AND
  hbond.hbond_id =  geom.hbond_id AND
  hbond.struct_id = don_site.struct_id AND
  hbond.don_id = don_site.site_id AND
  hbond.struct_id = acc_site.struct_id AND
  hbond.acc_id = acc_site.site_id;"

all_geom <- query_sample_sources(sample_sources, sele)

dens <- estimate_density_1d(
  data = all_geom,
  ids = c("sample_source", "acc_chem_type", "don_chem_type"),
  variable = "cosBAH")

# The 'hbdon_' part of the donor labels doesn't fit so strip them out
dens$don_chem_type <- sub("^hbdon_", '', dens$don_chem_type)

plot_id = "cosBAH_dAH2p4cut_chem_type"
p <- ggplot(data=dens, aes(x=acos(x)*180/pi, y=y, colour=sample_source, indicator=counts))
p <- p + geom_line()
p <- p + geom_indicator()
p <- p + facet_grid( don_chem_type ~ acc_chem_type )
p <- p + ggtitle("Hydrogen Bonds BAH Angle by Chemical Type w/ dist(A,H) < 2.4\n(normalized for equal volume per unit distance)")
p <- p + labs(x=expression(paste('Base -- Acceptor -- Hydrogen (degrees)')),
              y="FeatureDensity")
p <- p + theme_bw()
p <- p + scale_y_continuous(limits=c(-2.3,6))

save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
