# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

# Commend out plot because the secondary structure feature APSA is not integrated into the FeaturesReporter framework.

#check_setup()
#
#sele <- "
#CREATE TEMPORARY TABLE IF NOT EXISTS hbond_cosAHD AS SELECT
#  hbond.struct_id AS struct_id,
#  don_site.resNum AS don_resNum,
#  acc_site.resNum AS acc_resNum,
#  geom.cosAHD AS cosAHD
#FROM
#  hbonds AS hbond,
#  hbond_geom_coords AS geom,
#  hbond_sites AS don_site,
#  hbond_sites AS acc_site
#WHERE
#  hbond.struct_id = geom.struct_id AND
#  hbond.hbond_id =  geom.hbond_id AND
#  hbond.struct_id = don_site.struct_id AND
#  hbond.don_id = don_site.site_id AND
#  hbond.struct_id = acc_site.struct_id AND
#  hbond.acc_id = acc_site.site_id AND
#  don_site.HBChemType = 'hbdon_PBA' AND
#  acc_site.HBChemType = 'hbacc_PBA';
#
#CREATE INDEX hbond_cosAHD_don_site on
#  hbond_cosAHD ( struct_id, don_resNum );
#
#CREATE INDEX hbond_cosAHD_acc_site on
#  hbond_cosAHD ( struct_id, acc_resNum );
#
#ANALYZE;
#
#SELECT
#  don_apsa.secondary_struct AS don_ss,
#  acc_apsa.secondary_struct AS acc_ss,
#  hb.cosAHD AS cosAHD
#FROM
#  hbond_cosAHD AS hb,
#  apsa AS don_apsa,
#  apsa AS acc_apsa
#WHERE
#  don_apsa.struct_id = hb.struct_id AND
#  acc_apsa.struct_id = don_apsa.struct_id AND
#  don_apsa.resNum = hb.don_resNum  AND
#  acc_apsa.resNum = hb.acc_resNum;"
#
#
#all <-  query_sample_sources(
#  sample_sources,
#  sele)
#
#all$don_ss <- factor(all$don_ss)
#all$acc_ss <- factor(all$acc_ss)
#
#dens <- estimate_density_1d(
#  data = all,
#  ids = c("sample_source", "acc_ss", "don_ss"),
#  variable = "cosAHD")
#
#plot_id <- "backbone_backbone_hbond_cosAHD_trellis"
#p <- ggplot(data=dens, aes(x=acos(x)*360/pi, y=log(y+1), colour=sample_source, indicator=counts))
#p <- p + geom_line()
#p <- p + geom_indicator()
#p <- p + facet_grid( don_ss ~ acc_ss)
#p <- p + ggtitle("Backbone Backbone Hydrogen Bond AHD Angle by Donor Secondary Structure (APSA Method)\nnormalized for equal weight per unit distance")
#p <- p + labs(x=expression(paste('Acceptor -- Hydrogen -- Donor (Degrees)')),
#              y="log(FeatureDensity + 1)")
#p <- p + theme_bw()
#p <- p + theme(axis.text.y=theme_blank())
#
#save_plots(self, plot_id, sample_sources, output_dir, output_formats)
#
#
#dens <- estimate_density_1d(
#  data = all,
#  ids = c("sample_source", "acc_ss"),
#  variable = "cosAHD")
#
#plot_id <- "backbone_backbone_hbond_cosAHD_diversity_acceptor"
#p <- ggplot(data=dens, aes(x=acos(x)*360/pi, y=log(y+1), colour=acc_ss))
#p <- p + geom_line()
#p <- p + facet_wrap(  ~ sample_source )
#p <- p + ggtitle("Backbone Backbone Hydrogen Bond AHD Angle by Donor Secondary Structure (APSA Method)\nnormalized for equal weight per unit distance")
#p <- p + labs(x=expression(paste('Acceptor -- Hydrogen -- Donor (Degrees)')),
#              y="log(FeatureDensity + 1)")
#p <- p + theme_bw()
#p <- p + theme(axis.text.y=theme_blank())
#
#save_plots(self, plot_id, sample_sources, output_dir, output_formats)
#
#dens <- estimate_density_1d(
#  data = all,
#  ids = c("sample_source", "don_ss"),
#  variable = "cosAHD")
#
#
#plot_id <- "backbone_backbone_hbond_cosAHD_diversity_donor"
#p <- ggplot(data=dens, aes(x=acos(x)*360/pi, y=log(y+1), colour=don_ss))
#p <- p + geom_line()
#p <- p + facet_wrap(  ~ sample_source )
#p <- p + ggtitle("Backbone Backbone Hydrogen Bond AHD Angle by Donor Secondary Structure (APSA Method)\nnormalized for equal weight per unit distance")
#p <- p + labs(x=expression(paste('Acceptor -- Hydrogen -- Donor (Degrees)')),
#              y="log(FeatureDensity + 1)")
#p <- p + theme_bw()
#p <- p + theme(axis.text.y=theme_blank())
#
#save_plots(self, plot_id, sample_sources, output_dir, output_formats)
