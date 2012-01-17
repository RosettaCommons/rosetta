# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

check_setup()
feature_analyses <- c(feature_analyses, new("FeatureAnalysis",
id = "OHdonor",
filename = "scripts/analysis/plots/hbonds/OHdonor.R",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(){

sele <-"
SELECT
  geom.AHdist,
  geom.cosBAH,
  geom.cosAHD,
  geom.chi,
  acc_site.HBChemType AS acc_chem_type,
  don_site.HBChemType AS don_chem_type
FROM
  hbond_geom_coords AS geom,
  hbonds AS hbond,
  hbond_sites AS don_site,
  hbond_sites AS acc_site
WHERE
  hbond.struct_id = geom.struct_id AND
  hbond.hbond_id = geom.hbond_id AND
  hbond.struct_id = don_site.struct_id AND
  hbond.don_id = don_site.site_id AND
  hbond.struct_id = acc_site.struct_id AND
  hbond.acc_id = acc_site.site_id AND
  ( don_site.HBChemType = 'hbdon_AHX' OR don_site.HBChemType = 'hbdon_HXL' );";

all_geom <- query_sample_sources(sample_sources, sele)

plot_id <- "OHdonor_AHdist_all_acceptor_types"
dens <- estimate_density_1d(
  data = all_geom,
  ids = c("sample_source"),
  variable = "AHdist",
  weight_fun = radial_3d_normalization)
p <- ggplot(data=dens, aes(x=x, y=log(y+1), colour=sample_source, indicator=counts))
p <- p + geom_line()
p <- p + geom_indicator()
p <- p + opts(title = "Hydroxyl Donor Hydrogen Bonds A-H Distance\n(normalized for equal weight per unit distance)")
p <- p + labs(x=expression(paste('Acceptor -- Proton Distance (', ring(A), ')')),
              y="log(FeatureDensity + 1)")
p <- p + theme_bw()
p <- p + opts(axis.text.y=theme_blank())
p <- p + scale_y_continuous(limits=c(0,2.9), breaks=0:2)
p <- p + scale_x_continuous(limits=c(1.4,2.7), breaks=c(1.6, 1.9, 2.2, 2.6))
save_plots(plot_id, sample_sources, output_dir, output_formats)

plot_id <- "OHdonor_cosBAH_all_acceptor_types"
dens <- estimate_density_1d(
  data = all_geom,
  ids = c("sample_source"),
  variable = "cosBAH")
p <- ggplot(data=dens, aes(x=acos(x)*180/pi, y=log(y+1), colour=sample_source, indicator=counts))
p <- p + geom_line()
p <- p + geom_indicator()
p <- p + opts(title = "Hydroxyl Donor Hydrogen Bonds BAH Angle \n(normalized for equal weight per unit distance)")
p <- p + labs(x=expression(paste('Base -- Acceptor -- Hydrogen (degrees)')),
              y="log(FeatureDensity + 1)")
p <- p + theme_bw()
p <- p + opts(axis.text.y=theme_blank())
save_plots(plot_id, sample_sources, output_dir, output_formats)

plot_id <- "OHdonor_cosAHD_all_acceptor_types"
dens <- estimate_density_1d(
  data = all_geom,
  ids = c("sample_source"),
  variable = "cosAHD")
p <- ggplot(data=dens, aes(x=acos(x)*180/pi, y=log(y+1), colour=sample_source, indicator=counts))
p <- p + geom_line()
p <- p + geom_indicator()
p <- p + opts(title = "Hydroxyl Donor Hydrogen Bonds BAH Angle \n(normalized for equal weight per unit distance)")
p <- p + labs(x=expression(paste('Acceptor -- Hydrogen -- Donor (degrees)')),
              y="log(FeatureDensity + 1)")
p <- p + theme_bw()
p <- p + opts(axis.text.y=theme_blank())
save_plots(plot_id, sample_sources, output_dir, output_formats)

plot_id <- "OHdonor_chi_all_acceptor_types"
dens <- estimate_density_1d_wrap(
  data = all_geom,
  ids = c("sample_source"),
  variable = "chi")
p <- ggplot(data=dens, aes(x=acos(x)*360/pi, y=log(y+1), colour=sample_source, indicator=counts))
p <- p + geom_line()
p <- p + geom_indicator()
p <- p + opts(title = "Hydroxyl Donor Hydrogen Bonds BAH Angle \n(normalized for equal weight per unit distance)")
p <- p + labs(x=expression(paste('Acceptor Base -- Acceptor Torsion (degrees)')),
              y="log(FeatureDensity + 1)")
p <- p + theme_bw()
p <- p + opts(axis.text.y=theme_blank())
save_plots(plot_id, sample_sources, output_dir, output_formats)



sidechain_geom <- all_geom[all_geom$acc_chem_type != "hbacc_PBA",]

plot_id <- "OHdonor_AHdist_sidechain_acceptor_types"
dens <- estimate_density_1d(
  data = sidechain_geom,
  ids = c("sample_source"),
  variable = "AHdist",
  weight_fun = radial_3d_normalization)
p <- ggplot(data=dens, aes(x=x, y=log(y+1), colour=sample_source, indicator=counts))
p <- p + geom_line()
p <- p + geom_indicator()
p <- p + opts(title = "Hydroxyl Donor Hydrogen Bonds A-H Distance to Sidechain Acceptors\n(normalized for equal weight per unit distance)")
p <- p + labs(x=expression(paste('Acceptor -- Proton Distance (', ring(A), ')')),
              y="log(FeatureDensity + 1)")
p <- p + theme_bw()
p <- p + opts(axis.text.y=theme_blank())
p <- p + scale_y_continuous(limits=c(0,2.9), breaks=0:2)
p <- p + scale_x_continuous(limits=c(1.4,2.7), breaks=c(1.6, 1.9, 2.2, 2.6))
save_plots(plot_id, sample_sources, output_dir, output_formats)

plot_id <- "OHdonor_cosBAH_sidechain_acceptor_types"
dens <- estimate_density_1d(
  data = sidechain_geom,
  ids = c("sample_source"),
  variable = "cosBAH")
p <- ggplot(data=dens, aes(x=acos(x)*180/pi, y=log(y+1), colour=sample_source, indicator=counts))
p <- p + geom_line()
p <- p + geom_indicator()
p <- p + opts(title = "Hydroxyl Donor Hydrogen Bonds BAH Angle to Sidechain Acceptors\n(normalized for equal weight per unit distance)")
p <- p + labs(x=expression(paste('Base -- Acceptor -- Hydrogen (degrees)')),
              y="log(FeatureDensity + 1)")
p <- p + theme_bw()
p <- p + opts(axis.text.y=theme_blank())
save_plots(plot_id, sample_sources, output_dir, output_formats)

plot_id <- "OHdonor_cosAHD_sidechain_acceptor_types"
dens <- estimate_density_1d(
  data = sidechain_geom,
  ids = c("sample_source"),
  variable = "cosAHD")
p <- ggplot(data=dens, aes(x=acos(x)*180/pi, y=log(y+1), colour=sample_source, indicator=counts))
p <- p + geom_line()
p <- p + geom_indicator()
p <- p + opts(title = "Hydroxyl Donor Hydrogen Bonds AHD Angle to Sidechain Acceptors\n(normalized for equal weight per unit distance)")
p <- p + labs(x=expression(paste('Acceptor -- Hydrogen -- Donor (degrees)')),
              y="log(FeatureDensity + 1)")
p <- p + theme_bw()
p <- p + opts(axis.text.y=theme_blank())
save_plots(plot_id, sample_sources, output_dir, output_formats)

plot_id <- "OHdonor_chi_sidechain_acceptor_types"
dens <- estimate_density_1d_logspline(
  data = sidechain_geom,
  ids = c("sample_source"),
  variable = "chi")
p <- ggplot(data=dens, aes(x=acos(x)*360/pi, y=log(y+1), colour=sample_source, indicator=counts))
p <- p + geom_line()
p <- p + geom_indicator()
p <- p + opts(title = "Hydroxyl Donor Hydrogen Bonds chi Torsion Angle to Sidechain Acceptors\n(normalized for equal weight per unit distance)")
p <- p + labs(x=expression(paste('Acceptor Base -- Acceptor Torsion (degrees)')),
              y="log(FeatureDensity + 1)")
p <- p + theme_bw()
p <- p + opts(axis.text.y=theme_blank())
save_plots(plot_id, sample_sources, output_dir, output_formats)




})) # end FeatureAnalysis
