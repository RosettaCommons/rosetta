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
id = "chi_fit_cos_cosBAH_long_range",
author = "Matthew O'Meara",
brief_description = "",
long_description = "
These are example of asp/glu acceptors with his donors with long sequence
separation and the chi angle between 160 and 200 -> which is the syn
orbital sorted by hbond energy.",
feature_reporter_dependencies = c("StructureFeatures", "HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


examples_sele <-"
SELECT
  structure.tag,
  don_site_pdb.chain, don_site_pdb.resNum, don_site_pdb.iCode,
  acc_site_pdb.chain, acc_site_pdb.resNum, acc_site_pdb.iCode,
  geom.chi, geom.cosBAH,
  acc_site.HBChemType AS acc_chem_type,
  don_site.HBChemType AS don_chem_type,
  hbond.energy
FROM
  hbond_geom_coords AS geom,
  hbonds AS hbond,
  hbond_sites AS don_site,
  hbond_sites AS acc_site,
  hbond_sites_pdb AS don_site_pdb,
  hbond_sites_pdb AS acc_site_pdb,
  structures AS structure
WHERE
  hbond.struct_id = geom.struct_id AND hbond.hbond_id =  geom.hbond_id AND
  hbond.struct_id = don_site.struct_id AND hbond.don_id = don_site.site_id AND
  hbond.struct_id = acc_site.struct_id AND hbond.acc_id = acc_site.site_id AND
  hbond.struct_id = don_site_pdb.struct_id AND don_site_pdb.site_id = don_site.site_id AND
  hbond.struct_id = acc_site_pdb.struct_id AND acc_site_pdb.site_id = acc_site.site_id AND
  hbond.struct_id = structure.struct_id AND
  (don_site.HBChemType='hbdon_IME' OR don_site.HBChemType='hbond_IMD') AND
  acc_site.HBChemType='hbacc_CXL' AND
  ABS(don_site.resNum - acc_site.resNum) > 5 AND
  2.7925 < geom.chi AND geom.chi < 3.49065
ORDER BY
  hbond.energy;"


sele <-"
SELECT
  geom.chi, geom.cosBAH,
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
  hbond.acc_id = acc_site.site_id AND
  (acc_site.HBChemType='hbacc_PBA' OR
	acc_site.HBChemType='hbacc_CXL' OR
	acc_site.HBChemType='hbacc_CXA') AND
	ABS(don_site.resNum - acc_site.resNum) > 5;"

f <- query_sample_sources(sample_sources, sele)
f$chi <- (f$chi*180/pi)

## when the BAH angle is straigth (i.e. close to 0), the chi angle is less interesting
#f$BAH_deg <- acos(f$cosBAH)*180/pi
#
#print(summary(f))
#f <- f[f$BAH_deg > 25,]
#print(summary(f))

sc_f <- f[f$acc_chem_type != 'hbacc_PBA',]
dens <- estimate_density_1d_wrap(
  data = sc_f,
  ids = c("sample_source", "don_chem_type"),
  variable = "chi")

# Order the plots better and give more descriptive labels
dens$don_chem_type <- factor(dens$don_chem_type,
	levels=c(
		"hbdon_IMD",
  	"hbdon_GDE",
    "hbdon_HXL",
    "hbdon_PBA",
    "hbdon_IME",
    "hbdon_GDH",
    "hbdon_AHX",
    "hbdon_IND",
    "hbdon_CXA",
    "hbdon_AMO"),
	labels = c(
		"Donor IMD: his",
		"Donor GDE: arg",
		"Donor HXL: ser, thr",
		"Donor PBA: backbone",
		"Donor IME: his",
		"Donor GDH: arg",
		"Donor AHX: tyr",
		"Donor IND: trp",
		"Donor CXA: gln, asn",
		"Donor AMO: lys"))


cos_potential <- function(chi){ -1*(cos( 2 * chi * pi/180 ) + 2 ) / 3 }
potential <- data.frame( x=seq(0,360,length.out=200))
potential$y <- exp(-1*cos_potential(potential$x))/(12.5803 *180/pi)

plot_id = "chi_cosBAH_chem_type_sp2_long_range"
ggplot(data=dens) + theme_bw() +
	geom_line(data=potential, aes(x=x, y=y), size=1.1, color="darkgray") +
	geom_indicator(data=dens, aes(indicator=counts, colour=sample_source, group=sample_source)) +
	geom_line(data=dens, aes(x=x, y=y, colour=sample_source)) +
	facet_wrap( ~ don_chem_type) +
	ggtitle("Hydrogen Bonds CHI Angle for Sidechain sp2 Acceptors with sequence separation at least 6\n(normalized for equal volume per unit distance)") +
	scale_x_continuous('Acceptor Base -- Acceptor Torsion (degrees)', breaks=c(90,270)) +
	scale_y_continuous("Feature Density") +
	theme(legend.position=c(.58,.35)) +
	theme(legend.justification=c("left", "top"))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

#plot_id = "chi_cosBAH_chem_type_sp2_long_range_polar"
#p <- ggplot(data=dens)
#p <- p + geom_line(data=potential, aes(x=x, y=y), size=1.2)
#p <- p + geom_line(aes(x=x, y=y, colour=sample_source))
#p <- p + geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source))
#p <- p + facet_wrap( ~ don_chem_type)
#p <- p + ggtitle("Hydrogen Bonds CHI Angle for Sidechain sp2 Acceptors with sequence separation at least 6\n(normalized for equal volume per unit distance)")
#p <- p + theme_bw()
#p <- p + coord_polar(theta="x")
#p <- p + scale_x_continuous('Acceptor Base -- Acceptor Torsion (degrees)', breaks=c(90,270))
#
#save_plots(self, plot_id, sample_sources, output_dir, output_formats)


dens <- estimate_density_1d_wrap(
  data = f,
  ids = c("sample_source", "acc_chem_type", "don_chem_type"),
  variable = "chi")

plot_id = "chi_cosBAH_chem_type_sp2_by_acc_long_range"
ggplot(data=dens) + theme_bw() +
  geom_line(aes(x=x, y=y, colour=sample_source)) +
  geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
  facet_grid(don_chem_type ~ acc_chem_type) +
  ggtitle("Hydrogen Bonds CHI Angle by Acceptor types with sequence separation at least 6\n(normalized for equal volume per unit distance)") +
  scale_x_continuous('Acceptor Base -- Acceptor Torsion (degrees)', breaks=c(90,270)) +
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
