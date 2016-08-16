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
id = "cosB2AH_sp3_acceptors",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
  acc_atoms.base2_x AS b2x, acc_atoms.base2_y AS b2y, acc_atoms.base2_z AS b2z,
  acc_atoms.atm_x   AS ax,  acc_atoms.atm_y   AS ay,  acc_atoms.atm_z   AS az,
  don_atoms.atm_x   AS hx,  don_atoms.atm_y   AS hy,  don_atoms.atm_z   AS hz
FROM
  hbonds AS hb,
  hbond_sites AS acc,
  hbond_site_atoms AS don_atoms,
  hbond_site_atoms AS acc_atoms
WHERE
  acc.struct_id = hb.struct_id AND acc.site_id = hb.acc_id AND
  don_atoms.struct_id = hb.struct_id AND don_atoms.site_id = hb.don_id AND
  acc_atoms.struct_id = hb.struct_id AND acc_atoms.site_id = hb.acc_id AND
  (acc.HBChemType = 'hbacc_HXL' OR acc.HBChemType = 'hbacc_AHX');"

f <- query_sample_sources(sample_sources, sele)

f <- transform(f,
  B2AH_proj = vector_dotprod(
    vector_normalize(cbind(ax-b2x, ay-b2y, az-b2z)), vector_normalize(cbind(hx-ax, hy-ay, hz-az))))

dens <- estimate_density_1d(
  data = f, ids = c("sample_source"), variable = "B2AH_proj")

plot_id = "cosB2AH_sp3_acceptors"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=acos(x)*180/pi, y=y, colour=sample_source)) +
	geom_indicator(aes(colour=sample_source, indicator=counts, group=sample_source)) +
	ggtitle("Hydrogen Bonds B2AH Projection for sp3 Acceptors") +
	scale_x_continuous("cos( Base2 -- Acceptor--Hydrogen Angle" ) +
	scale_y_continuous("FeatureDensity")
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
