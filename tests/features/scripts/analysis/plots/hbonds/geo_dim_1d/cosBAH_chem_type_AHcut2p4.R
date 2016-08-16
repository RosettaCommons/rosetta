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
id = "cosBAH_chem_type_AHcut2p4",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
	acc_atoms.base_x AS bx, acc_atoms.base_y AS by, acc_atoms.base_z AS bz,
	acc_atoms.atm_x  AS ax, acc_atoms.atm_y  AS ay, acc_atoms.atm_z  AS az,
	don_atoms.atm_x  AS hx, don_atoms.atm_y  AS hy, don_atoms.atm_z  AS hz,
  acc.HBChemType AS acc_chem_type, don.HBChemType AS don_chem_type
FROM
  hbond_geom_coords AS geom,
  hbonds AS hb,
  hbond_sites AS don, hbond_sites AS acc,
  hbond_site_atoms AS don_atoms, hbond_site_atoms AS acc_atoms
WHERE
	geom.AHdist < 2.4 AND
  geom.struct_id = hb.struct_id AND geom.hbond_id = hb.hbond_id AND
  don.struct_id = hb.struct_id AND don.site_id = hb.don_id AND
  acc.struct_id = hb.struct_id AND acc.site_id = hb.acc_id AND
	don_atoms.struct_id = hb.struct_id AND don_atoms.site_id = hb.don_id AND
	acc_atoms.struct_id = hb.struct_id AND acc_atoms.site_id = hb.acc_id;"

f <- query_sample_sources(sample_sources, sele)

f <- transform(f,
	cosBAH = vector_dotprod(
		vector_normalize(cbind(ax-bx, ay-by, az-bz)),
		vector_normalize(cbind(hx-ax, hy-ay, hz-az))))


# This is deprecated please use the hbond_chem_types table for the lables instead
# Order the plots better and give more descriptive labels
f$don_chem_type <- factor(f$don_chem_type,
	levels = c("hbdon_IMD", "hbdon_IME", "hbdon_GDE", "hbdon_GDH",
		"hbdon_AHX", "hbdon_HXL", "hbdon_IND", "hbdon_AMO", "hbdon_CXA", "hbdon_PBA"),
	labels = c("dIMD: h", "dIME: h", "dGDE: r", "dGDH: r",
		"dAHX: y", "dHXL: s,t", "dIND: w", "dAMO: k", "dCXA: n,q", "dPBA: bb"))

# This is deprecated please use the hbond_chem_types table for the lables instead
# Order the plots better and give more descriptive labels
f$acc_chem_type <- factor(f$acc_chem_type,
	levels = c("hbacc_IMD", "hbacc_IME", "hbacc_AHX", "hbacc_HXL",
		"hbacc_CXA", "hbacc_CXL", "hbacc_PBA"),
	labels = c("aIMD: h", "aIME: h", "aAHX: y", "aHXL: s,t",
		"aCXA: n,q", "aCXL: d,e", "aPBA: bb"))


dens <- estimate_density_1d(
  data = f,
  ids = c("sample_source", "acc_chem_type", "don_chem_type"),
  variable = "cosBAH")

plot_id = "cosBAH_chem_type_AHdis2p4cut"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=acos(x)*180/pi, y=-log(y), colour=sample_source)) +
	geom_indicator(aes(indicator=counts)) +
	facet_grid( don_chem_type ~ acc_chem_type ) +
	ggtitle("Hydrogen Bonds BAH Angle by Chemical Type, AHdist < 2.4\n(normalized for equal volume per unit distance)") +
	labs(x=expression(paste('Base -- Acceptor -- Hydrogen (degrees)')),
	     y="log(FeatureDensity + 1)") +
	scale_y_continuous(limits=c(-2.3,6))
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
