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
id = "hdist_out_of_sp2plane",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
  hbond.energy,
  geom.AHdist, geom.chi,
  acc_site.HBChemType AS acc_chem_type,
  don_site.HBChemType AS don_chem_type,
	acc_atoms.base2_x AS ab2x, acc_atoms.base2_y AS ab2y, acc_atoms.base2_z AS ab2z, -- acceptor base 2 atom
 	acc_atoms.base_x  AS abx,  acc_atoms.base_y  AS aby,  acc_atoms.base_z  AS abz,  -- acceptor base atom
	acc_atoms.atm_x   AS ax,   acc_atoms.atm_y   AS ay,   acc_atoms.atm_z   AS az,   -- acceptor atom
	don_atoms.atm_x   AS hx,   don_atoms.atm_y   AS hy,   don_atoms.atm_z   AS hz    -- hydrogen atom
FROM
  hbond_geom_coords AS geom,
  hbonds AS hbond,
  hbond_sites AS don_site, hbond_sites AS acc_site,
  hbond_site_atoms AS don_atoms, hbond_site_atoms AS acc_atoms
WHERE
  hbond.struct_id = geom.struct_id AND hbond.hbond_id =  geom.hbond_id AND
  hbond.struct_id = acc_site.struct_id AND hbond.acc_id = acc_site.site_id AND
  hbond.struct_id = don_site.struct_id AND hbond.don_id = don_site.site_id AND
  abs( acc_site.resNum - don_site.resNum ) > 5 AND
  don_atoms.struct_id = hbond.struct_id AND don_atoms.site_id = hbond.don_id AND
  acc_atoms.struct_id = hbond.struct_id AND acc_atoms.site_id = hbond.acc_id AND
  (acc_site.HBChemType = 'hbacc_CXA' OR acc_site.HBChemType = 'hbacc_CXL' OR
		acc_site.HBChemType = 'hbacc_PBA');"
f <- query_sample_sources(sample_sources, sele)

ooplane <- function(ab2, ab, a, h){
	w1 <- vector_normalize(ab2-ab)
	w2 <- vector_normalize(a-ab)
	w3 <- vector_normalize(h-a)
	vector_dotprod(w3, vector_crossprod(w1, w2))
}

f$doop <- with(f, ooplane(
	cbind(ab2x, ab2y, ab2z), cbind(abx, aby, abz),
	cbind(ax, ay, az), cbind(hx, hy, hz)))

f$don_chem_type <- factor(f$don_chem_type,
	levels=c("hbdon_IMD", "hbdon_GDE", "hbdon_HXL", "hbdon_PBA",
		"hbdon_IME", "hbdon_GDH", "hbdon_AHX", "hbdon_IND",
		"hbdon_CXA", "hbdon_AMO"),
	labels = c("dIMD: h", "dGDE: r", "dHXL: s,t", "dPBA: bb",
		"dIME: h", "dGDH: r", "dAHX: y", "dIND: w",
		"dCXA: n,q", "dAMO: k"))

dens <- estimate_density_1d(f, c("sample_source", "don_chem_type"), "doop")

plot_id = "hdist_out_of_sp2plane"
ggplot(data=dens) + theme_bw() +
  geom_line(aes(x=x, y=y, colour=sample_source)) +
  geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
	facet_wrap( ~ don_chem_type ) +
  ggtitle("Hydrogen Distance out of the SP2 Plane") +
  scale_x_continuous('HDist out of plane') +
  scale_y_continuous('Feature Density') +
	theme(legend.position=c(.58,.35)) +
	theme(legend.justification=c("left", "top"))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
