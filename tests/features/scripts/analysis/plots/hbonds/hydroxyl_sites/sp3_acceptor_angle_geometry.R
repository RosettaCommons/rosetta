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
id = "sp3_acceptor_angle_geometry",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
	acc_atoms.base_x  AS bx,  acc_atoms.base_y  AS by,  acc_atoms.base_z  AS bz,
  acc_atoms.base2_x AS b2x, acc_atoms.base2_y AS b2y, acc_atoms.base2_z AS b2z,
  acc_atoms.atm_x   AS ax,  acc_atoms.atm_y   AS ay,  acc_atoms.atm_z   AS az,
  don_atoms.atm_x   AS hx,  don_atoms.atm_y   AS hy,  don_atoms.atm_z   AS hz,
  CASE acc.resNum - don.resNum
    WHEN -1 THEN '-1' WHEN -2 THEN '-2' WHEN -3 THEN '-3' WHEN -4 THEN '-4'
    WHEN 1 THEN '1,2,3,4' WHEN 2 THEN '1,2,3,4' WHEN 3 THEN '1,2,3,4' WHEN 4 THEN '1,2,3,4'
    ELSE 'long' END AS seq_sep
FROM
  hbonds AS hb,
  hbond_sites AS don,
  hbond_sites AS acc,
  hbond_site_atoms AS don_atoms,
  hbond_site_atoms AS acc_atoms
WHERE
  acc.struct_id = hb.struct_id AND acc.site_id = hb.acc_id AND
  don.struct_id = hb.struct_id AND don.site_id = hb.don_id AND
  don_atoms.struct_id = hb.struct_id AND don_atoms.site_id = hb.don_id AND
  acc_atoms.struct_id = hb.struct_id AND acc_atoms.site_id = hb.acc_id AND
  (acc.HBChemType = 'hbacc_HXL' OR acc.HBChemType = 'hbacc_AHX');"

f <- query_sample_sources(sample_sources, sele)


alt_chi_dihedral <- function(b2, b, a, h){
	alt_b <- (b + b2)/2
	alt_b2 <- vector_crossprod(b - b2, a - b) - alt_b
	vector_dihedral(alt_b2, alt_b, a, h)
}

f <- transform(f,
	vx = (bx + b2x)/2, vy = (by + b2y)/2, vz = (bz + b2z)/2)

cosVAH <- function(b2, b, a, h){
	alt_b <- (b + b2)/2
	vector_dotprod(vector_normalize(a-alt_b), vector_normalize(h-a))
}

f <- transform(f,
  cosBAH = vector_dotprod(
    vector_normalize(cbind(ax-bx, ay-by, az-bz)),
		vector_normalize(cbind(hx-ax, hy-ay, hz-az))),
	cosB2AH = vector_dotprod(
    vector_normalize(cbind(ax-b2x, ay-b2y, az-b2z)),
		vector_normalize(cbind(hx-ax, hy-ay, hz-az))),
	cosVAH = vector_dotprod(
		vector_normalize(cbind(ax-vx, ay-vy, az-vz)),
		vector_normalize(cbind(hx-ax, hy-ay, hz-az))),
	vchi = alt_chi_dihedral(
		cbind(b2x, b2y, b2z), cbind(bx, by, bz),
		cbind(ax, ay, az), cbind(hx, hy, hz)))

print(summary(f))

cosB2AH_dens <- estimate_density_1d(
  data = f, ids = c("sample_source", "seq_sep"), variable = "cosB2AH")
cosBAH_dens <- estimate_density_1d(
  data = f, ids = c("sample_source", "seq_sep"), variable = "cosBAH")
cosVAH_dens <- estimate_density_1d(
  data = f, ids = c("sample_source", "seq_sep"), variable = "cosVAH")

dens_angles <- rbind(
	data.frame(cosVAH_dens, cell="cosVAH"),
	data.frame(cosB2AH_dens, cell="cosB2AH"),
	data.frame(cosBAH_dens, cell="cosBAH"))


d_ply(dens_angles, .(sample_source), function(sub_da){
	ss_id <- sub_da$sample_source[1]
	plot_id <- paste("sp3_acceptor_angle_geometry_by_seq_sep_", ss_id, sep="")
	ggplot(data=sub_da) + theme_bw() +
		geom_line(aes(x=acos(x)*180/pi, y=y, colour=seq_sep)) +
		geom_indicator(aes(colour=seq_sep, indicator=counts, group=seq_sep)) +
		facet_wrap( ~cell, ncol=1) +
		ggtitle(paste("Hydrogen Bond sp3 Acceptors Angle Geometry\nss_id: ", ss_id, sep="")) +
		scale_x_continuous(paste('XXX -- Acceptor -- Donated Hydrogen Angle (degrees)')) +
		scale_y_continuous("FeatureDensity") +
		theme(legend.position="bottom", legend.direction="horizontal")
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
})

plot_id <- "sp3_acceptor_angle_geometry_by_seq_sep"
ggplot(data=dens_angles) + theme_bw() +
	geom_line(aes(x=acos(x)*180/pi, y=y, colour=seq_sep)) +
	geom_indicator(aes(colour=seq_sep, indicator=counts, group=seq_sep)) +
	facet_grid( sample_source ~ cell) +
	ggtitle("Hydrogen Bond sp3 Acceptors Angle Geometry") +
	scale_x_continuous(paste('XXX -- Acceptor -- Donated Hydrogen Angle (degrees)')) +
	scale_y_continuous("FeatureDensity") +
	theme(legend.position="bottom", legend.direction="horizontal")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

})) # end FeaturesAnalysis
