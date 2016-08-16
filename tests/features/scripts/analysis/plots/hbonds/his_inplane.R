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
id = "his_inplane",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("StructureFeatures", "HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
  structure.tag,
  don_site_pdb.chain, don_site_pdb.resNum, don_site_pdb.iCode,
  acc_site_pdb.chain, acc_site_pdb.resNum, acc_site_pdb.iCode,
  acc_site.resNum,
  hbond.energy,
  don_his_site.resNum,
  geom.AHdist, geom.chi,
  acc_atoms.base_x     AS p1x, acc_atoms.base_y     AS p1y, acc_atoms.base_z     AS p1z,
  acc_atoms.atm_x      AS p2x, acc_atoms.atm_y      AS p2y, acc_atoms.atm_z      AS p2z,
  don_his_atoms.base_x AS p3x, don_his_atoms.base_y AS p3y, don_his_atoms.base_z AS p3z,
  acc_his_atoms.atm_x  AS p4x, acc_his_atoms.atm_y  AS p4y, acc_his_atoms.atm_z  AS p4z,
  acc_site.HBChemType AS acc_chem_type
FROM
  hbond_geom_coords AS geom,
  hbonds AS hbond,
  hbond_sites AS don_his_site,
  hbond_sites AS acc_his_site,
  hbond_sites AS acc_site,
  hbond_site_atoms AS don_his_atoms,
  hbond_site_atoms AS acc_his_atoms,
  hbond_site_atoms AS acc_atoms,
  hbond_sites_pdb AS don_site_pdb,
  hbond_sites_pdb AS acc_site_pdb,
  structures AS structure
WHERE
  hbond.struct_id = geom.struct_id AND hbond.hbond_id =  geom.hbond_id AND
  hbond.struct_id = acc_site.struct_id AND hbond.acc_id = acc_site.site_id AND
  hbond.struct_id = don_his_site.struct_id AND hbond.don_id = don_his_site.site_id AND
  don_his_site.struct_id = acc_his_site.struct_id AND don_his_site.resNum = acc_his_site.resNum AND
  ( don_his_site.HBChemType = 'hbdon_IME' OR don_his_site.HBChemType = 'hbdon_IMD' ) AND
  ( acc_his_site.HBChemType = 'hbacc_IME' OR acc_his_site.HBChemType = 'hbacc_IMD' ) AND
  acc_site.HBChemType = 'hbacc_CXL' AND
  hbond.struct_id = don_site_pdb.struct_id AND don_site_pdb.site_id = don_his_site.site_id AND
  hbond.struct_id = acc_site_pdb.struct_id AND acc_site_pdb.site_id = acc_site.site_id AND
  structure.struct_id = hbond.struct_id AND
  don_his_site.struct_id = don_his_atoms.struct_id AND don_his_site.site_id == don_his_atoms.site_id AND
  acc_his_site.struct_id = acc_his_atoms.struct_id AND acc_his_site.site_id == acc_his_atoms.site_id AND
  acc_site.struct_id = acc_atoms.struct_id AND acc_site.site_id = acc_atoms.site_id;"


f <- query_sample_sources(sample_sources, sele)
print(summary(f))

# TODO REPLACE THIS MATH WITH THE VERSIONS IN scripts/methods/vector_math.R

normalize <- function(p){
  p/sqrt(sum(p*p))
}

print(normalize(c(1,1,1)))

cross <- function(p,q){
  c(x=(p[2] * q[3] ) - (p[3] * q[2] ),y=(p[3] * q[1] ) - (p[1] * q[3] ),z=(p[1] * q[2] ) - (p[2] * q[1] ))
}

dihedral <- function(df){
  p1 <- df[1,c("p1x", "p1y", "p1z")]
  p2 <- df[1,c("p2x", "p2y", "p2z")]
  p3 <- df[1,c("p3x", "p3y", "p3z")]
  p4 <- df[1,c("p4x", "p4y", "p4z")]

#  print(paste("p1", p1))
#  print(paste("p2", p2))
#  print(paste("p3", p3))
#  print(paste("p4", p4))

  a <- normalize(p2-p1)
  b <- normalize(p3-p2)
  c <- normalize(p4-p3)

#  print(paste("a", a))
#  print(paste("b", b))
#  print(paste("c", c))

  x = -sum(a*c) + (sum(a*b) * sum(b*c))
#  print(paste("x", x))

#  print(cross(b,c))

  y = sum(a*cross(b,c))
#  print(paste("y", y))

  angle <- atan2(y,x)
}

df <- data.frame(p1x=0, p1y=0, p1z=0, p2x=1, p2y=0, p2z=0, p3x=1, p3y=1, p3z=0, p4x=1, p4y=1, p4z=1)
print(paste("test dihedral", dihedral(df)))


p1 <- data.frame(x=c(0), y=c(0), z=c(0))
p2 <- data.frame(x=c(0), y=c(1), z=c(0))



f <- adply(f, 1, function(df){
  angle <- dihedral(df)
  data.frame(df, eta=angle*180/pi)
})

examples <- f[ f$eta < 270+20-180 & f$eta > 270 -20-180,]
write.csv(examples, "/tmp/eta_his_HXL_examples.csv")


dens <- estimate_density_1d_wrap(
  data = f,
  ids = c("sample_source"),
  variable = "eta")

print(summary(dens))

plot_id = "eta_his_CXL"
ggplot(data=dens) + theme_bw() +
  geom_line(aes(x=x, y=y, colour=sample_source)) +
  geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
  opts(title = "Hydrogen Bonds eta Angle HIS-Epsilon Donor to Charged Acceptor") +
  scale_x_continuous('Acceptor -- Donor Torsion (degrees)') +
  scale_y_continuous('Feature Density', limits=c(0,0.006))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

sele <-"
SELECT
  structure.tag,
  don_site_pdb.chain, don_site_pdb.resNum, don_site_pdb.iCode,
  acc_site_pdb.chain, acc_site_pdb.resNum, acc_site_pdb.iCode,
  acc_site.resNum,
  hbond.energy,
  don_his_site.resNum,
  geom.AHdist, geom.chi,
  acc_atoms.base_x     AS p1x, acc_atoms.base_y     AS p1y, acc_atoms.base_z     AS p1z,
  acc_atoms.atm_x      AS p2x, acc_atoms.atm_y      AS p2y, acc_atoms.atm_z      AS p2z,
  don_his_atoms.base_x AS p3x, don_his_atoms.base_y AS p3y, don_his_atoms.base_z AS p3z,
  acc_his_atoms.atm_x  AS p4x, acc_his_atoms.atm_y  AS p4y, acc_his_atoms.atm_z  AS p4z,
  acc_site.HBChemType AS acc_chem_type
FROM
  hbond_geom_coords AS geom,
  hbonds AS hbond,
  hbond_sites AS don_his_site,
  hbond_sites AS acc_his_site,
  hbond_sites AS acc_site,
  hbond_site_atoms AS don_his_atoms,
  hbond_site_atoms AS acc_his_atoms,
  hbond_site_atoms AS acc_atoms,
  hbond_sites_pdb AS don_site_pdb,
  hbond_sites_pdb AS acc_site_pdb,
  structures AS structure
WHERE
  hbond.struct_id = geom.struct_id AND hbond.hbond_id =  geom.hbond_id AND
  hbond.struct_id = acc_site.struct_id AND hbond.acc_id = acc_site.site_id AND
  hbond.struct_id = don_his_site.struct_id AND hbond.don_id = don_his_site.site_id AND
  don_his_site.struct_id = acc_his_site.struct_id AND don_his_site.resNum = acc_his_site.resNum AND
  ( don_his_site.HBChemType = 'hbdon_IME' OR don_his_site.HBChemType = 'hbdon_IMD' ) AND
  ( acc_his_site.HBChemType = 'hbacc_IME' OR acc_his_site.HBChemType = 'hbacc_IMD' ) AND
  acc_site.HBChemType = 'hbacc_CXL' AND
  hbond.struct_id = don_site_pdb.struct_id AND don_site_pdb.site_id = don_his_site.site_id AND
  hbond.struct_id = acc_site_pdb.struct_id AND acc_site_pdb.site_id = acc_site.site_id AND
  structure.struct_id = hbond.struct_id AND
  don_his_site.struct_id = don_his_atoms.struct_id AND don_his_site.site_id == don_his_atoms.site_id AND
  acc_his_site.struct_id = acc_his_atoms.struct_id AND acc_his_site.site_id == acc_his_atoms.site_id AND
  acc_site.struct_id = acc_atoms.struct_id AND acc_site.site_id = acc_atoms.site_id AND
  geom.chi > 2.61799388 AND geom.chi < 3.66519143;"  # 150-210 degrees

f <- query_sample_sources(sample_sources, sele)
print(summary(f))

df <- data.frame(p1x=0, p1y=0, p1z=0, p2x=1, p2y=0, p2z=0, p3x=1, p3y=1, p3z=0, p4x=1, p4y=1, p4z=1)
print(paste("test dihedral", dihedral(df)))

p1 <- data.frame(x=c(0), y=c(0), z=c(0))
p2 <- data.frame(x=c(0), y=c(1), z=c(0))



f <- adply(f, 1, function(df){
  angle <- dihedral(df)
  data.frame(df, eta=angle*180/pi)
})

examples <- f[ f$eta < 270+20-180 & f$eta > 270 -20-180,]
write.csv(examples, "/tmp/eta_his_HXL_examples.csv")


dens <- estimate_density_1d_wrap(
  data = f,
  ids = c("sample_source"),
  variable = "eta")

print(summary(dens))

plot_id = "eta_his_CXL_to_syn_orbital"
ggplot(data=dens) + theme_bw() +
  geom_line(aes(x=x, y=y, colour=sample_source)) +
  geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
  opts(title = "Hydrogen Bonds eta Angle HIS Donor to Charged Acceptor w/ 150<chi<210") +
  scale_x_continuous('Acceptor -- Donor Torsion (degrees)') +
  scale_y_continuous('Feature Density', limits=c(0,0.006))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
