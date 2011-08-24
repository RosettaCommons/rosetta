# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

check_setup()

plot_id <- "SSDists_Helix"

sele <-"
SELECT
  dist.N_N_dist,   dist.N_Ca_dist,  dist.N_C_dist,   dist.N_O_dist,
  dist.Ca_N_dist,  dist.Ca_Ca_dist, dist.Ca_C_dist,  dist.Ca_O_dist,
  dist.C_N_dist,   dist.C_Ca_dist,  dist.C_C_dist,   dist.C_O_dist,
  dist.O_N_dist,   dist.O_Ca_dist,  dist.O_C_dist,   dist.O_O_dist
FROM
  protein_backbone_atom_atom_pairs as dist,
	residue_secondary_structure as r1ss,
	residue_secondary_structure as r2ss
WHERE
  dist.struct_id = r1ss.struct_id AND
  dist.struct_id = r2ss.struct_id AND
  dist.resNum1 = r1ss.resNum AND
  dist.resNum2 = r2ss.resNum AND
  r1ss.dssp = 'H' AND r2ss.dssp = 'H'
LIMIT 10000;"

f <- query_sample_sources(sample_sources, sele)

print(summary(f))

n_n_dists <- function( df ) {
  data.frame( sample_source=c(df$sample_source), at1=c("N"), at2=c("N"), dist=c(df$N_N_dist) )
}
n_ca_dists <- function( df ) {
  data.frame( sample_source=c(df$sample_source,df$sample_source), at1=c("N","N"), at2=c("CA","CA"), dist=c(df$N_Ca_dist, df$Ca_N_dist) )
}
n_c_dists <- function( df ) {
  data.frame( sample_source=c(df$sample_source,df$sample_source), at1=c("N","N"), at2=c("C","C"), dist=c(df$N_C_dist, df$C_N_dist) )
}
n_o_dists <- function( df ) {
  data.frame( sample_source=c(df$sample_source,df$sample_source), at1=c("N","N"), at2=c("O","O"), dist=c(df$N_O_dist, df$O_N_dist) )
}
ca_ca_dists <- function( df ) {
  data.frame( sample_source=c(df$sample_source), at1=c("CA"), at2=c("CA"), dist=c(df$Ca_Ca_dist) )
}
ca_c_dists <- function( df ) {
  data.frame( sample_source=c(df$sample_source,df$sample_source), at1=c("CA","CA"), at2=c("C","C"), dist=c(df$Ca_C_dist, df$C_Ca_dist) )
}
ca_o_dists <- function( df ) {
  data.frame( sample_source=c(df$sample_source,df$sample_source), at1=c("CA","CA"), at2=c("O","O"), dist=c(df$Ca_O_dist, df$O_Ca_dist) )
}
c_c_dists <- function( df ) {
  data.frame( sample_source=c(df$sample_source), at1=c("C"), at2=c("C"), dist=c(df$C_C_dist) )
}
c_o_dists <- function( df ) {
  data.frame( sample_source=c(df$sample_source,df$sample_source), at1=c("C","C"), at2=c("O","O"), dist=c(df$C_O_dist, df$O_C_dist) )
}
o_o_dists <- function( df ) {
  data.frame( sample_source=c(df$sample_source), at1=c("O"), at2=c("O"), dist=c(df$O_O_dist) )
}

f1 = adply( f, 1, n_n_dists )
f2 = adply( f, 1, n_ca_dists )
f3 = adply( f, 1, n_c_dists )
f4 = adply( f, 1, n_o_dists )
f5 = adply( f, 1, ca_ca_dists )
f6 = adply( f, 1, ca_c_dists )
f7 = adply( f, 1, ca_o_dists )
f8 = adply( f, 1, c_c_dists )
f9 = adply( f, 1, c_o_dists )
f10 = adply( f, 1, o_o_dists )

g = rbind( f1, f2, f3, f4, f5, f6, f7, f8, f9, f10 )

print(summary(g))

dens <- estimate_density_1d(
  g, c("sample_source", "at1", "at2"),
  "dist", weight_fun = radial_3d_normalization)

p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x, y=y, colour=sample_source)) +
	geom_indicator(aes(indicator=counts, colour=sample_source)) +
	facet_grid(at1 ~ at2) +
	opts(title = "Backbone atom atom distances involving helical residues\nnormalized for equal weight per unit distance") +
	scale_y_log10("FeatureDensity", limits=c(1e-4,1e0)) +
	scale_x_continuous(expression(paste('Atom Atom Distances (', ring(A), ')')), limits=c(1.5,8), breaks=2:8)

if(nrow(sample_sources) <= 3){
	p <- p + opts(legend.position="bottom", legend.direction="horizontal")
}

save_plots(plot_id, sample_sources, output_dir, output_formats)


#sele <-
#SELECT
#  r1conf.secstruct, r2conf.secstruct,
#  dist.N_N_dist,   dist.N_Ca_dist,  dist.N_C_dist,   dist.N_O_dist,
#  dist.Ca_N_dist,  dist.Ca_Ca_dist, dist.Ca_C_dist,  dist.Ca_O_dist,
#  dist.C_N_dist,   dist.C_Ca_dist,  dist.C_C_dist,   dist.C_O_dist,
#  dist.O_N_dist,   dist.O_Ca_dist,  dist.O_C_dist,   dist.O_O_dist
#FROM
#  hbond_sites as acc_site1,
#  hbond_sites as don_site1,
#  hbonds as hb1,
#  hbond_sites as acc_site2,
#  hbond_sites as don_site2,
#  hbonds as hb2,
#  protein_backbone_atom_atom_pairs as dist,
#	protein_residue_conformation as r1conf,
#	protein_residue_conformation as r2conf
#WHERE
#  acc_site1.struct_id = don_site1.struct_id AND
#  hb1.struct_id = acc_site1.struct_id AND
#  hb1.acc_id = acc_site1.site_id AND
#  hb1.don_id = don_site1.site_id AND
#  acc_site1.resNum = don_site1.resNum + 4 AND
#
#  acc_site2.struct_id = don_site2.struct_id AND
#  hb2.struct_id = acc_site2.struct_id AND
#  hb2.acc_id = acc_site2.site_id AND
#  hb2.don_id = don_site2.site_id AND
#  acc_site2.resNum = don_site2.resNum + 4 AND
#
#  don_site1.resNum = acc_site2.resnum AND
#  hb1.struct_id = hb2.struct_id AND
#
#  dist.struct_id = hb1.struct_id AND
#  ( dist.resNum1 = don_site1.resNum OR dist.resNum2 = don_site1.resNum )
#
#LIMIT 10000;

