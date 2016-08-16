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
id = "bb_by_pos_seq_sep",
author = "Matthew O'Meara",
brief_description = "",
long_description = "
   Following IUPAC-IUB 1970 Rule 6.3, a pi-helix is a regular 4.4_{16} helix consisting of at least 2 i<-i+5 hydrogen bonds. The first helical residue is the first whose CO group is regularly hydrogen bonded to NH along the helix and the last helical residue is the last whose NH is regularly hydrogen bonded to the CO along the helix.

  In establishing where a pi-helix begins and ends, DSSP takes 'regular' to mean not bifurcated, while Fodje, and Al-Kardaghi allow bifurcation as long as the strongest one is forming the i<i+5 interaction.

M. N. Fodje, S. Al-Kardaghi, Occurrence, conformational features and amino acid propensities for the pi-helix. Protein Engineering vol. 15 no. 5 pp. 353-358, 2002.

IUPAC-IUB Comm. on Biochem. Nomenclature, IUPAC-IUB Commission on Biochemical Nomenclature. Abbreviations and symbols for the description of the conformation of polypeptide chains. Tentative rules (1969), Biochemistry, 1970, 9 (18), pp 3471-3479.

W. Kabsch, C. Sander, Dictionary of protein secondary structure: pattern recognition of hydrogen-bonded and geometrical features. Biopolymers. 1983 22 2577-2637. PMID: 6667333; UI: 84128824.

Questions:
  * What is the dssp distribution of backbone residues participating in an i<-i+5 or i-5<-i hydrogen bond?
  * What fraction of backbone-backbone i<-i+5 hydrogen bonds are part of pi-helices?
  * Do the hydrogen bonds that form pi-helices have specific geometry?
  * Does backbone-backbone i<-i+5 residues have specific geometry?
  * Does Rosetta form pi-helices too often or to rarely?
  * Does Rosetta recapitulate the hydrogen bond geometry that form pi-helices?
  * Does Rosetta recapitulate the geometry of backbone-backbone i<-i+5 hydrogen bonds?",

feature_reporter_dependencies = c("ResidueSecondaryStructureFeatures", "HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
	geom.AHdist, geom.cosBAH, geom.cosAHD, geom.chi,
	don_dssp.dssp AS don_dssp, acc_dssp.dssp AS acc_dssp,
	don.resNum - acc.resNum AS seq_sep,
	CASE don.resNum - acc.resNum
		WHEN -4 THEN 'M4' WHEN -3 THEN 'M3'
		WHEN -2 THEN 'M2' WHEN -1 THEN 'M1'
		WHEN 0  THEN '0'  WHEN 1 THEN 'P1'
		WHEN 2  THEN 'P2' WHEN 3 THEN 'P4'
		WHEN 4  THEN 'P4' WHEN 5 THEN 'P5'
		WHEN 6  THEN 'P6' WHEN 7 THEN 'P7'
		ELSE 'long' END AS seq_sep_class
FROM
	hbonds AS hb,
	hbond_geom_coords AS geom,
	hbond_sites AS don, hbond_sites AS acc,
	hbond_sites_pdb AS don_pdb, hbond_sites_pdb AS acc_pdb,
	residue_secondary_structure AS don_dssp,
	residue_secondary_structure AS acc_dssp
WHERE
	geom.struct_id = hb.struct_id AND geom.hbond_id = hb.hbond_id AND
	don.struct_id = hb.struct_id AND don.site_id = hb.don_id AND
	don.HBChemType = 'hbdon_PBA' AND
	acc.struct_id = hb.struct_id AND acc.site_id = hb.acc_id AND
	acc.HBChemType = 'hbacc_PBA' AND
	don_pdb.struct_id = hb.struct_id AND don_pdb.site_id = hb.don_id AND
	don_pdb.heavy_atom_temperature < 30 AND
	acc_pdb.struct_id = hb.struct_id AND acc_pdb.site_id = hb.acc_id AND
	acc_pdb.heavy_atom_temperature < 30 AND
	don_dssp.struct_id = hb.struct_id AND don_dssp.resNum = don.resNum AND
	acc_dssp.struct_id = hb.struct_id AND acc_dssp.resNum = acc.resNum AND
	don.resNum - acc.resNum >= 4;"

f <- query_sample_sources(sample_sources, sele)

dens <- estimate_density_1d(
	f, c("sample_source", "seq_sep_class"),
	"AHdist", weight_fun = radial_3d_normalization, adjust=.5)

plot_id <- "hbond_backbone_backbone_seq_sep5_vs_long_range_AHdist"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x, y=y, colour=seq_sep_class)) +
	geom_indicator(aes(indicator=counts, colour=seq_sep_class, group=seq_sep_class)) +
	facet_wrap( ~ sample_source) +
	ggtitle("Backbone-Backbone H-Bonds A-H Distance by Sequence Separation Classes\nB-Factor < 30 normalized for equal weight per unit distance") +
	scale_y_continuous("FeatureDensity", limits=c(0,6), breaks=c(1,3,5)) +
	scale_x_continuous(expression(paste('Acceptor -- Proton Distance (', ring(A), ')')), limits=c(1.4,2.7), breaks=c(1.6, 1.9, 2.2, 2.6))

if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)




dens <- estimate_density_1d_reflect_boundary(
	f, c("sample_source", "seq_sep_class"),
	"cosAHD", reflect_right=T, right_boundary=1, adjust=.5)

plot_id <- "hbond_backbone_backbone_seq_sep5_vs_long_range_cosAHD"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=180-acos(x)*180/pi, y=y, colour=seq_sep_class)) +
	geom_indicator(aes(indicator=counts, colour=seq_sep_class, group=seq_sep_class)) +
	facet_wrap( ~ sample_source) +
	ggtitle("Backbone-Backbone H-Bonds AHD Angle by Sequence Separation Classes\nB-Factor < 30 normalized for equal weight per unit distance") +
	scale_y_continuous("FeatureDensity", limits=c(0,20), breaks=c(0,5,10,15)) +
	scale_x_continuous("Acceptor -- Hydrogen -- Donor (degrees)", trans="reverse")

if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)



dens <- estimate_density_1d(
	f, c("sample_source", "seq_sep_class"), "cosBAH", adjust=.5)

plot_id <- "hbond_backbone_backbone_seq_sep5_vs_long_range_cosBAH"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=acos(x)*180/pi, y=y, colour=seq_sep_class)) +
	geom_indicator(aes(indicator=counts, colour=seq_sep_class, group=seq_sep_class)) +
	facet_wrap( ~ sample_source) +
	ggtitle("Backbone-Backbone H-Bonds BAH Angle by Sequence Separation Classes\nB-Factor < 30 normalized for equal weight per unit distance") +
	scale_x_continuous(paste('Base -- Acceptor -- Hydrogen (degrees)')) +
	scale_y_continuous("FeatureDensity")

if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)


f$chi_deg <- f$chi*180/pi
dens <- estimate_density_1d_wrap(
	f, c("sample_source", "seq_sep_class"), "chi_deg", adjust=.5)

plot_id <- "hbond_backbone_backbone_seq_sep5_vs_long_range_chi"
p <- ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=x, y=y, colour=seq_sep_class)) +
	geom_indicator(aes(indicator=counts, colour=seq_sep_class, group=seq_sep_class)) +
	facet_wrap( ~ sample_source) +
	ggtitle("Backbone-Backbone H-Bonds CHI Angle by Sequence Separation Classes\nB-Factor < 30 normalized for equal weight per unit distance") +
	scale_x_continuous('Acceptor Base -- Acceptor Torsion (degrees)', breaks=c(90,270)) +
	scale_y_continuous("FeatureDensity")

if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)



f <- transform(f,
	capx = 2*sin(acos(cosBAH)/2)*cos(chi),
	capy = 2*sin(acos(cosBAH)/2)*sin(chi))

capx_limits <- range(f$capx)
capy_limits <- range(f$capy)
f <- ddply(f, c("sample_source", "seq_sep_class"),
	transform, counts = length(sample_source))

plot_id <- "hbond_backbone_backbone_seq_sep5_vs_long_range_chi_BAH"
d_ply(f, .(sample_source), function(sub_f){
	ss_id <- sub_f[1,"sample_source"]

	p <- ggplot(data=sub_f) + theme_bw() +
		theme(panel.background=element_rect(fill="#00007F", colour="#00007F")) +
		stat_bin2d(aes(x=capx, y=capy, fill=log(..density..)), binwidth=c(.01, .01)) +
		polar_equal_area_grids_bw() +
		geom_indicator(aes(indicator=counts, colour="white")) +
		coord_equal(ratio=1) +
		facet_wrap( ~ seq_sep_class) +
		ggtitle(paste("Backbone-Backbone H-Bonds CHI Angle by Sequence Separation Classes\nB-Factor < 30 equal coordinate projection Sample Source: ", ss_id, sep="")) +
		scale_x_continuous('Longitude: Acceptor Base -- Acceptor Torsion (degrees) ', limits=capx_limits, breaks=c()) +
		scale_y_continuous("Latitude: Acceptor Base -- Acceptor -- Hydrogen (degrees)", limits=capy_limits, breaks=c())

	if(nrow(sample_sources) <= 3){
		p <- p + theme(legend.position="bottom", legend.direction="horizontal")
	}

	save_plots(self, plot_id, sample_sources[sample_sources$sample_source == ss_id,], output_dir, output_formats)
})

})) # end FeaturesAnalysis
