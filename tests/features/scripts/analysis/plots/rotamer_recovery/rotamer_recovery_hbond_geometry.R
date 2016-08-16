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
id = "rotamer_recovery_hbond_geometry",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("ResidueFeatures", "RotamerRecoveryFeatures", "ResidueSecondaryStructureFeatures", "PdbDataFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


sele <-"
CREATE TEMPORARY TABLE IF NOT EXISTS chi_recovered (
	struct_id INTEGER,
	resNum INTEGER,
	first,
	PRIMARY KEY(struct_id, resNum));


INSERT INTO chi_recovered
SELECT
	res.struct_id AS struct_id, res.resNum,
	rr.divergence < 2 AS second
FROM
	residues AS res,
	residue_pdb_confidence AS bb,
	rotamer_recovery AS rr,
	residue_secondary_structure AS ss, dssp_codes AS dssp_code
WHERE
	rr.struct_id = res.struct_id AND rr.resNum = res.resNum AND
	bb.struct_id = res.struct_id AND bb.residue_number = res.resNum AND
	bb.max_sc_temperature < 30 AND
	ss.struct_id = res.struct_id AND ss.resNum = res.resNum AND
	dssp_code.code = ss.dssp AND res.name3 = 'GLU';

CREATE INDEX IF NOT EXISTS hbonds_struct_id_acc_id ON hbonds(
	struct_id, acc_id);

SELECT
	rr.first,
	geo.AHdist, geo.cosBAH, geo.cosAHD, geo.chi,
	acc.HBChemType AS acc_chem_type, don.HBChemType AS don_chem_type,
	CASE don.resNum - acc.resNum
	WHEN -1 THEN '-1' WHEN -2 THEN '-2' WHEN -3 THEN '-3' WHEN -4 THEN '-4'
	WHEN 1 THEN '1' WHEN 2 THEN '2' WHEN 3 THEN '3' WHEN 4 THEN '4'
	ELSE 'long' END AS seq_sep
FROM
	hbonds AS hb,
	hbond_sites AS acc,
	chi_recovered AS rr,
	hbond_sites AS don,
	hbond_geom_coords AS geo
WHERE
	geo.struct_id = hb.struct_id AND geo.hbond_id = hb.hbond_id AND
	acc.struct_id = hb.struct_id AND acc.site_id = hb.acc_id AND
	don.struct_id = hb.struct_id AND don.site_id = hb.don_id AND
	acc.HBChemType = 'hbacc_CXL' AND
	rr.struct_id = acc.struct_id AND rr.resNum = acc.resNum;"

#ANALYZE;
#
#SELECT
#	hb.struct_id,
#	geo.AHdist, geo.cosBAH, geo.cosAHD, geo.chi,
#	acc.resNum,
#	acc.HBChemType AS acc_chem_type, don.HBChemType AS don_chem_type,
#	CASE don.resNum - acc.resNum
#	WHEN -1 THEN '-1' WHEN -2 THEN '-2' WHEN -3 THEN '-3' WHEN -4 THEN '-4'
#	WHEN 1 THEN '1' WHEN 2 THEN '2' WHEN 3 THEN '3' WHEN 4 THEN '4'
#	ELSE 'long' END AS seq_sep
#FROM
#	hbonds AS hb,
#	hbond_sites AS acc,
#	hbond_sites AS don,
#	hbond_geom_coords AS geo
#WHERE
#	geo.struct_id = hb.struct_id AND geo.hbond_id = hb.hbond_id AND
#	acc.struct_id = hb.struct_id AND acc.site_id = hb.acc_id AND
#	don.struct_id = hb.struct_id AND don.site_id = hb.don_id AND
#	acc.HBChemType = 'hbacc_CLX' AND
#	rr.struct_id = acc.struct_id AND rr.resNum = acc.resNum;"
#

f <- query_sample_sources(sample_sources, sele)

# Order the plots better and give more descriptive labels
f$don_chem_type <- factor(f$don_chem_type,
	levels = c("hbdon_IMD", "hbdon_IME", "hbdon_GDE", "hbdon_GDH",
		"hbdon_AHX", "hbdon_HXL", "hbdon_IND", "hbdon_AMO", "hbdon_CXA", "hbdon_PBA"),
	labels = c("dIMD: h", "dIME: h", "dGDE: r", "dGDH: r",
		"dAHX: y", "dHXL: s,t", "dIND: w", "dAMO: k", "dCXA: n,q", "dPBA: bb"))

# Order the plots better and give more descriptive labels
f$acc_chem_type <- factor(f$acc_chem_type,
	levels = c("hbacc_IMD", "hbacc_IME", "hbacc_AHX", "hbacc_HXL",
		"hbacc_CXA", "hbacc_CXL", "hbacc_PBA"),
	labels = c("aIMD: h", "aIME: h", "aAHX: y", "aHXL: s,t",
		"aCXA: n,q", "aCXL: d,e", "aPBA: bb"))

f$first <- factor(f$first, levels=c(1,0), labels=c("Recovered", "Not Recovered"))
plot_parts <- list(
	theme_bw(),
	geom_line(aes(x=x, y=y, colour=sample_source, size=first)),
	geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)),
	scale_y_continuous("FeatureDensity)", limits=c(0,6), breaks=c(1,3,5)),
	scale_size_manual("first chi angle", values=c(.5, 1)))


plot_id <- "hbond_GLU_helical_rotamer_recovery_AHdist_by_don_chem_type"
dens <- estimate_density_1d(
	f, c("sample_source", "don_chem_type", "first"),
	"AHdist", weight_fun = radial_3d_normalization)

p <- ggplot(data=dens) + plot_parts +
	ggtitle("HBond A-H Distance for GLU by Donor Type, B-Factor < 30\nnormalized for equal weight per unit distance") +
	facet_wrap(~don_chem_type) +
	scale_x_continuous(expression(paste('Acceptor -- Proton Distance (', ring(A), ')')), limits=c(1.4,2.7), breaks=c(1.6, 1.9, 2.2, 2.6))

if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)


plot_id <- "hbond_GLU_helical_rotamer_recovery_AHdist_by_seq_sep"
dens <- estimate_density_1d(
	f, c("sample_source", "seq_sep", "first"),
	"AHdist", weight_fun = radial_3d_normalization)

p <- ggplot(data=dens) + plot_parts +
	ggtitle("HBond A-H Distance for GLU by Sequence Separation, B-Factor < 30\nnormalized for equal weight per unit distance") +
	facet_wrap(~seq_sep) +
	scale_x_continuous(expression(paste('Acceptor -- Proton Distance (', ring(A), ')')), limits=c(1.4,2.7), breaks=c(1.6, 1.9, 2.2, 2.6))

if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)


plot_id <- "hbond_GLU_helical_rotamer_recovery_AHdist_by_seq_sep_don_chem_type"
dens <- estimate_density_1d(
	f, c("sample_source", "seq_sep", "first", "don_chem_type"),
	"AHdist", weight_fun = radial_3d_normalization)

p <- ggplot(data=dens) + plot_parts +
	ggtitle("HBond A-H Distance for GLU by Sequence Separation and Donor Type, B-Factor < 30\nnormalized for equal weight per unit distance") +
	facet_grid(seq_sep ~ don_chem_type) +
	scale_x_continuous(expression(paste('Acceptor -- Proton Distance (', ring(A), ')')), limits=c(1.4,2.7), breaks=c(1.6, 1.9, 2.2, 2.6))

if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)


f$first_indicator <- f$first == "Recovered"
g <- cast(f, seq_sep + sample_source ~ don_chem_type, margins=T, value="first_indicator", mean)
ascii(g, header=F, digits=4)

#|============================================================================================================================================================== 
#1.1+| | seq_sep | sample_source             | dIMD: h | dIME: h | dGDE: r | dGDH: r | dAHX: y | dHXL: s,t | dIND: w | dAMO: k | dCXA: n,q | dPBA: bb | (all)  
#| 1  | 1       | top8000_r45890_111114     |         | 0.8889  | 0.5122  | 0.6957  | 0.6923  | 0.5294    | 0.6667  | 0.7143  | 0.5714    | 0.7072   | 0.6874 
#| 2  | 1       | top8000_olf_r45890_111114 |         | 0.8333  | 0.5625  | 0.6415  | 0.8333  | 0.6250    |         | 0.5897  | 0.5676    | 0.6544   | 0.6395 
#
#| 4  | -1      | top8000_r45890_111114     | 0.5000  | 1.0000  | 0.7931  | 0.8158  | 0.4000  | 0.5514    | 1.0000  | 0.7386  | 0.7193    | 0.6071   | 0.7019 
#| 5  | -1      | top8000_olf_r45890_111114 | 0.4444  | 1.0000  | 0.8100  | 0.8302  | 0.6667  | 0.5204    |         | 0.6944  | 0.8108    | 0.5938   | 0.7118 
#
#| 7  | 2       | top8000_r45890_111114     | 0.6364  | 0.8827  | 0.8477  | 0.8369  | 0.9157  | 0.7517    | 0.9008  | 0.8079  | 0.8503    | 0.7696   | 0.8228 
#| 8  | 2       | top8000_olf_r45890_111114 | 0.6667  | 0.8609  | 0.8694  | 0.8489  | 0.8812  | 0.7448    | 0.9196  | 0.7877  | 0.9026    | 0.7731   | 0.8254 
#
#| 10 | -2      | top8000_r45890_111114     | 0.7838  | 0.8710  | 0.7855  | 0.7711  | 0.8647  | 0.6983    | 0.8615  | 0.7660  | 0.7800    | 0.8056   | 0.7813 
#| 11 | -2      | top8000_olf_r45890_111114 | 0.6944  | 0.8276  | 0.7739  | 0.7594  | 0.8284  | 0.7034    | 0.7937  | 0.6944  | 0.6818    | 0.8415   | 0.7667 
#
#| 13 | 3       | top8000_r45890_111114     | 0.5385  | 0.8462  | 0.6482  | 0.6172  | 0.8667  | 0.7983    | 0.7500  | 0.6112  | 0.6839    | 0.6555   | 0.6362 
#| 14 | 3       | top8000_olf_r45890_111114 | 0.7333  | 0.6000  | 0.6417  | 0.6347  | 0.7143  | 0.8182    | 0.7143  | 0.7151  | 0.6963    | 0.6796   | 0.6707 
#
#| 16 | -3      | top8000_r45890_111114     | 0.8750  | 1.0000  | 0.9083  | 0.8601  | 0.3333  | 0.8529    | 0.8333  | 0.8478  | 0.7760    | 0.9098   | 0.8873 
#| 17 | -3      | top8000_olf_r45890_111114 | 0.7727  | 0.6000  | 0.8769  | 0.8472  | 0.3125  | 0.7808    | 1.0000  | 0.7870  | 0.7840    | 0.8573   | 0.8413 
#
#| 19 | 4       | top8000_r45890_111114     | 0.7727  | 0.7288  | 0.7699  | 0.7795  | 0.8140  | 0.5379    | 0.7778  | 0.7562  | 0.5708    | 0.9106   | 0.7329 
#| 20 | 4       | top8000_olf_r45890_111114 | 0.8947  | 0.7347  | 0.7651  | 0.7541  | 0.7750  | 0.7130    | 0.7188  | 0.7453  | 0.6336    | 0.9115   | 0.7426 
#
#| 22 | -4      | top8000_r45890_111114     | 0.7333  | 0.7812  | 0.8864  | 0.8507  | 0.7436  | 0.8794    | 0.9057  | 0.9156  | 0.9115    | 0.9253   | 0.8931 
#| 23 | -4      | top8000_olf_r45890_111114 | 0.6923  | 0.6429  | 0.8386  | 0.7996  | 0.7059  | 0.8271    | 0.7500  | 0.8445  | 0.8392    | 0.8766   | 0.8375 
#
#| 25 | long    | top8000_r45890_111114     | 0.8801  | 0.8584  | 0.8452  | 0.8357  | 0.8450  | 0.8746    | 0.8654  | 0.8302  | 0.8447    | 0.8998   | 0.8559 
#| 26 | long    | top8000_olf_r45890_111114 | 0.8310  | 0.8189  | 0.8347  | 0.8164  | 0.8143  | 0.8491    | 0.8164  | 0.7785  | 0.8190    | 0.8877   | 0.8324 
#
#| 28 | (all)   | (all)                     | 0.8318  | 0.8371  | 0.8116  | 0.8094  | 0.8284  | 0.8360    | 0.8450  | 0.7903  | 0.8068    | 0.8742   | 0.8272 
#|============================================================================================================================================================== 


})) # end FeaturesAnalysis
