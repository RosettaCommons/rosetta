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
id = "hbond_score_vs_score",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("StructureFeatures", "HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <- "
CREATE TEMPORARY TABLE ref_hbonds (
	tag TEXT,
	don_resNum INTEGER,
	don_atmNum INTEGER,
	acc_resNum INTEGER,
	acc_atmNum INTEGER,
	don_chem_type TEXT,
	acc_chem_type TEXT,
	energy REAL,
	PRIMARY KEY(tag, don_resNum, don_atmNum, acc_resNum, acc_atmNum));

INSERT INTO ref_hbonds SELECT
	struct.tag,
	don.resNum AS don_resNum, don.atmNum AS don_atmNum,
	acc.resNum AS acc_resNum, acc.atmNum AS acc_atmNum,
	don.HBChemType AS don_chem_type,
	acc.HBChemType AS acc_chem_type,
	hb.energy AS energy
FROM
	ref.structures AS struct,
	ref.hbonds AS hb,
	ref.hbond_sites AS don, hbond_sites AS acc,
	ref.hbond_sites_pdb AS don_pdb, ref.hbond_sites_pdb AS acc_pdb
WHERE
	hb.struct_id = struct.struct_id AND
	don.struct_id = hb.struct_id AND don.site_id = hb.don_id AND
	acc.struct_id = hb.struct_id AND acc.site_id = hb.acc_id AND
	don_pdb.struct_id = hb.struct_id AND don_pdb.site_id = hb.don_id AND
	don_pdb.heavy_atom_temperature < 30 AND
	acc_pdb.struct_id = hb.struct_id AND acc_pdb.site_id = hb.acc_id AND
	acc_pdb.heavy_atom_temperature < 30;

CREATE UNIQUE INDEX IF NOT EXISTS new.hbond_sites_struct_id_resNum_atmNum ON hbond_sites(
	struct_id, resNum, atmNum);

CREATE UNIQUE INDEX IF NOT EXISTS new.hbonds_struct_id_acc_id_don_id ON hbonds(
	struct_id, acc_id, don_id);

CREATE UNIQUE INDEX IF NOT EXISTS new.structures_tag ON structures( tag );


SELECT
	ref_hb.*,
	new_hb.energy AS new_energy
FROM
	ref_hbonds AS ref_hb,
	new.structures AS new_struct,
	new.hbonds AS new_hb,
	new.hbond_sites AS new_don, new.hbond_sites AS new_acc
WHERE
	new_struct.tag = ref_hb.tag AND
	new_don.struct_id = new_struct.struct_id AND
	new_don.resNum = ref_hb.don_resNum AND
	new_don.atmNum = ref_hb.don_atmNum AND

	new_acc.struct_id = new_struct.struct_id AND
	new_acc.resNum = ref_hb.acc_resNum AND
	new_acc.atmNum = ref_hb.acc_atmNum AND

	new_hb.struct_id = new_struct.struct_id AND
	new_hb.acc_id = new_acc.site_id AND
	new_hb.don_id = new_don.site_id;"

f <- query_sample_sources_against_ref(sample_sources, sele)

f$energy_weighted <- f$energy
f$new_energy_weighted <- f$new_energy

#hbond_sr_bb weighting
f <- transform(f,
	energy_weighted = ifelse(
		don_chem_type == 'hbdon_PBA' &
		acc_chem_type == 'hbacc_PBA' &
		acc_resNum - don_resNum < 5 &
		acc_resNum - don_resNum > -5,
		energy_weighted*.585, energy_weighted))

#hbond_lr_bb weighting
f <- transform(f,
	energy_weighted = ifelse(
		don_chem_type == 'hbdon_PBA' &
		acc_chem_type == 'hbacc_PBA' &
		acc_resNum - don_resNum >= 5 &
		acc_resNum - don_resNum <= -5,
		energy_weighted*1.17, energy_weighted ))

#hbond_bb_sc weighting
f <- transform(f,
	energy_weighted = ifelse(
		(don_chem_type == 'hbdon_PBA' & acc_chem_type != 'hbacc_PBA') |
		(don_chem_type != 'hbdon_PBA' & acc_chem_type == 'hbacc_PBA'),
		energy_weighted*1.17, energy_weighted))

#hbond_sc weighting
f <- transform(f,
	energy_weighted = ifelse(
		don_chem_type != 'hbdon_PBA' & acc_chem_type != 'hbacc_PBA',
		energy_weighted*1.1, energy_weighted))



#hbond_sr_bb hbond_lr_bb hbond_bb_sc weighting
f <- transform(f,
	new_energy_weighted = ifelse(
		don_chem_type == 'hbdon_PBA' | acc_chem_type == 'hbacc_PBA',
		new_energy_weighted*1.17, new_energy_weighted))

#hbond_sc weighting
f <- transform(f,
	new_energy_weighted = ifelse(
		don_chem_type != 'hbdon_PBA' & acc_chem_type != 'hbacc_PBA',
		new_energy_weighted*1.1, new_energy_weighted))



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

diagonal <- data.frame(
	x=seq(-1.755, 0, length.out=200),
	y=seq(-1.755, 0, length.out=200))


plot_id <- "hbond_score_vs_score_rosetta_vs_sp2"
p <- ggplot(data=f, aes(x=energy_weighted, y=new_energy_weighted)) +
	theme_bw() +
	geom_point(size=.4) +
	geom_line(data=diagonal, aes(x=x, y=y), color="darkgray") +
	stat_density2d(size=.2) +
	facet_grid(acc_chem_type ~ don_chem_type) +
	ggtitle(paste("Weighted HBond Energy Rosetta vs Sp2 Potential:") ) +
	scale_x_continuous("Weighted Rosetta HBond Energy", limits=c(-1.755, 0), breaks=c(-1.5, -1, -.5)) +
	scale_y_continuous("Weighted Sp2 Potential HBond Energy", limits=c(-3, 0), breaks=c(-3, -2, -1))
if(nrow(sample_sources) <= 3){
	p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
