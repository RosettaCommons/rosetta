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
id = "rotamer_angles",
filename = "scripts/analysis/plots/backbone_geometry/rotamer_angles.R",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("ResidueFeatures", "ProteinResidueConformationFeatures", "PdbDataFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


sele <-"
CREATE TABLE IF NOT EXISTS nchi_in_res (
        name3 TEXT,
        nchi INTEGER,
        PRIMARY KEY (name3));
INSERT OR IGNORE INTO nchi_in_res VALUES('ARG', 4);
INSERT OR IGNORE INTO nchi_in_res VALUES('LYS', 4);
INSERT OR IGNORE INTO nchi_in_res VALUES('MET', 3);
INSERT OR IGNORE INTO nchi_in_res VALUES('GLN', 3);
INSERT OR IGNORE INTO nchi_in_res VALUES('GLU', 3);
INSERT OR IGNORE INTO nchi_in_res VALUES('TYR', 2);
INSERT OR IGNORE INTO nchi_in_res VALUES('ILE', 2);
INSERT OR IGNORE INTO nchi_in_res VALUES('ASP', 2);
INSERT OR IGNORE INTO nchi_in_res VALUES('TRP', 2);
INSERT OR IGNORE INTO nchi_in_res VALUES('PHE', 2);
INSERT OR IGNORE INTO nchi_in_res VALUES('HIS', 2);
INSERT OR IGNORE INTO nchi_in_res VALUES('ASN', 2);
INSERT OR IGNORE INTO nchi_in_res VALUES('THR', 1);
INSERT OR IGNORE INTO nchi_in_res VALUES('SER', 1);
INSERT OR IGNORE INTO nchi_in_res VALUES('PRO', 1);
INSERT OR IGNORE INTO nchi_in_res VALUES('CYS', 1);
INSERT OR IGNORE INTO nchi_in_res VALUES('VAL', 1);
INSERT OR IGNORE INTO nchi_in_res VALUES('LEU', 1);
INSERT OR IGNORE INTO nchi_in_res VALUES('ALA', 0);
INSERT OR IGNORE INTO nchi_in_res VALUES('GLY', 0);

SELECT
	res.name3 AS res_type,
	res_dofs.chi1, res_dofs.chi2, res_dofs.chi3, res_dofs.chi4,
	nchi_in_res.nchi
FROM
	residues AS res,
	residue_pdb_confidence AS res_conf,
	residue_secondary_structure AS ss,
	protein_residue_conformation AS res_dofs,
	nchi_in_res
WHERE
	res_conf.struct_id = res.struct_id AND res_conf.residue_number = res.resNum AND
	res_conf.max_temperature < 30 AND
	ss.struct_id = res.struct_id AND ss.resNum == res.resNum AND
	res_dofs.struct_id = res.struct_id AND res_dofs.seqpos == res.resNum AND
	nchi_in_res.name3 = res.name3;"

f <- query_sample_sources(sample_sources, sele)

m_f <- melt(f, measure.vars=c("chi1", "chi2", "chi3", "chi4"), variable_name="chi_angle")

dens <- estimate_density_1d(
	m_f, c("sample_source", "res_type", "chi_angle", "nchi"), "value", xlim=c(-180, 180))

d_ply(dens, .(chi_angle), function(sub_f){
	chi <- sub_f$chi_angle[1]
	print(as.numeric(chi))
	plot_id <- paste("rotamer_geometry", chi, sep="_")
	p <- ggplot(
		data=dens[dens$chi_angle == chi & dens$nchi >= as.numeric(chi),]) +
		theme_bw() +
		geom_line(aes(x=x, y=y, colour=sample_source)) +
		geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
		facet_wrap( ~ res_type) +
		ggtitle(("Rotamer ", chi, ", BFact < 30", sep="")) +
		scale_x_continuous("Dihedral Angle")
	if(nrow(sample_sources) <= 3){
		p <- p + theme(legend.position="bottom", legend.direction="horizontal")
	}

	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
})

})) # end FeaturesAnalysis
