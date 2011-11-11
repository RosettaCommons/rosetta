 # -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

check_setup()

plot_id <- "rotamer_recovery_by_res_type"

sele <-"
SELECT
	rr.divergence AS divergence,
	res.name3 AS res_type
FROM
	rotamer_recovery AS rr,
	residues AS res,
	residue_pdb_confidence AS res_confidence
WHERE
	res.resNum = rr.resNum AND
	res.struct_id = rr.struct_id AND
	res_confidence.struct_id = res.struct_id AND
	res_confidence.residue_number = res.resNum AND
	res_confidence.max_temperature < 20;"

f <- query_sample_sources(sample_sources, sele)

dens <- estimate_density_1d(
	data = f,
	ids = c("sample_source", "res_type"),
	variable = "divergence")

ggplot(data=dens) + theme_bw() +
	geom_line(aes(x=log(x+1), y=log(y+1), colour=sample_source)) +
	geom_indicator(aes(indicator=counts, colour=sample_source)) +
	facet_wrap( ~ res_type ) +
	opts(title = "Rotamer Recovery by Residue Type, B-Factor < 20") +
	labs(x="log(Automorphic RMSD + 1)", y="log(FeatureDensity + 1)") +
save_plots(plot_id, sample_sources, output_dir, output_formats)
