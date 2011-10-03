# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


check_setup()

sele <- "
SELECT
  pdb.heavy_atom_temperature AS temperature,
	label.label AS chem_type
FROM
	hbond_sites AS site,
	hbond_chem_types AS label,
	hbond_sites_pdb AS pdb
WHERE
	site.struct_id = pdb.struct_id AND site.site_id = pdb.site_id AND
	site.HBChemType = label.chem_type AND
	pdb.heavy_atom_temperature != 0 AND pdb.heavy_atom_temperature < 100;"

# Execute the SQL query on each sample source.
f <- query_sample_sources(sample_sources, sele)

dens <- estimate_density_1d(f, c("sample_source", "chem_type"), "temperature")

# Generate a lattice of density plots for each chemical type type
plot_id <- "polar_site_temperature_factor_by_chem_type"
p <- ggplot(dens) + theme_bw() +
	geom_line(aes(x, log(y+1), colour=sample_source)) +
	geom_indicator(aes(indicator=counts, colour=sample_source)) +
	facet_wrap(~chem_type) +
	opts(title = "Polar Site Heavy Atom Temperature Factor By Chemical Type") +
	scale_x_continuous("Temperature Factor")
	scale_y_continuous("log(FeatureDensity + 1)")
if(nrow(sample_sources) <= 3){
	p <- p + opts(legend.position="bottom", legend.direction="horizontal")
}

save_plots(plot_id, sample_sources, output_dir, output_formats)
