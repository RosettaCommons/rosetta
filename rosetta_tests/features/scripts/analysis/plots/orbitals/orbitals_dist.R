# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

check_setup()

sele <-"
SELECT
  dist,
  orbName1
FROM
  orbital_polar_hydrogen_interactions
WHERE
  dist < 4;"

all_geom <- query_sample_sources(sample_sources, sele)

dens <- estimate_density_1d(
  data = all_geom,
  ids = c("sample_source", "orbName1"),
  variable = "dist",
  weight_fun = radial_3d_normalization)

plot_id <- "orbitals_polar_hydrogen_dist"
p <- ggplot(data=dens, aes(x=x, y=log(y+1), colour=orbName1, indicator=counts))
p <- p + geom_line()
p <- p + geom_indicator()
p <- p + opts(title = "Orbital Polar Hydrogen Interaction Distance\nnormalized for equal weight per unit distance")
p <- p + labs(x=expression(paste('Acceptor -- Donor Distance (', ring(A), ')')),
              y="log(FeatureDensity + 1)")
p <- p + theme_bw()

save_plots(plot_id, sample_sources, output_dir, output_formats)
