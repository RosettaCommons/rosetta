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
  residues AS res
WHERE
  res.resNum = rr.resNum AND
  res.struct_id = rr.struct_id;"

all_geom <-  query_sample_sources(sample_sources, sele)

dens <- estimate_density_1d(
  data = all_geom,
  ids = c("sample_source", "res_type"),
  variable = "divergence")

p <- ggplot(data=dens, aes(x=log(x+1), y=log(y+1), color=sample_source, indicator=counts))
p <- p + geom_line()
p <- p + geom_indicator()
p <- p + facet_wrap( ~ res_type )
p <- p + opts(title = "Rotamer Recovery by Score Type")
p <- p + labs(x="log(Automorphic RMSD + 1)",
              y="log(FeatureDensity + 1)")
p <- p + theme_bw()

save_plots(plot_id, sample_sources, output_dir, output_formats)
