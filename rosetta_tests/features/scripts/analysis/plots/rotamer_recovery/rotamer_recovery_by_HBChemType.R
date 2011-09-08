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
  rr.divergence AS divergence,
  hb_site.HBChemType AS hb_chem_type
FROM
  rotamer_recovery AS rr,
  hbond_sites AS hb_site
WHERE
  hb_site.resNum = rr.resNum AND
  hb_site.struct_id = rr.struct_id;"

all_geom <- query_sample_sources(sample_sources, sele)

dens <- estimate_density_1d(
  data = all_geom,
  ids = c("sample_source", "hb_chem_type"),
  variable = "divergence")

plot_id <- "rotamer_recovery_by_HBChemType"
p <- ggplot(data=dens, aes(x=log(x+1), y=log(y+1), color=sample_source, indicator=counts))
p <- p + geom_line() + geom_indicator()
p <- p + facet_wrap( ~ hb_chem_type )
p <- p + opts(title = "Rotamer Recovery by Hydrogen Bond Chemical Type")
p <- p + labs(x="<- better      log(Automorphic RMSD + 1)      worse ->",
              y="log(FeatureDensity + 1)")
p <- p + theme_bw()
save_plots(plot_id, sample_sources, output_dir, output_formats)
