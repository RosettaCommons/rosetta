# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

check_setup()
feature_analyses <- c(feature_analyses, new("FeaturesAnalysis",
id = "orbitals_dist",
filename = "scripts/analysis/plots/orbitals_dist.R",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("OrbitalFeatures"),
run=function(self){



sele <-"
SELECT
  dist,
	htype2,
  orbName1
FROM
  HPOL_orbital
WHERE
   dist < 5;"

all_geom <- query_sample_sources(sample_sources, sele)

all_geom$orb_type <- factor(all_geom$orbName1,
	levels = c("C.pi.sp2", "N.pi.sp2", "N.p.sp2", "O.pi.sp2", "O.p.sp2",
            "O.p.sp3", "S.p.sp3", "O.pi.sp2.bb", "O.p.sp2.bb"))
all_geom$htype2 <- factor(all_geom$htype2)

plot_parts <- list(
  theme_bw(),
  geom_line(aes(colour=orbName1)),
  geom_indicator(aes(indicator=counts, colour=orbName1)),
  facet_grid(htype2 ~ orb_type),
  scale_y_continuous("log(FeatureDensity + 1)"))




dens <- estimate_density_1d(
  data = all_geom,
  ids = c("sample_source", "orbName1"),
  variable = "dist",
  weight_fun = radial_3d_normalization)

plot_id <- "HPOL_orbHdist"
#p <- ggplot(data=dens, aes(x=x, y=log(y+1), colour=orbName1, indicator=counts)) + plot_parts
p <- ggplot(data=dens, aes(x=x, y=log(y+1))) + plot_parts
p <- p + geom_line()
p <- p + geom_indicator()
p <- p + opts(title = "HPOL_sc_H_sc_orb Interaction Distance\nnormalized for equal weight per unit distance")
p <- p + labs(x=expression(paste('Acceptor -- Donor Distance (', ring(A), ')')),
              y="log(FeatureDensity + 1)")
p <- p + theme_bw()

save_plots(self, plot_id, sample_sources, output_dir, output_formats)

})) # end FeaturesAnalysis
