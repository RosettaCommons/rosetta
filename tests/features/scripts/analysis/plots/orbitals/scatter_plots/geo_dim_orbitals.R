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
id = "geo_dim_orbitals",
author = "Matthew O'Meara, Steven Combs",
brief_description = "",
feature_reporter_dependencies = c("OrbitalFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
  OrbHdist, AOH_angle, DHO_angle,
  htype2, orbName1
FROM
  HPOL_orbital
WHERE
  OrbHdist < 4;"


f <- query_sample_sources(sample_sources, sele)


f$orb_type <- factor(f$orbName1,
  levels = c("C.pi.sp2", "N.pi.sp2", "N.p.sp2", "O.pi.sp2", "O.p.sp2",
            "O.p.sp3", "O.pi.sp2.bb", "O.p.sp2.bb"))
f$htype2 <- factor(f$htype2)


plot_parts_scatter <- list(
  theme_bw(),
  aes(y=OrbHdist),
  geom_point(size=.1),
  stat_density2d(size=.2),
  facet_grid(orb_type ~ htype2 ),
  labs(y=expression(paste('(Orbital -- Hydrogen Distance) (', ring(A), ')'))),
  scale_x_continuous(limits=c(-1,1)),
  scale_y_continuous(limits=c(0,4)))

###AOH
d_ply(sample_sources, .(sample_source), function(ss){
	ss_id <- ss$sample_source[1]
	plot_id <- paste("geo_dim_AOH_2d_scatter_OrbHpol", ss_id, sep="_")
	ggplot(f, aes(x=AOH_angle, y=OrbHdist)) + plot_parts_scatter +
		labs(x="cos(Acceptor -- Orbital -- Hpol)") +
		ggtitle(paste("Acceptor Orbital Hpol 2D plot\nss_id: ", ss_id, sep=""))
	save_plots(plot_id, ss, output_dir, output_formats)
})

###DHO
d_ply(sample_sources, .(sample_source), function(ss){
  ss_id <- ss$sample_source[1]
  plot_id <- paste("geo_dim_DHO_2d_scatter_OrbHpol", ss_id, sep="_")
  ggplot(f, aes(x=DHO_angle, y=OrbHdist)) + plot_parts_scatter +
    labs(x="cos(Donor -- Hpol -- Orbital)") +
    ggtitle(paste("Donor Hpol Orbital 2D plot\nss_id: ", ss_id, sep=""))
  save_plots(plot_id, ss, output_dir, output_formats)
})


######################HARO
sele <-"
SELECT
  OrbHdist, AOH_angle, DHO_angle,
  htype2, orbName1
FROM
  HARO_orbital
WHERE
  OrbHdist < 5;"
f <- query_sample_sources(sample_sources, sele)

f$orb_type <- factor(f$orbName1,
  levels = c("C.pi.sp2", "N.pi.sp2", "N.p.sp2", "O.pi.sp2", "O.p.sp2",
            "O.p.sp3", "O.pi.sp2.bb", "O.p.sp2.bb"))
f$htype2 <- factor(f$htype2)

###AOH
d_ply(sample_sources, .(sample_source), function(ss){
  ss_id <- ss$sample_source[1]
  plot_id <- paste("geo_dim_AOH_2d_scatter_OrbHaro", ss_id, sep="_")
  ggplot(f, aes(x=AOH_angle, y=OrbHdist)) + plot_parts_scatter +
    labs(x="cos(Acceptor -- Orbital -- Haro)") +
    ggtitle(paste("Acceptor Orbital Haro 2D plot\nss_id: ", ss_id, sep=""))
  save_plots(plot_id, ss, output_dir, output_formats)
})

###DHO
d_ply(sample_sources, .(sample_source), function(ss){
  ss_id <- ss$sample_source[1]
  plot_id <- paste("geo_dim_DHO_2d_scatter_OrbHaro", ss_id, sep="_")
  ggplot(f, aes(x=DHO_angle, y=OrbHdist)) + plot_parts_scatter +
    labs(x="cos(Donor -- Haro -- Orbital)") +
    ggtitle(paste("Donor Haro Orbital 2D plot\nss_id: ", ss_id, sep=""))
  save_plots(plot_id, ss, output_dir, output_formats)
})
})) # end FeaturesAnalysis

