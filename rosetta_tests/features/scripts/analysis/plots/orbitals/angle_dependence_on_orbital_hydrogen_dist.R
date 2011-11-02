# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

check_setup()


#########################   POLAR   #############################################################

sele <-"
SELECT
	dist, AOH_angle, DHO_angle,
	htype2, orbName1
FROM
	HPOL_orbital
WHERE
	dist < 5;"


f <- query_sample_sources(sample_sources, sele)
f$OrbHdist_quantile <- cut(f$dist, quantile(f$dist, seq(0,1,by=.5)))

f$orb_type <- factor(f$orbName1,
	levels = c("C.pi.sp2", "N.pi.sp2", "N.p.sp2", "O.pi.sp2", "O.p.sp2",
						"O.p.sp3", "S.p.sp3", "O.pi.sp2.bb", "O.p.sp2.bb"))
f$htype2 <- factor(f$htype2)

plot_parts <- list(
	theme_bw(),
	geom_line(aes(colour=OrbHdist_quantile)),
	geom_indicator(aes(indicator=counts, colour=OrbHdist_quantile)),
	facet_grid(htype2 ~ orb_type),
	scale_y_continuous("log(FeatureDensity + 1)"))

n_pts <- 500
xlim=range(f$AOH_angle)


d_ply(f, .variables=("sample_source"), function(sub_f){

function(df){
  lgs <- logspline(df$AOH_angle, lbound=-1.0, ubound=1.0)
  x <- seq(from=xlim[1], to=xlim[2], length.out=n_pts)
  y <- dlogspline(x, lgs)
  d <-data.frame(x=x, y=y, counts=nrow(df))
  return(d)

}
	ss_id <- sub_f[1, "sample_source"]
	
	plot_id <- paste("HPOL_cosAOrbH_OrbHdist_quantiles_", ss_id, "_by_orbital_type",sep="")
	dens <- estimate_density_1d(sub_f,
		c("orb_type", "htype2", "OrbHdist_quantile"), "AOH_angle")
	ggplot(data=dens, aes(x=x, y=log(y+1))) + plot_parts +
		scale_x_continuous('Aceptor -- Orbital -- Hydrogen (cos(angle))') +
		opts(title = paste("HPOL Cos(aceptor, orbital, hydrogen) and OrbHdist quantile\nss_id: ", ss_id, sep=""))
	save_plots(plot_id, sample_sources, output_dir, output_formats)


})



#########################   POLAR   #############################################################
########################Donor - Hydrogen - Orbital##############################################
###############################################################################################
xlim=range(f$DHO_angle)
d_ply(f, .variables=("sample_source"), function(sub_f){
function(df){
  lgs <- logspline(df$DHO_angle, lbound=-1.0, ubound=1.0)
  x <- seq(from=xlim[1], to=xlim[2], length.out=n_pts)
  y <- dlogspline(x, lgs)
  d <-data.frame(x=x, y=y, counts=nrow(df))
  return(d)

}
  ss_id <- sub_f[1, "sample_source"]

  plot_id <- paste("cosDHorb_OrbHdist_quantiles_", ss_id, "_by_orbital_type",sep="")
  dens <- estimate_density_1d(sub_f,
    c("orb_type", "htype2", "OrbHdist_quantile"), "DHO_angle")
  ggplot(data=dens, aes(x=x, y=log(y+1))) + plot_parts +
    scale_x_continuous('Donor -- Hydrogen -- Orbital (cos(angle))') +
    opts(title = paste("Cos(donor, hydrogen, orbital) and OrbHdist quantile\nss_id: ", ss_id, sep=""))
  save_plots(plot_id, sample_sources, output_dir, output_formats)


})

###########################   AROMATIC     #########################################

sele <-"
SELECT
  dist, AOH_angle, DHO_angle,
  htype2, orbName1
FROM
  HARO_orbital
WHERE
  dist < 5;"

f <- query_sample_sources(sample_sources, sele)
f$OrbHdist_quantile <- cut(f$dist, quantile(f$dist, seq(0,1,by=.5)))

f$orb_type <- factor(f$orbName1,
  levels = c("C.pi.sp2", "N.pi.sp2", "N.p.sp2", "O.pi.sp2", "O.p.sp2",
            "O.p.sp3", "S.p.sp3", "O.pi.sp2.bb", "O.p.sp2.bb"))
f$htype2 <- factor(f$htype2)

plot_parts <- list(
  theme_bw(),
  geom_line(aes(colour=OrbHdist_quantile)),
  geom_indicator(aes(indicator=counts, colour=OrbHdist_quantile)),
  facet_grid(htype2 ~ orb_type),
  scale_y_continuous("log(FeatureDensity + 1)"))

n_pts <- 500
xlim=range(f$AOH_angle)


d_ply(f, .variables=("sample_source"), function(sub_f){

function(df){
  lgs <- logspline(df$AOH_angle, lbound=-1.0, ubound=1.0)
  x <- seq(from=xlim[1], to=xlim[2], length.out=n_pts)
  y <- dlogspline(x, lgs)
  d <-data.frame(x=x, y=y, counts=nrow(df))
  return(d)

}
  ss_id <- sub_f[1, "sample_source"]

  plot_id <- paste("HARO_cosAOrbH_OrbHdist_quantiles_", ss_id, "_by_orbital_type",sep="")
  dens <- estimate_density_1d(sub_f,
    c("orb_type", "htype2", "OrbHdist_quantile"), "AOH_angle")
  ggplot(data=dens, aes(x=x, y=log(y+1))) + plot_parts +
    scale_x_continuous('Aceptor -- Orbital -- Hydrogen (cos(angle))') +
    opts(title = paste("HARO Cos(aceptor, orbital, hydrogen) and OrbHdist quantile\nss_id: ", ss_id, sep=""))
  save_plots(plot_id, sample_sources, output_dir, output_formats)


})

#########################   HARO   #############################################################
########################Donor - Hydrogen - Orbital##############################################
###############################################################################################
xlim=range(f$DHO_angle)
d_ply(f, .variables=("sample_source"), function(sub_f){
function(df){
  lgs <- logspline(df$DHO_angle, lbound=-1.0, ubound=1.0)
  x <- seq(from=xlim[1], to=xlim[2], length.out=n_pts)
  y <- dlogspline(x, lgs)
  d <-data.frame(x=x, y=y, counts=nrow(df))
  return(d)

}
  ss_id <- sub_f[1, "sample_source"]

  plot_id <- paste("HARO_cosDHorb_OrbHdist_quantiles_", ss_id, "_by_orbital_type",sep="")
  dens <- estimate_density_1d(sub_f,
    c("orb_type", "htype2", "OrbHdist_quantile"), "DHO_angle")
  ggplot(data=dens, aes(x=x, y=log(y+1))) + plot_parts +
    scale_x_continuous('Donor -- Hydrogen -- Orbital (cos(angle))') +
    opts(title = paste("HARO Cos(donor, hydrogen, orbital) and OrbHdist quantile\nss_id: ", ss_id, sep=""))
  save_plots(plot_id, sample_sources, output_dir, output_formats)


})



