# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

check_setup()

require(logspline)

###### POLAR HYDROGENS ######
sele <-"
SELECT
	angle,  -- this is really cos(angle)
  orbName1
FROM
	orbital_polar_hydrogen_interactions;"

all_geom <- query_sample_sources(sample_sources, sele)

n_pts <- 50
xlim=range(all_geom$angle)
dens <- ddply(all_geom, .variables=c("sample_source", "orbName1"),
function(df){
  lgs <- logspline(df$angle, lbound=xlim[1], ubound=xlim[2])
	x <- seq(from=xlim[1], to=xlim[2], length.out=n_pts)
	y <- dlogspline(x, lgs)
  d <-data.frame(x=x, y=y, counts=nrow(df))
  return(d)
})
###############################################


plot_id <- "orbitals_polar_hydogens_cos_angle"

p <- ggplot(data=dens, aes(x=x, y=log(y+1), colour=orbName1, indicator=counts))
p <- p + geom_line()
p <- p + geom_indicator()
p <- p + opts(title = "Orbital Polar Hydrogen Interaction Cos( Angle )\nnormalized for equal weight per unit distance")
p <- p + labs(x=expression(paste('Cos( Base--Orbital--Donor )')),
              y="log(FeatureDensity + 1)")
p <- p + theme_bw()

save_plots(plot_id, sample_sources, output_dir, output_formats)


###### AROMATIC HYDROGENS ######
sele <-"
SELECT
	angle,  -- this is really cos(angle)
  orbName1
FROM
	orbital_aromatic_hydrogen_interactions;"

all_geom <- query_sample_sources(sample_sources, sele)

n_pts <- 50
xlim=range(all_geom$angle)
dens <- ddply(all_geom, .variables=c("sample_source", "orbName1"),
function(df){
  lgs <- logspline(df$angle, lbound=xlim[1], ubound=xlim[2])
	x <- seq(from=xlim[1], to=xlim[2], length.out=n_pts)
	y <- dlogspline(x, lgs)
  d <-data.frame(x=x, y=y, counts=nrow(df))
  return(d)
})
###############################################


plot_id <- "orbitals_aromatic_hydogens_cos_angle"

p <- ggplot(data=dens, aes(x=x, y=log(y+1), colour=orbName1, indicator=counts))
p <- p + geom_line()
p <- p + geom_indicator()
p <- p + opts(title = "Orbital Aromatic Hydrogen Interaction Cos( Angle )\nnormalized for equal weight per unit distance")
p <- p + labs(x=expression(paste('Cos( Base--Orbital--Donor )')),
              y="log(FeatureDensity + 1)")
p <- p + theme_bw()

save_plots(plot_id, sample_sources, output_dir, output_formats)

###### ORBITALS ######


sele <-"
SELECT
  angle,  -- this is really cos(angle)
  orbName1
FROM
  orbital_orbital_interactions;"

all_geom <- query_sample_sources(sample_sources, sele)

n_pts <- 50
xlim=range(all_geom$angle)
dens <- ddply(all_geom, .variables=c("sample_source", "orbName1"),
function(df){
  lgs <- logspline(df$angle, lbound=xlim[1], ubound=xlim[2])
  x <- seq(from=xlim[1], to=xlim[2], length.out=n_pts)
  y <- dlogspline(x, lgs)
  d <-data.frame(x=x, y=y, counts=nrow(df))
  return(d)
})
###############################################


plot_id <- "orbitals_orbitals_cos_angle"

p <- ggplot(data=dens, aes(x=x, y=log(y+1), colour=orbName1, indicator=counts))
p <- p + geom_line()
p <- p + geom_indicator()
p <- p + opts(title = "Orbital Orbital Interaction Cos( Angle )\nnormalized for equal weight per unit distance")
p <- p + labs(x=expression(paste('Cos( Base--Orbital--Donor )')),
              y="log(FeatureDensity + 1)")
p <- p + theme_bw()

save_plots(plot_id, sample_sources, output_dir, output_formats)
