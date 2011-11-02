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

###### HPOL_sc_H_sc_orb ######
sele <-"
SELECT
	AOH_angle,  -- this is really cos(angle)
  orbName1
FROM
	HPOL_orbital;"

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
