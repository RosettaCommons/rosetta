# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


radial_3d_normalization <- function(x){ 1/(x^2*sum(1/x^2))}


uniform_normalization <- function(x){ rep(1/length(x),length(x)) }


no_normalization <- function(x) rep(1,length(x))


# Equal Area Coordinate Grids
major_long_coords <- transform(
	expand.grid(long=seq(0, 2*pi, length.out=200), lat=c(pi/6, pi/2)),
	capx = 2*sin(lat/2)*cos(long),
	capy = 2*sin(lat/2)*sin(long))

minor_long_coords <- transform(
	expand.grid(long=seq(0, 2*pi, length.out=200), lat=c(pi/3, pi*2/3)),
	capx = 2*sin(lat/2)*cos(long),
	capy = 2*sin(lat/2)*sin(long))

major_lat_coords <- transform(
	expand.grid(long=seq(0,2*pi, length.out=5), lat=c(pi/3, pi*2/3, length.out=200)),
	capx = 2*sin(lat/2)*cos(long),
	capy = 2*sin(lat/2)*sin(long))

minor_lat_coords <- transform(
	expand.grid(long=seq(0,2*pi, length.out=5) + pi/4, lat=c(pi/3, pi*2/3, length.out=200)),
	capx = 2*sin(lat/2)*cos(long),
	capy = 2*sin(lat/2)*sin(long))


polar_equal_area_grids_bw <- function(scale=1, ...) {
	list(
		geom_path(
			data=minor_long_coords,
			aes(x=capx, y=capy, group=lat),
			size=scale * .5, colour="grey98", ...),
		geom_path(
			data=minor_lat_coords,
			aes(x=capx, y=capy, group=long),
			size=scale * .5, colour="grey98", ...),
		geom_path(
			data=major_long_coords,
			aes(x=capx, y=capy, group=lat),
			size=scale * .2, colour="grey90", ...),
		geom_path(
			data=major_lat_coords,
			aes(x=capx, y=capy, group=long),
			size=scale * .2, colour="grey90", ...),
		opts(
			panel.grid.major = theme_blank(),
			panel.grid.minor = theme_blank()))
}
