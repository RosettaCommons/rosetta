# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
require(proto)

radial_3d_normalization <- function(x){ 1/(x^2*sum(1/x^2))}

conical_3d_normalization <- function(x){1*abs(cos(x))}


uniform_normalization <- function(x){ rep(1/length(x),length(x)) }


no_normalization <- function(x) rep(1,length(x))


#Adapted from
#http://stackoverflow.com/questions/7660893/boxed-geom-text-with-ggplot2
GeomTextBoxed <- proto(ggplot2:::GeomText, {
	objname <- "text_boxed"
	draw_groups <- function(., ...) .$draw(...)
	draw <- function(.,
		data, scales, coordinates, ..., parse = FALSE,
		expand = 1.2, bgcol = "grey50", bgfill = NA, bgalpha = .8) {
		lab <- data$label
		if (parse) {
			lab <- parse(text = lab)
		}
		with(coord_transform(coordinates, data, scales), {
			sizes <- llply(1:nrow(data),
				function(i) with(data[i, ], {
					grobs <- textGrob(lab[i], default.units="native", rot=angle, gp=gpar(fontsize=size * .pt))
					list(w = grobWidth(grobs), h = grobHeight(grobs))
				}))

			gList(rectGrob(x, y,
				width = do.call(unit.c, lapply(sizes, "[[", "w")) * expand,
				height = do.call(unit.c, lapply(sizes, "[[", "h")) * expand,
				gp = gpar(col = alpha(bgcol, bgalpha), fill = alpha(bgfill, bgalpha))),
				.super$draw(., data, scales, coordinates, ..., parse))
		})
	}

  draw_legend <- function(., data, ...) {
		data <- aesdefaults(data, .$default_aes(), list(...))
		with(data,
			textGrob("a", 0.5, 0.5, rot = angle,
				gp=gpar(col=alpha(colour, alpha), fontsize = size * .pt)))
	}

  icon <- function(.) textGrob("text", rot=45, gp=gpar(cex=1.2))
  default_stat <- function(.) StatIdentity
  required_aes <- c("x", "y", "label")
  default_aes <- function(.)
		aes(colour="black", size=5 , angle=0, hjust=0.5,
			vjust=0.5, alpha = 1, family="", fontface=1, lineheight=1.2)
  guide_geom <- function(x) "text"
})

geom_text_boxed <- function(
	mapping = NULL,
	data = NULL,
	stat = "identity",
	position = "identity",
	parse = FALSE,
	...
) {
	GeomTextBoxed$new(
		mapping = mapping,
		data = data,
		stat = stat,
		position = position,
		parse = parse,
		...)
}


# These the coordinates for the Lambert-Azmuthal plots:
#   longitude (around) 30, 60, 90, 120 degrees
#   latitude  (in-out) 30, 60, 90, 120 degrees (starting from the positive x-axis)

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

long_labels <- transform(
	expand.grid(long=pi*3/2, lat=c(pi/6, pi/3, pi/2, pi*2/3)),
	label = as.character(round(180/pi*lat, 0)),
	capx = 2*sin(lat/2)*cos(long),
	capy = 2*sin(lat/2)*sin(long))

polar_equal_area_grids_bw <- function(scale=1, label_scale=1, bgcolor="#00007F", ...) {
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
		geom_text_boxed(
			data=long_labels,
			aes(x=capx, y=capy, label=paste(label, sep="")),
			parse=T,
			size=label_scale*3, colour="white", bgcol=NA, bgfill=bgcolor),
		opts(
			panel.grid.major = theme_blank(),
			panel.grid.minor = theme_blank()))
}
