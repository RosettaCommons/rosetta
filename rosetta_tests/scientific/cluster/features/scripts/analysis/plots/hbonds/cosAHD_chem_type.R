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
  geom.cosAHD,
  acc_site.HBChemType AS acc_chem_type,
  don_site.HBChemType AS don_chem_type
FROM
  hbond_geom_coords AS geom,
  hbonds AS hbond,
  hbond_sites AS don_site,
  hbond_sites AS acc_site
WHERE
  hbond.struct_id = geom.struct_id AND
  hbond.hbond_id =  geom.hbond_id AND
  hbond.struct_id = don_site.struct_id AND
  hbond.don_id = don_site.site_id AND
  hbond.struct_id = acc_site.struct_id AND
  hbond.acc_id = acc_site.site_id;"

all_geom <- query_sample_sources(sample_sources, sele)

#### ADD NULL MODEL ####
n_pts=10000
null.geom.raw <- data.frame(x=rnorm(n_pts),y=rnorm(n_pts), z=rnorm(n_pts))

null.geom <- data.frame(
  cosAHD=(null.geom.raw$x/sqrt(null.geom.raw$x^2+null.geom.raw$y^2+null.geom.raw$z^2)),
  sample_source="spherical_null")

all_geom_with_null <- ddply(all_geom, .variables=c("acc_chem_type", "don_chem_type"), function(df){
  null.geom$don_chem_type <- df$don_chem_type[1]
  null.geom$acc_chem_type <- df$acc_chem_type[1]
  rbind(df, null.geom)
})
#########################

##### make density histogram ###############
n_pts <- 50
dens <- ddply(all_geom_with_null, .variables=c("sample_source", "acc_chem_type", "don_chem_type"),
function(df){
  breaks <- seq(from=-1, to=1, length.out=n_pts)
  h <- hist(df$cosAHD, breaks=breaks, freq=FALSE, plot=FALSE)
  d <-data.frame(x=h$mids, y=h$density, counts=nrow(df))
  return(d)
})


n_pts <- 50
xlim=range(all_geom_with_null$cosAHD)
dens <- ddply(all_geom_with_null, .variables=c("sample_source", "acc_chem_type", "don_chem_type"),
function(df){
  lgs <- logspline(df$cosAHD, lbound=xlim[1], ubound=xlim[2])
	x <- seq(from=xlim[1], to=xlim[2], length.out=n_pts)
	y <- dlogspline(x, lgs)
  d <-data.frame(x=x, y=y, counts=nrow(df))
  return(d)
})



###########################################



# The 'hbdon_' part of the donor labels doesn't fit so strip them out
dens$don_chem_type <- sub("^hbdon_", '', dens$don_chem_type)

plot_id = "cosAHD_chem_type"

p <- ggplot(data=dens, aes(x=acos(x)*180/pi, y=-log(y), colour=sample_source, indicator=counts))
p <- p + geom_line()
p <- p + geom_indicator()
p <- p + facet_grid( don_chem_type ~ acc_chem_type)
p <- p + opts(title = "Hydrogen Bonds AHD Angle by Chemical Type\n(normalized for equal volume per unit distance)")
p <- p + labs(x=expression(paste('Acceptor -- Hydrogen -- Donor (degrees)')),
              y="-log(FeatureDensity)")
p <- p + theme_bw()
p <- p + opts(axis.text.y=theme_blank())
p <- p + scale_y_continuous(limits=c(-2,2), breaks=c(-2,-1,0,1,2))
p <- p + scale_x_continuous(limits=c(0,90), breaks=c(0,30,60,90), transform="reverse")

plot_id = "cosAHD_chem_type"
save_plots(plot_id, sample_sources, output_dir, output_formats)


##########################################################

#make density histogram
n_pts <- 30
xlim <- range(all_geom$cosAHD)
dens <- ddply(all_geom, .variables=c("acc_chem_type"),
function(df){
  breaks <- seq(from=-1, to=1, length.out=n_pts)
  h <- hist(df$cosAHD, breaks=breaks, plot=FALSE)
  d <-data.frame(x=h$mids, y=h$density, counts=nrow(df))
  return(d)
})

p <- ggplot(data=dens, aes(x=acos(x)*180/pi, y=-log(y), colour=acc_chem_type, indicator=counts))
p <- p + geom_line()
p <- p + geom_indicator()
p <- p + opts(title = "Hydrogen Bonds AHD Angle by Chemical Type\n(normalized for equal volume per unit distance)")
p <- p + labs(x=expression(paste('Acceptor -- Hydrogen -- Donor (degrees)')),
              y="-log(FeatureDensity)")
p <- p + theme_bw()
p <- p + scale_y_continuous(limits=c(-2,6))
p <- p + scale_x_continuous(limits=c(0,90), breaks=c(0,30,60,90), transform="reverse")

plot_id = "cosAHD_by_acc_chem_type"
save_plots(plot_id, sample_sources, output_dir, output_formats)

##########################################################




#make density histogram
n_pts <- 30

null.geom_with_donor <- null.geom
null.geom_with_donor$don_chem_type <- "spherical_null"
null.geom_with_donor$acc_chem_type <- "spherical_null"

all_geom_with_null <- rbind(all_geom, null.geom_with_donor)
xlim <- range(all_geom$cosAHD)
dens <- ddply(all_geom_with_null, .variables=c("don_chem_type"),
function(df){
  breaks <- seq(from=-1, to=1, length.out=n_pts)
	print(summary(df$cosAHD))
  h <- hist(df$cosAHD, breaks=breaks, plot=FALSE)
  d <-data.frame(x=h$mids, y=h$density, counts=nrow(df))
  return(d)
})

p <- ggplot(data=dens, aes(x=acos(x)*180/pi, y=-log(y), colour=don_chem_type, indicator=counts))
p <- p + geom_line()
p <- p + geom_indicator()
p <- p + opts(title = "Hydrogen Bonds AHD Angle by Chemical Type\n(normalized for equal volume per unit distance)")
p <- p + labs(x=expression(paste('Acceptor -- Hydrogen -- Donor (degrees)')),
              y="-log(FeatureDensity)")
p <- p + theme_bw()
p <- p + scale_x_continuous(limits=c(-2,6))
p <- p + scale_x_continuous(limits=c(0,90), breaks=c(0,30,60,90), transform="reverse")

plot_id = "cosAHD_chem_type_by_don_chem_type"
save_plots(plot_id, sample_sources, output_dir, output_formats)
