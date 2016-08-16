# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


require(ggplot2)
require(proto)

GeomIndicator <- proto(ggplot2:::Geom, {
	objname <- "indicator"
	draw <- function(., data, scales, coordinates, ...){

		if("xpos" %in% names(data)){
			if(data$xpos[1] == "left"){
				xpos <- .07
			} else if( data$xpos[1] == "center"){
				xpos <- .5
			} else if( data$xpos[1] == "right"){
				xpos <- .97
			} else if( is.numeric(data$xpos[1]) && 0 <= data$xpos[1] && data$xpos[1] <= 100){
				xpos <- data$xpos[1]
			} else {
				stop(paste("In geom_indicator(), unrecognized value xpos=\"", data$xpos[1],"\".	Please use 'left', 'right' or 'center', or a value from 0 to 1.", sep=""))
			}
		} else {
			xpos <- .97
		}

		if("ypos" %in% names(data)){
			if(data$ypos[1] == "top"){
				ypos <- .97
			} else if( data$ypos[1] == "center"){
				ypos <- .5
			} else if( data$ypos[1] == "bottom"){
				ypos <- .03
			} else if( is.numeric(data$ypos[1]) && 0 <= data$ypos[1] && data$ypos[1] <= 100){
				ypos <- data$ypos[1]
			} else {
				stop(paste("In geom_indicator(), unrecognized value ypos=\"", data$ypos[1],"\".	Please use 'top', 'bottom' or 'center', or a value from 0 to 1.", sep=""))
			}
		} else {
			ypos <- .97
		}

		if(!is.null(data$xjust[1])){
			if(data$xjust[1] %in% c("left", "center", "right")){
				xjust <- data$xjust[1]
			} else {
				stop(paste("In geom_indicator(), unrecognized value xjust=\"", data$xjust[1],"\". Please use 'left', 'right' or 'center'.", sep=""))
			}
		} else {
			 if(xpos < 1/3){
				 xjust <- "left"
			 } else if(xpos >= 1/3 && xpos < 2/3){
				 xjust <- "center"
			 } else {
				 xjust <- "right"
			 }
		}

		if(!is.null(data$yjust[1])){
			if(data$yjust[1] %in% c("top", "center", "bottom")){
				yjust <- data$yjust[1]
			} else {
				stop(paste("In geom_indicator(), unrecognized value yjust=\"", data$yjust[1],"\". Please use 'top', 'bottom' or 'center'.", sep=""))
			}
		} else {
			 if(ypos < 1/3){
				 yjust <- "bottom"
			 } else if(ypos >= 1/3 && ypos < 2/3){
				 yjust <- "center"
			 } else {
				 yjust <- "top"
			 }
		}

		indicator <- data$indicator[1]
		size <- data$size[1]
		if(!is.na(indicator) && !is.null(indicator)){
				if(is.character(indicator)){
		indicator_display_value <- indicator
				} else {
					indicator_display_value <- prettyNum(data$indicator[1], big.mark=",")
				}
			level <- data$group[1] - 1
				colour <- data$colour[1]
			textGrob(indicator_display_value,
						 unit(xpos, "npc"), unit(ypos, "npc") - unit(level, "line"),
						 just=c(xjust, yjust),
						 gp=gpar(col=colour, fontsize=size*12/5, cex=.75))
		}
	}
	desc_params <- list()
	default_stat <- function(.) StatIdentity
	icon <- function(.) textGrob("text", gp=gpar(cex=.75))
	desc <- "Count Instances"
	required_aes <- c("indicator")
	default_aes <- function(.) aes(colour="black", xpos="right", ypos="top", xjust=NULL, yjust=NULL, size=5, group=1)
	guide_geom <- function(x) "blank"
	example <- function(.){
		data <- rbind(
			data.frame(x= c(rnorm(200)+5,rnorm(170)+2),strata=factor(1),slice =factor(1)),
			data.frame(x= c(rnorm(500)+5,rnorm(32 )+2),strata=factor(2),slice =factor(1)),
			data.frame(x= c(rnorm(356)+5,rnorm(120)+2),strata=factor(2),slice =factor(2)))
		p <- ggplot(data=data, aes(x=x, colour=strata))
		p <- p + geom_density()
		p <- p + geom_indicator()
		p <- p + facet_wrap( ~ slice )
		print(p)
	}
})


if(compareVersion(sessionInfo()$otherPkgs$ggplot2$Version, "0.9.0") < 0){
	geom_indicator <- GeomIndicator$build_accessor()
} else {
	geom_indicator <- function (
		mapping = NULL,
		data = NULL,
		stat = "identity",
		position = "identity",
		parse = FALSE,
		...
	) {
		GeomIndicator$new(
			mapping = mapping,
			data = data,
			stat = stat,
			position = position,
			parse = parse,
			...)
	}
}
