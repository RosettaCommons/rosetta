# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


date_code <- function(d=NA){
	# reference http://www.r-cookbook.com/node/17
	if(is.na(d)) d <- Sys.Date()
	pattern <- '20([[:digit:]]{2})-([[:digit:]]{2})-([[:digit:]]{2})'
	paste(
		sub(pattern, '\\1', d), sub(pattern, '\\2', d), sub(pattern, '\\3', d),
		sep="")
}




# Save the last ggplot() object created. For each output format,
# generate a plot and put in the output directory
save_plots <- function(
	features_analysis,
	plot_id,
	sample_sources,
	output_dir,
	output_formats,
	...
) {
	tryCatch(plot_id, error=function(e){
		stop(paste(
			"ERROR: Unable to save the plot because ",
			"the 'plot_id' is not specified.\n", e, sep=""))
	})

	tryCatch(features_analysis, error=function(e){
		stop(paste(
			"ERROR: Unable to save the plot '", plot_id,"' ",
			"because the specified 'features_analysis' is not valid.\n",
			e, sep=""))
	})

	tryCatch(sample_sources, error=function(e){
		stop(paste(
			"ERROR: Unable to save the plot '", plot_id, "' ",
			"because the specified 'sample_sources' is not valid.\n",
			e, sep=""))
	})

	if(nrow(sample_sources)==0){
		stop(paste(
			"ERROR: Unable to save the plot '", plot_id, "' ",
			"because no sample_sources were specified.\n", e, sep=""))
	}

	tryCatch(output_dir, error=function(e){
		stop(paste(
			"ERROR: Unable to save the plot '", plot_id, "' ",
			"because the specified 'output_dir' ",
			"is not a valid variable.\n",
			e, sep=""))
	})

	tryCatch(output_formats, error=function(e){
		stop(paste(
			"ERROR: Unable to save the plot '", plot_id, "' ",
			"because the 'output_formats' parameter is not valid.\n",
			e, sep=""))
	})
	plot_formats <- output_formats[output_formats$type == "plot",]

	if(nrow(plot_formats)==0){
		stop(paste(
			"ERROR: Unable to save the plot '", plot_id, "' ",
			"because no output formats were specified.", sep=""))
	}

	a_ply(plot_formats, 1, function(fmt){
		full_output_dir <- file.path(output_dir, features_analysis@id, fmt$id)
		if(!file.exists(full_output_dir)){
			dir.create(full_output_dir, recursive=TRUE)
		}
		date <- date_code()
		fname <- paste(plot_id, date, sep="_")
		full_path <- file.path(full_output_dir, paste(fname, fmt$extension, sep=""))
		cat("Saving Plot: ", full_path)
		timing <- system.time({
			tryCatch({
				plot <- last_plot()
				if(as.character(fmt$id) == "output_minimal_raster" ||
					 as.character(fmt$id) == "output_minimal_vector" ||
					 as.character(fmt$id) == "output_minimal_pdf"){
					plot$labels$title <- NULL
					plot <- plot + theme(legend.position="none")
				}

				if(fmt$add_footer){
					ggsave_with_footer(
						filename=full_path,
						plot=plot,
						width=fmt$width,
						height=fmt$height,
						dpi=fmt$dpi,
						scale=fmt$scale,
						footer_text=analysis_script,
						...)
				} else {
					ggsave(
						filename=full_path,
						plot=plot,
						width=fmt$width,
						height=fmt$height,
						dpi=fmt$dpi,
						scale=fmt$scale,
						...)
				}
			}, error=function(e){
				cat("\n")
				cat(paste(
					"ERROR: Generating and saving the plot:\n",
					e, sep=""))
			})
		})
		cat(" ... ", as.character(round(timing[3],2)), "s\n", sep="")
	})
}

# Egregious code duplication of ggsave to accomodate adding a footer
ggsave_with_footer <- function (
	filename = default_name(plot),
	plot = last_plot(),
	device = default_device(filename),
	path = NULL,
	scale = 1,
	width = par("din")[1],
	height = par("din")[2], dpi = 300,
	keep = plot$options$keep,
	drop = plot$options$drop,
	footer_text=NULL,
	...)
{
	if (!inherits(plot, "ggplot"))
		stop("plot should be a ggplot2 plot")
	eps <- ps <- function(..., width, height) grDevices::postscript(...,
		width = width, height = height, onefile = FALSE, horizontal = FALSE,
		paper = "special")
	tex <- function(..., width, height) grDevices::pictex(...,
		width = width, height = height)
	pdf <- function(..., version = "1.4") grDevices::pdf(...,
		version = version)
	svg <- function(...) grDevices::svg(...)
	wmf <- function(..., width, height) grDevices::win.metafile(...,
		width = width, height = height)
	png <- function(..., width, height) grDevices::png(..., width = width,
		height = height, res = dpi, units = "in")
	jpg <- jpeg <- function(..., width, height) grDevices::jpeg(...,
		width = width, height = height, res = dpi, units = "in")
	bmp <- function(..., width, height) grDevices::bmp(..., width = width,
		height = height, res = dpi, units = "in")
	tiff <- function(..., width, height) grDevices::tiff(...,
		width = width, height = height, res = dpi, units = "in")
	default_name <- function(plot) {
		paste(digest.ggplot(plot), ".pdf", sep = "")
	}
	default_device <- function(filename) {
		pieces <- strsplit(filename, "\\.")[[1]]
		ext <- tolower(pieces[length(pieces)])
		match.fun(ext)
	}
	if (missing(width) || missing(height)) {
		message("Saving ", prettyNum(width * scale, digits = 3),
			"\" x ", prettyNum(height * scale, digits = 3), "\" image")
	}
	width <- width * scale
	height <- height * scale
	if (!is.null(path)) {
		filename <- file.path(path, filename)
	}
	device(file = filename, width = width, height = height, ...)
	on.exit(capture.output(dev.off()))
	print(plot, keep = keep, drop = drop)
	add_footer_text(footer_text)
	invisible()
}

# add footer text to the lower right corner
add_footer_text <- function(text){

	# Broken by ggplot2 update to version 0.9.0.  FIXME!

#	seekViewport("background")
#	pushViewport(viewport(name="footer", x=.99, y=.01, just=c(1,0), width=.4, height=.04))
#
#	#A box to help debug where the viewport will be
#	#	grid.rect(gp=gpar(col="red"))
#
#	grid.text(as.character(text), x=1, hjust=1, gp=gpar(fontsize=5, col="lightgray"))
#	upViewport(0)
}

#Utility function to add facet_grid, position legend, and save plots
plot_field = function(p, plot_id, grid = NULL){
    
  if (! is.null(grid)){
    p <- p+ facet_grid(facets=grid)
  }
  if(nrow(sample_sources) <= 3){
    p <- p + theme(legend.position="bottom", legend.direction="horizontal")
  }
  save_plots(self, plot_id, sample_sources, output_dir, output_formats)
}

