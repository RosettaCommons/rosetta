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
id = "anchor_transform_correlation",
filename = "scripts/analysis/plots/loops/anchor_transform_correlation.R",
author = "Brian D. Weitzner",
brief_description = "",
long_description = "
This features analysis computes the correlation between the six dimensions of a 
coordinate frame transformation pairwise
",

feature_reporter_dependencies = c("loop_anchor_features"),
run=function(self, sample_sources, output_dir, output_formats){

library(grid)
  
# maximum number of rows to select 
limit <- 10^7 / 2

# minimum number of hits to make it worthwhile to plot
min_number_of_examples <- 10

number_ticks <- function(n) {function(limits){
  x <- diff(range(limits)) / 20
  seq(limits[1] + x, limits[2] - x, length=n)}}

label_ticks <- function(n) {function(limits){
  x <- diff(range(limits)) / 20
  signif(seq(from=(limits[1] + x), to=(limits[2] - x), length=n), digits=2)}}

plotmatrix <- function (input_data, mapping = aes(), 
                        columns = names(input_data), colour = "black") 
{
  data <- input_data[,columns]
  grid <- expand.grid(x = 1:ncol(data), y = 1:ncol(data))
  grid <- subset(grid, x != y)
  all <- do.call("rbind", lapply(1:nrow(grid), function(i) {
    xcol <- grid[i, "x"]
    ycol <- grid[i, "y"]
    # NOTE: xvar and yvar refer to the position of the panel in the grid
    # while x and y refer to the data being plotted in that panel.
    # Becase of this, xvar and x refer to different variables.
    data.frame(xvar = names(data)[ycol], yvar = names(data)[xcol], 
               x = data[, xcol], y = data[, ycol], input_data)
  }))
  all$xvar <- factor(all$xvar, levels = names(data))
  all$yvar <- factor(all$yvar, levels = names(data))

  densities <- do.call("rbind", lapply(1:ncol(data), function(i) {
    data.frame(xvar = names(data)[i], yvar = names(data)[i], 
               x = data[, i], 
               subset(input_data, select = !(names(input_data) %in% columns)))
  }))

  faceted_density <- mapply(function(d){
    stat_density(aes(x = x, 
                     y = ..scaled.. * diff(range(x))+ min(x)), 
                 data = d, 
                 position = "identity", geom = "line")
  }, dlply(densities, .(xvar, yvar)))
  
  mapping <- plyr::defaults(mapping, aes_string(x = "x", y = "y"))
  class(mapping) <- "uneval"

  p <- ggplot(all, mapping) + facet_grid(xvar ~ yvar, scales="free", 
                                         labeller=label_parsed) +
    scale_color_manual(values=c("darkgrey", "black")) +
    scale_fill_manual(values=c("darkgrey", "black")) +
    stat_bin2d(aes(alpha=..density.., fill=sample_source), show_guide = FALSE) + 
    faceted_density +
    scale_y_continuous("Degree of freedom", breaks=number_ticks(4), 
                       labels=label_ticks(4), expand=c(0.1,0)) +
    scale_x_continuous("Degree of freedom", breaks=number_ticks(4), 
                       labels=label_ticks(4), expand=c(0.1,0)) +
    coord_equal(ratio=1) +
    ggtitle(paste("Scatterplot matrix for loop anchor transform", 
                  "degrees of freedom for", input_data$length, 
                  "residue loops.", sep = " ")) +
    theme_bw() +
    theme(strip.background=element_blank(), 
          strip.text.y=element_text(angle=0),
          axis.text.x=element_text(size=6, family="serif"),
          axis.text.y=element_text(size=6, family="serif"),
          plot.title=element_text(size=10, family="serif"),
          axis.title.x=element_text(size=9, family="serif"),
          axis.title.y=element_text(size=9, family="serif"),
          legend.title=element_blank(),
          legend.text=element_text(size=9, family="serif"),
          legend.key.size=unit(0.025, units="npc"),
          panel.grid.minor=element_blank(), 
          panel.grid.major=element_blank(),
          panel.border = element_rect(colour = "black"),
          panel.margin = unit(0, units="npc"))
    
    
  if(nrow(sample_sources) <= 3){
    p <- p + theme(legend.position="bottom", legend.direction="horizontal")
  }
}

sele <- paste("SELECT residue_begin, residue_end, x, y, z, phi, psi, theta",
              "FROM loop_anchor_transforms LIMIT", limit)

f <- query_sample_sources(sample_sources, sele)

# Add length information to data frame
f$length <- as.factor(f$residue_end - f$residue_begin + 1)

# Scale phi to be on the domain [0, 2pi)
f$phi <- ifelse( f$phi < 0, f$phi + 2*pi, f$phi)

# Split data frame into sub-frames by length of loop
# Plot the results for each loop length and save the plots to disk
d_ply(f, .(length), function(df){
  plot_and_save_to_disk <- TRUE
  for (sample_source in sample_sources$sample_source) {
    if (nrow(df[df$sample_source == sample_source,]) < min_number_of_examples){
      plot_and_save_to_disk <- FALSE
    }
  }
  
  if (plot_and_save_to_disk) {
    # data columns are all columns that do not contain position, length or
    # sample source information
    data_columns <- !(names(df) %in% c("residue_begin","residue_end","length", 
                                       "sample_source"))
    p <- plotmatrix(df, aes(colour=sample_source), subset(names(df), 
                                                          data_columns))
    
    # use cairo to preserve transparency in EPS format
    save_plots(self, paste("anchor_transform_scatterplot_matrix_for", 
                           df$length[1], "residue_loops", sep = "_"), 
               sample_sources, output_dir, 
               output_formats[output_formats$extension == ".eps",], 
               device=cairo_ps)
    
    # don't use cairo for non-EPS
    save_plots(self, paste("anchor_transform_scatterplot_matrix_for", 
                           df$length[1], "residue_loops", sep = "_"), 
               sample_sources, output_dir, 
               output_formats[output_formats$extension != ".eps",])
    
  }
})

})) # end FeaturesAnalysis
