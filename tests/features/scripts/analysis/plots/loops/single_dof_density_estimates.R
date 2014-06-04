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
id = "single_dof_density_estimates",
filename = "scripts/analysis/plots/loops/single_dof_density_estimates.R",
author = "Brian D. Weitzner",
brief_description = "",
long_description = "
This features analysis computes the density estimate of each degree of freedom
in a three-dimensional transformation at each segment length.
",

feature_reporter_dependencies = c("loop_anchor_features"),
run=function(self, sample_sources, output_dir, output_formats){

library(grid)
  
# maximum number of rows to select 
limit <- 10^7 / 2

# minimum number of hits to make it worthwhile to plot
min_number_of_examples <- 10

plotmatrix <- function (input_data, mapping = aes(), 
                        columns = names(input_data), colour = "black") 
{
  # reformat the data frame to have columns for 
  # 1) the degree of freedom of the value (dof)
  # 2) the value of the observation (x)
  # 3) the dataset the observation came from (sample_source)
  densities <- do.call("rbind", lapply(columns, function(i) {  
    data.frame(dof=rep(i, nrow(input_data)),
               x=input_data[[i]],
               sample_source=input_data$sample_source)
  }))
  
  # plot the densities across 3 columns
  p <- ggplot(densities, mapping) + facet_wrap(~ dof, ncol=3, scales="free") + 
    stat_density(aes(x=x, y = ..scaled..,  fill=sample_source), 
                 position = "identity", alpha=0.2) +
    scale_x_continuous("Degree of Freddom", expand = c(0, 0)) + 
    scale_y_continuous("", expand = c(0.05, -0.04)) +
    scale_color_manual(values=c("darkgrey", "black")) +
    scale_fill_manual(values=c("darkgrey", "black")) +
    ggtitle(paste("Density estimates for each loop anchor transform", 
                  "degree of freedom for", input_data$length, 
                  "residue loops.", sep = " ")) +
    theme_classic() +
    theme(strip.background=element_blank(), 
          axis.text.x=element_text(size=6, family="serif"),
          axis.text.y=element_text(size=6, family="serif"),
          plot.title=element_text(size=10, family="serif"),
          axis.title.x=element_text(size=9, family="serif"),
          axis.title.y=element_text(size=9, family="serif"),
          legend.title=element_blank(),
          legend.key.size=unit(0.025, units="npc"),
          legend.text=element_text(size=9, family="serif"),
          axis.line = element_line(size=0.1,),
          axis.ticks = element_line(size=0.1))
          
  
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
    save_plots(self, paste("anchor_transform_dof_density_estimates_for", 
                           df$length[1], "residue_loops", sep = "_"), 
               sample_sources, output_dir, 
               output_formats[output_formats$extension == ".eps",], 
               device=cairo_ps)
    
    # don't use cairo for non-EPS
    save_plots(self, paste("anchor_transform_dof_density_estimates_for", 
                           df$length[1], "residue_loops", sep = "_"), 
               sample_sources, output_dir, 
               output_formats[output_formats$extension != ".eps",])
  }
})

})) # end FeaturesAnalysis
