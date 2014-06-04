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
id = "tau_qq_plot",
filename = "scripts/analysis/plots/loops/tau_qq_plot.R",
author = "Brian D. Weitzner",
brief_description = "",
long_description = "
This features analysis computes the density estimate of each degree of freedom
in a three-dimensional transformation at each segment length.
",

feature_reporter_dependencies = c("loop_anchor_features"),
run=function(self, sample_sources, output_dir, output_formats){
# maximum number of rows to select 
limit <- 10^7

# minimum number of hits to make it worthwhile to plot
min_number_of_examples <- 10

# calculate parametrs for reference line in QQ plot
# taken and modified from:
# http://stackoverflow.com/questions/4357031/qqnorm-and-qqline-in-ggplot2
qq_params <- function(df, data_col) {
  
  y <- quantile(df[[data_col]][!is.na(df[[data_col]])], c(0.25, 0.75))
  x <- qnorm(c(0.25, 0.75))
  slope <- diff(y)/diff(x)
  int <- y[1L] - slope * x[1L]
  
  return(c(slope, int))
}

plot_and_save <- function(df, title, filename, limit, ratio=0, 
                          data_col="tau101") {
  abline_params <- qq_params(df, data_col)
  
  mapping <- aes_string(sample=data_col)
  
  p <- ggplot(df, mapping) + ggtitle(title) +
    geom_abline(slope=abline_params[1],
                intercept=abline_params[2], 
                size=1.2, alpha=0.9) +
    stat_qq(aes(alpha=0.5)) +
    scale_x_continuous("Theoretical quantiles", limit=c(-3, 3)) +
    scale_y_continuous(expression(paste(tau, " quantiles")), limit=limit) +
    coord_fixed(ratio=ifelse(ratio==0, 1/abline_params[1], ratio))
  
  if(nrow(sample_sources) <= 3){
    p <- p + theme(legend.position="bottom", legend.direction="horizontal")
  }
  
  # use cairo to preserve transparency in EPS format
  save_plots(self, paste(filename, sep = "_"), 
             sample_sources, output_dir, 
             output_formats[output_formats$extension == ".eps",], 
             device=cairo_ps)
  
  # don't use cairo for non-EPS
  save_plots(self, paste(filename, sep = "_"), 
             sample_sources, output_dir, 
             output_formats[output_formats$extension != ".eps",])
}
  
sele <- paste("SELECT residue_begin, residue_end, tau101",
              "FROM loop_anchor_transforms LIMIT", limit)

f <- query_sample_sources(sample_sources, sele)

plot_and_save_to_disk <- TRUE
for (sample_source in sample_sources$sample_source) {
  if (nrow(f[f$sample_source == sample_source,]) < min_number_of_examples){
    plot_and_save_to_disk <- FALSE
  }
}
  
if (plot_and_save_to_disk) {
  ratio <- 1 / qq_params(f, "tau101")[1]
  plot_and_save(f, expression(paste("All values of ", tau)), "tau_qq_plot",
                range(f$tau101))
  
  tau.cutoff = f$tau101[order(f$tau101)[round(length(f$tau101) * 0.8)]]
  # only plot the lower 80% of tau values, but use range of the unmodified 
  # data.frame to put both plots on the same plot area.
  plot_and_save(subset(f, tau101 <= tau.cutoff),
                expression(paste("Lower 80% of ", tau)), "lower_tau_qq_plot",
                range(f$tau101), ratio)
}

})) # end FeaturesAnalysis
