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
This features analysis computes the correlation between the six dimensions of a coordinate frame transformation pairwise
",

feature_reporter_dependencies = c("loop_anchor_features"),
run=function(self, sample_sources, output_dir, output_formats){

plotmatrix <- function (input_data, mapping = aes(), columns = names(input_data), colour = "black") 
{
  data <- input_data[,columns]		
  grid <- expand.grid(x = 1:ncol(data), y = 1:ncol(data))
  grid <- subset(grid, x != y)
  all <- do.call("rbind", lapply(1:nrow(grid), function(i) {
	 xcol <- grid[i, "x"]
	 ycol <- grid[i, "y"]
	 data.frame(xvar = names(data)[ycol], yvar = names(data)[xcol], 
	 x = data[, xcol], y = data[, ycol], input_data)
  }))
  all$xvar <- factor(all$xvar, levels = names(data))
  all$yvar <- factor(all$yvar, levels = names(data))

  densities <- do.call("rbind", lapply(1:ncol(data), function(i) {
	 data.frame(xvar = names(data)[i], yvar = names(data)[i], 
	 x = data[, i], subset(input_data, select = !(names(f) %in% columns)))
  }))
  mapping <- plyr::defaults(mapping, aes_string(x = "x", y = "y"))
  class(mapping) <- "uneval"
  p <- ggplot(all, mapping) + facet_grid(xvar ~ yvar, scales = "free") + 
	 geom_point(na.rm = TRUE, size=0.3) + stat_density(aes(x = x, 
	 y = ..scaled.. * diff(range(x)) + min(x)), data = densities, 
	 position = "identity", geom = "line") + 
	 opts(title = paste("Scatterplot matrix for loop anchor transform degrees of freedom for", loop_length, "residue loops.", sep = " " )) +
	 scale_y_continuous("Degree of freedom") +
	 scale_x_continuous("Degree of freedom") +
	 coord_equal(ratio=1)

  if(nrow(sample_sources) <= 3){
	 p <- p + opts(legend.position="bottom", legend.direction="horizontal")
  }
}

statement <- c("SELECT x, y, z, phi, psi, theta, alpha, omega FROM loop_anchor_transforms", "SELECT x, y, z, phi, psi, theta, alpha, omega FROM loop_anchor_transforms_three_res", "SELECT loop_anchor_transforms.x, loop_anchor_transforms.y, loop_anchor_transforms.z, loop_anchor_transforms.phi, loop_anchor_transforms.psi, loop_anchor_transforms.theta, loop_anchor_transforms_three_res.alpha, loop_anchor_transforms_three_res.omega FROM loop_anchor_transforms INNER JOIN loop_anchor_transforms_three_res ON loop_anchor_transforms.struct_id = loop_anchor_transforms_three_res.struct_id AND loop_anchor_transforms.residue_begin = loop_anchor_transforms_three_res.residue_begin AND loop_anchor_transforms.residue_end = loop_anchor_transforms_three_res.residue_end")
file_names <- c("one_res", "three_res", "hybrid")
scale_phi <- c(TRUE, FALSE, TRUE)

for (j in 1:3) {
  means_data_frame <- data.frame()
  std_devs_data_frame <- data.frame()
  for (i in 1:30) {
	 loop_length <- i
	 # Loop anchor transfrom parameters for all candidate loops of a given length
	 sele <-paste(statement[j], "WHERE residue_end - residue_begin = ", loop_length, " LIMIT 20000;", sep = " ")

	 f <- query_sample_sources(sample_sources, sele)

	 if(!length(f$sample_source)){ next }

	 if(scale_phi[j]){
		f$phi <- ifelse( f$phi < 0, f$phi + 2*pi, f$phi)
	 } else{
		f$theta <- ifelse( f$theta > pi/2, f$theta - pi, f$theta)
	 }

	 data_columns <- !(names(f) %in% c("sample_source"))

	 p <- plotmatrix(f, aes(colour=sample_source), subset(names(f), data_columns))

	 save_plots(self, paste(file_names[j], "anchor_transform_scatterplot_matrix_for", loop_length,"residue_loops", sep = "_"), sample_sources, output_dir, output_formats)

	 sample_source_names <- levels(unique(f$sample_source))

	 # setup the means data frame on first pass
	 if (length(means_data_frame) != length(f[data_columns])) {
		means_data_frame <- data.frame(t(rep(NA, length(f[data_columns]))))
		colnames(means_data_frame) <- colnames(f[data_columns])
		means_data_frame <- means_data_frame[-1,]
	 }
	 # setup the std devs data frame on first pass
	 if (length(std_devs_data_frame) != length(f[data_columns])) {
		std_devs_data_frame <- data.frame(t(rep(NA, length(f[data_columns]))))
		colnames(std_devs_data_frame) <- colnames(f[data_columns])
		std_devs_data_frame <- std_devs_data_frame[-1,]
	 }

	 for (sample_source_name in sample_source_names) {
		data_from_current_sample_source <- f[f$sample_source == sample_source_name,]
		number_of_loop_anchor_transforms <- ifelse(length(data_from_current_sample_source), length(data_from_current_sample_source[,1]), 0)

		cat(paste("\n",loop_length, " residue loops from ", sample_source_name, ": ", number_of_loop_anchor_transforms, " data points\n", sep=""))
		if (number_of_loop_anchor_transforms) {
		  loop_anchor_transform_data <- data_from_current_sample_source[data_columns]
		  print(formatC(cov.wt(loop_anchor_transform_data)$cov, digits=4, width= 7, format='f'))

		  if(substr(sample_source_name, start=1, stop=10) == "antibodies") {
			 data_frame_dimensions <- length(loop_anchor_transform_data)
			 means_data_frame[loop_length,] <- ifelse(rep(number_of_loop_anchor_transforms, data_frame_dimensions), mean(loop_anchor_transform_data), t(rep(NA, data_frame_dimensions)))
			 std_devs_data_frame[loop_length,] <- ifelse(rep(number_of_loop_anchor_transforms, data_frame_dimensions), sd(loop_anchor_transform_data), t(rep(NA, data_frame_dimensions)))
		  }
		}
	 }
	 cat("\n")
  }

  antibody_lat_means_file_name <- paste(output_dir, paste("LAT_means_", file_names[j], ".csv", sep = ""), sep = "/")
  cat("Antibody Means: written to", antibody_lat_means_file_name, "\n")
  write.csv(means_data_frame, file = antibody_lat_means_file_name)
  cat("\n")

  antibody_lat_std_devs_file_name <- paste(output_dir, paste("LAT_std_devs_", file_names[j], ".csv", sep = ""), sep = "/")
  cat("Antibody Standard Deviations: written to", antibody_lat_std_devs_file_name, "\n")
  write.csv(std_devs_data_frame, file = antibody_lat_std_devs_file_name)
}

})) # end FeaturesAnalysis
