# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


feature_analyses <- c(feature_analyses, new("FeaturesAnalysis",
id = "EXAMPLE_PLOT",
author = "Matthew O'Meara",

brief_description = "A simple demonstration of the Grammar of Graphics plotting functionality.  Make a bar graph of the number of residues in each supplied sample source",

long_description = "
 This script can be run like this:

   mini/test/scientific/cluster/features/compare_sample_sources.R --script scripts/analysis/plots/EXAMPLE_PLOT.R features_<sample_source_id1>.db3 ...

 It will generate:

   rosetta_tests/features/build/output_web_raster/EXAMPLE_PLOT_<date_code>_<sample_source_id1>[_<sample_source_id2> ...].png",

keywords = c("bar_chart"),
feature_reporter_dependencies = c("ResidueFeatures"),

run=function(self, sample_sources, output_dir, output_formats){

# The SQL query is applied to each sample source and the resulting
# tables are appended together with an additional sample_source column
# and returned as a single data.frame object.
sele <- "SELECT count(*) AS num_residues FROM residues;"
features <- query_sample_sources(sample_sources, sele)

# Add code to analyze the features data.frame here such as generating
# densities, computing additional geometric measurements, etc.

# The ggplot plotting system maps the columns of the data.frame to aesthetic parts of the plot
# see http://had.co.nz/ggplot2/
plot <- ggplot(data=features)

# theme_bw is a clean black and white theme
plot <- plot + theme_bw()

#geom_bar is a bar graph layer, here the sample_source and
#num_residues are mapped to the x, y, and color aesthetics.
plot <- plot + geom_bar(aes(x=sample_source, y=num_residues, fill=sample_source))

# Set title and axis labels
plot <- plot + ggtitle("Dataset Size") + labs(x="Number of Residues", y="Count")

# The resulting filename out of save_plots is <plot_id>_<date_code>_<sample_source1>_..._.<extension>
plot_id <- "EXAMPLE_PLOT"

# See './compare_sample_sources.R --help' about specifying the output_dir and output_formats.

save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
