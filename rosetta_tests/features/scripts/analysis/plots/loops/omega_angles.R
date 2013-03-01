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
id = "omega_angles",
filename = "scripts/analysis/plots/loops/omega_angles.R",
author = "Brian D. Weitzner",
brief_description = "",
long_description = "
This features analysis computes the correlation between the six dimensions of a coordinate frame transformation pairwise
",

feature_reporter_dependencies = c("loop_anchor_features"),
run=function(self, sample_sources, output_dir, output_formats){


# Loop anchor transfrom parameters for all candidate loops of a given length
sele <-paste("SELECT omega FROM loop_anchor_transforms;")

f <- query_sample_sources(sample_sources, sele)
p <- ggplot(f, aes(x=omega)) +
    geom_histogram(aes(y=..density..), fill="grey", colour="grey50") + 
    geom_density(colour="black") +
    ggtitle((paste(omega, " in antibodies")))
save_plots(self, paste("omega_scatterplot_matrix", sep = "_"), sample_sources, output_dir, output_formats)

})) # end FeaturesAnalysis