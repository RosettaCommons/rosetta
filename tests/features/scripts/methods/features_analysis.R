# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

setClass("FeaturesAnalysis",
	representation(
		# Unique identifier for the features analysis
		# e.g. if unspecified use the filename instead
		id = "character",

		# For feature analyses that come from scripts, this is
		# the filename of the script
		filename = "character",

		# The author of the script
		author = "character",

		# A breif description of the features this script analyses
		brief_description = "character",

		# A description of the purpose of the analysis script
		# and any specific transformations, or methods for
		# plotting etc.
		long_description = "character",

		# Keywords can be used to locate analysis scripts
		keywords = "vector",

		# Which FeatureReporter objects need to be extracted
		# into the feature database for each sample source in
		# order to run this features analysis
		# See https://wiki.rosettacommons.org/index.php/FeaturesDatabaseSchema
		feature_reporter_dependencies = "vector",

		# The features analysis to be run
		run="function"),
	prototype(
		id = "",
		filename="",
		author="",
		brief_description = "Class for features analysis script.",
		long_description = "Instantiate this class to create feature a features analysis.  This is usually done in a features analysis script. See scripts/analysis/plots/EXAMPLE_PLOT.R",
		keywords = NA,
		feature_reporter_dependencies = NA,
		run=function(self, sample_sources, output_dir, output_formats){print("Implement this function to preform the features analysis script")}))



# Here is an example of how to create a FeaturesAnalysis object
#z <- new("FeaturesAnalysis",
#	id ="test_features_analysis",
#	brief_description="breif_description",
#	long_description="long_description",
#	feature_reporter_dependencies=c("abc"),
#	run=function(self, sample_sources, output_dir, output_formats){print("Running test features analysis")})


