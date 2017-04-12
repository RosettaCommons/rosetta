# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

check_setup()
feature_analyses <- c(feature_analyses, new("FeaturesAnalysis",
id = "percent_overcommitted",
author = "Andrew Leaver-Fay",
brief_description = "Compute the percentage of structures in the sample source that have more
than two donors to an acceptor, or more than one acceptor to a donor.",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT DISTINCT
  struct.tag as tag,
  hbond.struct_id AS struct_id
FROM
  hbonds as hbond,
	structures as struct
WHERE
  ( hbond.accRank > 2 OR hbond.donRank > 1 ) AND
	hbond.struct_id == struct.struct_id;"

f <- query_sample_sources(sample_sources, sele)

f <- ddply(f, .(sample_source),
  transform, counts = length(sample_source))

print(f)


})) # end FeaturesAnalysis
