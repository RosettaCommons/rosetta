# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

library(plyr)
library(reshape2)

feature_analyses <- c(feature_analyses, methods::new("FeaturesAnalysis",
id = "sequence_recovery_summary",
author = "Matthew O'Meara",
brief_description = "",
long_description = "",

# Assume that sample_sources[1] is the database of natives
# and that sample_sources[2] is the database of redesigns

feature_reporter_dependencies = c("ResidueFeatures", "PdbDataFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
  ss.dssp AS sec_struct,
  res.name3 AS res_type
FROM
  residue_secondary_structure as ss,
  residues as res
WHERE
  ss.struct_id == res.struct_id AND ss.resNum == res.resNum;
"


f <- query_sample_sources(sample_sources, sele)
tf <- as.data.frame(table(f))
print(summary(tf))
print(class(tf))
tf <- merge( tf, ddply( tf, .(sample_source,sec_struct), function(subf) { data.frame(total = sum(subf$Freq)) } ) )
print(summary(tf))
print(class(tf))
tf <- ddply(tf,.(sample_source,sec_struct,res_type), function(subf) { data.frame(Freq = subf$Freq / subf$total)})
tf$Freq <- round(as.numeric(tf$Freq),3)
tf_out <- acast(tf, sample_source~res_type~sec_struct, value.var="Freq")
print(tf_out)

tf2 <- ddply( as.data.frame(table(f)), .(sample_source,res_type), function(subf) { data.frame(counts = sum(subf$Freq)) })
tf2 <- merge( tf2, ddply( tf2, .(sample_source),function(subf) { data.frame(total= sum(subf$counts)) }))
tf2 <- ddply( tf2, .(sample_source,res_type), function(subf) { data.frame( freq = subf$counts / subf$ total )})
tf2$freq <- round(as.numeric(tf2$freq),3)
tf2_out <- acast(tf2, sample_source~res_type,value.var="freq")

print("overall AA recovery rates")
print(tf2_out)

table_id <- "amino_acid_by_ss_breakdown"
table_title <- "Which amino acids are designed into which secondary structures"
save_tables( self, tf, table_id, sample_sources, output_dir, output_formats,
  caption=table_title, caption.placement="top" )

#sele <- "
#DROP TABLE nchi;
#DROP TABLE chi_recovered;"
#query_sample_sources(sample_sources, sele, warn_zero_rows=F)

#f$avg_recovery <- round(as.numeric(f$avg_recovery), 2)
#
#table_id <- "sequence_recovery_summary"
#table_title <- "Average Rotamer Recovery, BFactor < 30\n"
#save_tables(self,
#	f, table_id,
#	sample_sources, output_dir, output_formats,
#	caption=table_title, caption.placement="top")

})) # end FeaturesAnalysis
