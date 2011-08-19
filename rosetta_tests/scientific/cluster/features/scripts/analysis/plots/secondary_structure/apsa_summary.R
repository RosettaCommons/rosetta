# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

check_setup()

plot_id <- "apsa_summary"

sele <- "
SELECT
  secondary_struct,
  count(secondary_struct) AS count
FROM
  apsa
GROUP BY
  secondary_struct;"

all <-  query_sample_sources(sample_sources, sele)

p <- ggplot(data=all, aes(x=secondary_struct, y=count, fill=sample_source))
p <- p + geom_bar(position="dodge") + coord_flip()
p <- p + opts(title = "Secondary Structure Diversity\nAutomated Protein Structure Analysis Method")
p <- p + labs(x = "Type of Secondary Structure",
              y = "Count")
p <- p + theme_bw()

save_plots(plot_id, sample_sources, output_dir, output_formats)
