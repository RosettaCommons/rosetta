# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

check_setup()

sele <-"

SELECT
  acc_site.HBChemType as acc_chem_type,
  count(acc_site.HBChemType) AS acc_chem_type_count
FROM    
  hbonds AS hbond,
  hbond_sites AS acc_site
WHERE
  hbond.struct_id == acc_site.struct_id AND hbond.acc_id == acc_site.site_id
GROUP BY
  acc_site.HBChemType;"

f <- query_sample_sources(sample_sources, sele)

plot_id <- "hbonds_num_by_acc_chem_type"
ggplot(f, aes(acc_chem_type, log(acc_chem_type_count))) + theme_bw() + 
  geom_bar(aes(fill=sample_source), position="dodge") +
  opts(title = "HBond Counts for each Acceptor Type") +
  labs(x="Acceptor Type", y="Counts")
save_plots(plot_id, sample_sources, output_dir, output_formats)
