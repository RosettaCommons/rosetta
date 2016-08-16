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
id = "bifurcated_at_acceptor",
author = "Matthew O'Meara",
brief_description = "The rank of a hydrogen bond at donor site or acceptor site is rank of the relative Rosetta HBond energy of the hydrogen bond at the site.",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
CREATE TEMPORARY TABLE hbs AS SELECT
  acc_site.struct_id AS struct_id,
  acc_site.site_id AS acc_site_id,
  geom.AHDist AS AHdist,
  acc_site.HBChemType AS acc_chem_type,
  don_site.HBChemType AS don_chem_type,
  hbond.accRank AS acc_rank
FROM
  hbonds AS hbond,
  hbond_sites AS acc_site,
  hbond_sites AS don_site,
  hbond_geom_coords AS geom
WHERE
  hbond.struct_id = acc_site.struct_id AND hbond.acc_id = acc_site.site_id AND
  hbond.struct_id = don_site.struct_id AND hbond.don_id = don_site.site_id AND
  hbond.struct_id = geom.struct_id AND hbond.hbond_id = geom.hbond_id;

CREATE INDEX hbs_index ON hbs(struct_id, acc_site_id);

SELECT
  hb1.AHdist AS AHdist1,
  hb2.AHdist AS AHdist2,
  hb1.don_chem_type AS don_chem_type1,
  hb2.don_chem_type AS don_chem_type2
FROM
  hbs AS hb1,
  hbs AS hb2
WHERE
  hb1.struct_id = hb2.struct_id AND
  hb1.acc_site_id = hb2.acc_site_id AND
  hb1.acc_rank > hb2.acc_rank;"


f <- query_sample_sources(sample_sources, sele)

f <- ddply(f, .(sample_source, don_chem_type1, don_chem_type2),
  transform, counts = length(sample_source))

# This is deprecated please use the hbond_chem_types table for the lables instead
# Order the plots better and give more descriptive labels
f$don_chem_type1 <- factor(f$don_chem_type1,
	levels = c("hbdon_IMD", "hbdon_IME", "hbdon_GDE", "hbdon_GDH",
		"hbdon_AHX", "hbdon_HXL", "hbdon_IND", "hbdon_AMO", "hbdon_CXA", "hbdon_PBA"),
	labels = c("dIMD: h", "dIME: h", "dGDE: r", "dGDH: r",
		"dAHX: y", "dHXL: s,t", "dIND: w", "dAMO: k", "dCXA: n,q", "dPBA: bb"))


# This is deprecated please use the hbond_chem_types table for the lables instead
# Order the plots better and give more descriptive labels
f$don_chem_type2 <- factor(f$don_chem_type2,
	levels = c("hbdon_IMD", "hbdon_IME", "hbdon_GDE", "hbdon_GDH",
		"hbdon_AHX", "hbdon_HXL", "hbdon_IND", "hbdon_AMO", "hbdon_CXA", "hbdon_PBA"),
	labels = c("dIMD: h", "dIME: h", "dGDE: r", "dGDH: r",
		"dAHX: y", "dHXL: s,t", "dIND: w", "dAMO: k", "dCXA: n,q", "dPBA: bb"))


plot_id <- "hbond_bifurcated_at_acceptor_AHdist"
d_ply(f, .(sample_source), function(sub_f){
	ss_id <- sub_f$sample_source[1]
	ggplot(sub_f) + theme_bw() +
		geom_point(aes(x=AHdist1, y=AHdist2), size=.4) +
		stat_density2d(aes(x=AHdist1, y=AHdist2), size=.2) +
		geom_indicator(aes(indicator=counts), colour="red") + # TODO: fix colour
		facet_grid(don_chem_type1 ~ don_chem_type2) +
		scale_x_continuous(
			expression(paste('Acceptor -- Hydrogen Distance (', ring(A), ') Worse Energy HBond')),
			limits=c(1.5, 3), breaks=c(1.5, 2.0, 2.5))+
		scale_y_continuous(
			expression(paste('Acceptor -- Hydrogen Distance (', ring(A), ') Better Energy HBond')),
			limits=c(1.5, 3), breaks=c(1.5, 2.0, 2.5))+
                  ggtitle(paste("HBond Pairs Bifurcated at Acceptor   ss_id:",ss_id))
	save_plots(self, plot_id, sample_sources[sample_sources$sample_source == ss_id,], output_dir, output_formats)
})


})) # end FeaturesAnalysis
