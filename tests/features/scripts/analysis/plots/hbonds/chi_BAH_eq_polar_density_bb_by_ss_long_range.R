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
id = "chi_BAH_eq_polar_density_bb_by_ss_long_range",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
	geom.cosBAH,
	geom.chi,
	don_ss.dssp AS don_ss, acc_ss.dssp AS acc_ss
FROM
	hbond_geom_coords AS geom,
	hbonds AS hbond,
	hbond_sites AS don_site,
	hbond_sites AS acc_site,
	residue_secondary_structure AS don_ss,
	residue_secondary_structure AS acc_ss
WHERE
	hbond.struct_id  = geom.struct_id     AND hbond.hbond_id = geom.hbond_id    AND
	hbond.struct_id  = don_site.struct_id AND hbond.don_id   = don_site.site_id AND
	hbond.struct_id  = acc_site.struct_id AND hbond.acc_id   = acc_site.site_id AND
	don_ss.struct_id = don_site.struct_id AND don_ss.resNum  = don_site.resNum  AND
	acc_ss.struct_id = acc_site.struct_id AND acc_ss.resNum  = acc_site.resNum  AND
	ABS(don_site.resNum - acc_site.resNum) > 5 AND
	don_site.HBChemType = 'hbdon_PBA' AND acc_site.HBChemType = 'hbacc_PBA';"
f <- query_sample_sources(sample_sources, sele)

# The coding for loop and irregular is " ", change this to "C"
levels(f$acc_ss)[levels(f$acc_ss) == " "] <- "C"
levels(f$don_ss)[levels(f$don_ss) == " "] <- "C"

####                                                             ####
# This relabelling is deprecated use the 'dssp_codes' table instead!#
###                                                              ####
# Add descriptive names to the DSSP codes
f$acc_ss_name <- factor(f$acc_ss,
	levels = c("H", "E", "T", "G", "B", "S", "I", "C"),
	labels = c('H: a-Helix', 'E: b-Sheet', 'T: HB Turn', 'G: 3/10 Helix',
		'B: b-Bridge', 'S: Bend', 'I: pi-Helix','C: Irregular'))
f$don_ss_name <- factor(f$don_ss,
	levels = c("H", "E", "T", "G", "B", "S", "I", "C"),
	labels = c('H: a-Helix', 'E: b-Sheet', 'T: HB Turn', 'G: 3/10 Helix',
		'B: b-Bridge', 'S: Bend', 'I: pi-Helix','C: Irregular'))

f <- f[!is.na(f$don_ss_name),]
f <- f[!is.na(f$acc_ss_name),]

#equal area projection
f <- transform(f,
	capx = 2*sin(acos(cosBAH)/2)*cos(chi),
	capy = 2*sin(acos(cosBAH)/2)*sin(chi))
#capx_limits <- range(f$capx); capy_limits <- range(f$capy)
capx_limits <- c(-1.5,1.5); capy_limits <- capx_limits;

f <- ddply(f, c("sample_source", "don_ss_name", "acc_ss_name"),
	transform, counts = length(sample_source))

d_ply(sample_sources, .(sample_source), function(sample_source){
	ss_id <- sample_source$sample_source[1]
	plot_id = paste("chi_sinBAH_eq_polar_density_bb_by_ss_long_range", ss_id, sep="_")

	ggplot(data=subset(f, sample_source == ss_id)) + theme_bw() +
		stat_bin2d(aes(x=capx, y=capy, fill=log(..density..)), binwidth=c(.06,.06)) +
		polar_equal_area_grids_bw() +

		geom_indicator(aes(indicator=counts)) +
		facet_grid(don_ss_name ~ acc_ss_name) +
		ggtitle(paste("Backbone-Backbone HBonds with Sequence Separation > 5: CHI vs BAH Angles by DSSP\nEqual Coordinate Projection   Sample Source: ", ss_id, sep="")) +
		scale_x_continuous('2*sin(BAH/2) * cos(CHI)', limits=capx_limits, breaks=c(-1, 0, 1)) +
		scale_y_continuous('2*sin(BAH/2) * sin(CHI)', limits=capy_limits, breaks=c(-1, 0, 1)) +
		scale_fill_gradientn('log(Density)', colours=jet.colors(15))
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
})


})) # end FeaturesAnalysis
