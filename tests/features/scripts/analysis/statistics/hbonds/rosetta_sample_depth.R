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
id = "rotamer_recovery_by_secondary_structure",
author = "Matthew O'Meara",
brief_description = "",
long_description = "
The purpose of this script is to assess how accurately
distributions generated from rosetta decoys represent
distributions of from rosetta decoys in general. The strategy is
cross validation to assess the likelihood of sampling a set of
training features from an feature distribution estimated from a
training set of features.

A consern is that decoys generated for the same sequence are
similar. Therefore the cross validation grouping will be done over
the input sequences.",

feature_reporter_dependencies = c("StructureFeatures", "HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
	geom.AHdist, geom.cosBAH, geom.cosAHD, geom.chi,
	acc_site.HBChemType AS acc_chem_type,
	don_site.HBChemType AS don_chem_type,
	native_tag.native_tag,
	structure.tag
FROM
	hbond_geom_coords AS geom,
	hbonds AS hbond,
	hbond_sites AS don_site,
	hbond_sites AS acc_site,
	native_tags AS native_tag,
	structures as structure
WHERE
	hbond.struct_id = geom.struct_id AND hbond.hbond_id   = geom.hbond_id AND
	hbond.struct_id = don_site.struct_id AND hbond.don_id = don_site.site_id AND
	hbond.struct_id = acc_site.struct_id AND hbond.acc_id = acc_site.site_id AND
	don_site.HBChemType != 'hbdon_PBA' AND
	acc_site.HBChemType != 'hbacc_PBA' AND
	abs(don_site.resNum - acc_site.resNum) > 5 AND
	structure.struct_id = hbond.struct_id AND
	native_tag.struct_id = hbond.struct_id AND
	structure.struct_id = hbond.struct_id;"

f <- query_sample_sources(sample_sources, sele)

f$don_chem_type <- factor(f$don_chem_type,
	levels = c("hbdon_IMD", "hbdon_IME", "hbdon_GDE", "hbdon_GDH",
		"hbdon_AHX", "hbdon_HXL", "hbdon_IND", "hbdon_AMO", "hbdon_CXA", "hbdon_PBA"),
	labels = c("dIMD: h", "dIME: h", "dGDE: r", "dGDH: r",
		"dAHX: y", "dHXL: s,t", "dIND: w", "dAMO: k", "dCXA: n,q", "dPBA: bb"))

f$acc_chem_type <- factor(f$acc_chem_type,
	levels = c("hbacc_IMD", "hbacc_IME", "hbacc_AHX", "hbacc_HXL",
		"hbacc_CXA", "hbacc_CXL", "hbacc_PBA"),
	labels = c("aIMD: h", "aIME: h", "aAHX: y", "aHXL: s,t",
		"aCXA: n,q", "aCXL: d,e", "aPBA: bb"))

nstruct <- seq(5, 120, 5)
z <- adply(nstruct, 1, function(n){
	sub_f <- f[as.numeric(substr(f$tag, 10, 100)) <= n,]
	z <- cross_validate_statistics(sub_f,
		"native_tag", c("sample_source"), "AHdist", histogram_kl_divergence)
	z$n <- n
	z
})

plot_id <- "hbond_AHdist_rosetta_nstruct_cross_validation"
p <- ggplot(z) + theme_bw() +
	geom_boxplot(aes(x=n, y=KL.Divergence, group=n)) +
#	geom_smooth(aes(ymin=KL.mean - KL.var, ymax=KL.mean+KL.var), stat="identity") +
	ggtitle(("Kullback-Leibler divergence of full cross validation of number of Rosetta predictions per sequence\nComputed for HBond A-H distances for sidechain-sidechain bonds with seq. sep. > 5")) +
	scale_x_continuous("Number of Rosetta predictions per Sequence") +
	scale_y_continuous("KL-Divergence (50 bin)")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

nstruct <- seq(5, 120, 5)
z <- adply(nstruct, 1, function(n){
	sub_f <- f[as.numeric(substr(f$tag, 10, 100)) <= n,]
	z <- cross_validate_statistics(sub_f,
		"native_tag", c("sample_source"), "AHdist", histogram_kl_divergence)
	z$n <- n
	z
})

plot_id <- "hbond_cosBAH_rosetta_nstruct_cross_validation"
p <- ggplot(z) + theme_bw() +
	geom_boxplot(aes(x=n, y=KL.Divergence, group=n)) +
#	geom_smooth(aes(ymin=KL.mean - KL.var, ymax=KL.mean+KL.var), stat="identity") +
	ggtitle(("Kullback-Leibler divergence of full cross validation of number of Rosetta predictions per sequence\nComputed for HBond A-H distances for sidechain-sidechain bonds with seq. sep. > 5")) +
	scale_x_continuous("Number of Rosetta predictions per Sequence") +
	scale_y_continuous("KL-Divergence (50 bin)")
save_plots(self, plot_id, sample_sources, output_dir, output_formats)




plot_id <- "hbond_AHdist_by_chem_type_nstruct_cross_validation"
d_ply(f, .(sample_source), function(sub_f){
	ss_id <- f$sample_source[1]
	ss <- sample_sources[sample_sources$sample_source == ss_id, ]


	z <- cross_validate_statistics(f,
		"native_tag", c("sample_source", "don_chem_type", "acc_chem_type"), "AHdist", histogram_kl_divergence)

	z$chem_type <- paste(z$don_chem_type, z$acc_chem_type, sep="_")
	p <- ggplot(data=z) + theme_bw() +
		geom_line(aes(x=chem_type, y=KL.Divergence, colour=native_tag)) +
		ggtitle(paste("Hydrogen Bonds A-H Distance by Cross Validation Group\nnormalized for equal weight per unit distance ss_id:", ss_id)) +
		scale_y_continuous("KL Divergence (50 bins)")
		scale_x_continuous("Donor and Acceptor Chemical Type")
	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
})






plot_id <- "hbond_AHdist_by_chem_type_nstruct_cross_validation"
d_ply(f, .(sample_source), function(sub_f){
	ss_id <- f$sample_source[1]
	ss <- sample_sources[sample_sources$sample_source == ss_id, ]

	native_tags <- unique(f$native_tag)
	z <- 1:length(cv_groups) %% 4
	names(z) <- native_tags
	f$cv_group <- z[f$native_tag]

	dens <- estimate_density_1d(
		f, c("cv_group", "don_chem_type", "acc_chem_type"),
		"AHdist", weight_fun = radial_3d_normalization)

	p <- ggplot(data=dens) + theme_bw() +
		geom_line(aes(x=x, y=y, colour=cv_group)) +
		geom_indicator(aes(indicator=counts, colour=cv_group, group=cv_group)) +
		ggtitle(paste("Hydrogen Bonds A-H Distance by Cross Validation Group\nnormalized for equal weight per unit distance ss_id:", ss_id)) +
		facet_grid(don_chem_type ~ acc_chem_type) +
		scale_y_continuous("FeatureDensity", limits=c(0,6), breaks=c(1,3,5)) +
		scale_x_continuous(expression(paste('Acceptor -- Proton Distance (', ring(A), ')')), limits=c(1.4,2.7), breaks=c(1.6, 1.9, 2.2, 2.6))

	if(nrow(sample_sources) <= 3){
		p <- p + theme(legend.position="bottom", legend.direction="horizontal")
	}

	save_plots(self, plot_id, sample_sources, output_dir, output_formats)
})

})) # end FeaturesAnalysis
