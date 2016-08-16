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
id = "rotamer_bins",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("ResidueFeatures", "RotamerFeatures", "ResidueType"),
run=function(self, sample_sources, output_dir, output_formats){


sele <-"
SELECT
	res.name3 AS res_type,
	rot.rotamer_bin AS rotamer_bin,
	count(*) AS rot_bin_count
FROM
	residues AS res,
	residue_rotamers AS rot
WHERE
	rot.struct_id = res.struct_id AND
	rot.residue_number = res.resNum
GROUP BY
	res_type, rotamer_bin
ORDER BY
	res_type, rot_bin_count DESC;"

f <- query_sample_sources(sample_sources, sele)
f$rotamer_bin <- factor(f$rotamer_bin)

f <- ddply(f, .(res_type), function(df){
	if(sum(df$rot_bin_count) < 50){
		return(data.frame())
	} else {
		return(df)
	}
})

f <- ddply(f, .(sample_source, res_type), transform,
	rot_bin_fraction=rot_bin_count/sum(rot_bin_count))

first_sample_source_id <- as.character(sample_sources$sample_source[1])
zz <- ddply(f, .(res_type), function(df) {
	df_wide <- cast(df, res_type + rotamer_bin ~ sample_source, value="rot_bin_fraction")
	df_wide <- df_wide[order(df_wide[,first_sample_source_id], decreasing=T),]
	df_wide$rank <- seq_along(df_wide$rotamer_bin)
	attr(df_wide, "idvars") <- c("res_type", "rotamer_bin", "rank")
	df_long <- melt(df_wide)
	names(df_long) <- c(
		"res_type", "rotamer_bin", "rank", "rot_bin_fraction", "sample_source")
	df_long
})



plot_id <- "rotamer_bin_counts"
p <- ggplot(data=zz) + theme_bw() +
	geom_line(aes(x=rank, y=rot_bin_fraction, colour=sample_source)) +
	facet_wrap(~res_type, scales="free") +
	ggtitle(("Rotamer Bin Counts", sep="")) +
	scale_y_continuous("Rotamer Bin Density") +
	scale_x_continuous("Rotamer Bin Rank")
if(nrow(sample_sources) <= 3){
  p <- p + theme(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
