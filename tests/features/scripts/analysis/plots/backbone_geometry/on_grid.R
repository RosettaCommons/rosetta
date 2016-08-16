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
id = "on_grid",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("ResidueFeatures", "ProteinBackboneTorsionAngleFeatures", "ResidueSecondaryStructureFeatures", "PdbDataFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
	bb.phi, bb.psi
FROM
	residue_pdb_confidence AS res_conf,
	protein_backbone_torsion_angles AS bb
WHERE
	bb.struct_id = res_conf.struct_id AND
	bb.resNum == res_conf.residue_number AND
	res_conf.max_temperature < 30 AND
	bb.phi != 0 AND bb.psi != 0;"

f <- query_sample_sources(sample_sources, sele)

f <- rbind(
	data.frame(
		sample_source="Uniform Null",
		phi=runif(1000000, 0, 2.5),
		psi=runif(1000000, 0, 2.5)),
	f)

f$off_grid <- with(f,
	pmin(phi %% 5, (5-phi) %% 5, psi %% 5, (5-psi) %% 5))

normalization <- function(x){
  a <- pmin(abs(x) %% 5, (5-abs(x)) %% 5)
  w <- 1/(20-8*a)
  w/sum(w)
}

dens <- estimate_density_1d_reflect_boundary(
	f,
	c("sample_source"), "off_grid",
	weight_fun=normalization,
	reflect_left=T, left_boundary=0,
	reflect_right=T, right_boundary=2.5,
	adjust = .1)

at_boundary <- ddply(f, .(sample_source), function(df){
	x <- nrow(df[df$off_grid < .05 ,])/nrow(df)
	data.frame(
		sample_source = df$sample_source[1],
		at_boundary = x,
		at_boundary_string = paste("Percent < .05A  ", as.character(round(x*100, 0)), "%"), sep="")
})

cat("at_boundary: off_grid < .05 degrees")
print(at_boundary)

#plot_id <- "on_grid"
#ggplot(data=dens) +
#	theme_bw() +
#	geom_line(aes(x=x, y=y, colour=sample_source)) +
##	geom_vline(aes(x=.05)) +
#	geom_indicator(
#		data=at_boundary,
#		aes(indicator=at_boundary_string,
#		group=sample_source, colour=sample_source)) +
#	ggtitle(paste(
#		"Min Distance From Grid Boundary; B-Factor < 30", sep=""),
#		legend.position = "none") +
#	scale_y_continuous("Feature Density") +
#	scale_x_continuous("Min Distance From Grid (degrees)", limits=c(0,2.5)) +
#save_plots(self, plot_id, sample_sources, output_dir, output_formats)
#
#plot_id <- "on_grid_clean"
#ggplot(data=dens) +
#	theme_bw() +
#	geom_line(aes(x=x, y=y, colour=sample_source), size=1.5) +
#	scale_y_continuous("Feature Density") +
#	theme(legend.position=c(.8, .8)) +
#	scale_colour_grey("Energy Function") +
#	scale_x_continuous("Min Distance From Grid (degrees)", limits=c(0,.25)) +
#save_plots(self, plot_id, sample_sources, output_dir, output_formats)



# In a uniform distribution, the probability of "off_grid" is
# P(x) = (20-8x)/25

# The quantile function for off_grid is the inverse of the CDF:
# p = Integrate[(20 - 8x)25, (x, 0, q)] = -4(q-5)q/25
# q = -5/2(sqrt(1-p)-1)

q_theory <- function(p){ -5/2 * ( sqrt(1-p) - 1)}

f <- ddply(f, .(sample_source), function(df){
	df$counts <- nrow(df)
	df
})

sub_f <- ddply(f, .(sample_source), function(df){
	df[sample(nrow(df), 50000),]
})

#plot_id <- "on_grid_quantile_quantile"
#ggplot() +
#	theme_bw() +
#	stat_qq(data=sub_f, aes(sample=off_grid, colour=sample_source),
#		geom="line", dist=q_theory) +
#	geom_indicator(data=f,
#		aes(indicator=counts, group=sample_source, colour=sample_source),
#		xpos="left") +
#	ggtitle("Min Distance From Grid Boundary, Quantile vs Quantile; B-Factor < 30") +
#	scale_x_continuous(
#		expression(paste("Uniform Distribution Quantile (degrees)", sep="")),
#		limits=c(0, 2.5)) +
#	scale_y_continuous(
#		expression(paste("Sample Distribution Quantile (degrees)", sep="")),
#		limits=c(0, 2.5))
#save_plots(self, plot_id, sample_sources, output_dir, output_formats)


plot_id <- "on_grid_quantile_quantile_clean"
ggplot() +
	theme_bw() +
	stat_qq(data=sub_f, aes(sample=off_grid, colour=sample_source),
		geom="line", dist=q_theory, size=1.5) +
	theme(
		legend.position = c(.25,.8)) +
	scale_colour_grey("Energy Function") +
	scale_x_continuous("Uniform Distribution Quantile (degrees)",
		limits=c(0, 2.5)) +
	scale_y_continuous("Sample Distribution Quantile (degrees)",
		limits=c(0, 2.5))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)




})) # end FeaturesAnalysis
