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
id = "fade_bump",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

sele <- "
SELECT
  geom.AHdist AS AHdist,
  don_site.HBChemType AS don_chem_type,
  acc_site.HBChemType AS acc_chem_type
FROM
  hbonds AS hbond,
  hbond_geom_coords AS geom,
  hbond_sites AS don_site,
  hbond_sites AS acc_site
WHERE
  hbond.struct_id = geom.struct_id AND hbond.hbond_id =  geom.hbond_id AND
  hbond.struct_id = don_site.struct_id AND hbond.don_id = don_site.site_id AND
  hbond.struct_id = acc_site.struct_id AND hbond.acc_id = acc_site.site_id;"

f <-  query_sample_sources(sample_sources, sele)

plot_parts <- list(
   theme_bw(),
   geom_line(aes(x, y, colour=sample_source), size=1.5),
   geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)),
   labs(x=expression(paste('Acceptor -- Proton Distance (', ring(A), ')')),
     y="FeatureDensity"))

get_density <- function(data, ids){
  estimate_density_1d(
    data = data,
    ids = ids,
    variable = "AHdist",
    n_pts=300,
    weight_fun = radial_3d_normalization,
    histogram=TRUE)
}

################################
plot_id <- "fade_first_bump_all"
################################
dens <- get_density(f[1.8 <= f$AHdist & f$AHdist <= 2.0, ], c("sample_source"))
ggplot(data=dens) +
  geom_vline(xintercept = 1.9, colour="darkgray") +
  plot_parts +
  ggtitle("Hydrogen Bond A-H Distance First Fade Bump\nNormalized for equal Weight per Unit Angle") +
  scale_x_continuous( breaks=c(1.85, 1.9, 1.95))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


##########################################
plot_id <- "fade_first_bump_don_chem_type"
##########################################
dens <- get_density(f[1.8 <= f$AHdist & f$AHdist <= 2.0, ], c("sample_source", "don_chem_type"))
ggplot(data=dens) +
  geom_vline(xintercept = 1.9, colour="darkgray") +
  plot_parts +
  facet_wrap(~don_chem_type) +
  ggtitle("Hydrogen Bond A-H Distance First Fade Bump by Donor Chem Type\nNormalized for equal Weight per Unit Angle") +
  scale_x_continuous( breaks=c(1.85, 1.9, 1.95))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


##########################################
plot_id <- "fade_first_bump_acc_chem_type"
##########################################
dens <- get_density(f[1.8 <= f$AHdist & f$AHdist <= 2.0, ], c("sample_source",  "acc_chem_type"))
ggplot(data=dens) +
  geom_vline(xintercept = 1.9, colour="darkgray") +
  plot_parts +
  facet_wrap(~acc_chem_type) +
  ggtitle("Hydrogen Bond A-H Distance First Fade Bump by Acceptor Chem Type\nNormalized for equal Weight per Unit Angle") +
  scale_x_continuous( breaks=c(1.85, 1.9, 1.95))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


######################################
plot_id <- "fade_first_bump_chem_type"
######################################
dens <- get_density(f[1.8 <= f$AHdist & f$AHdist <= 2.0, ], c("sample_source", "don_chem_type", "acc_chem_type"))
ggplot(data=dens) +
  geom_vline(xintercept = 1.9, colour="darkgray") +
  plot_parts +
  facet_grid(don_chem_type ~ acc_chem_type) +
  ggtitle("Hydrogen Bond A-H Distance First Fade Bump by Chem Type\nNormalized for equal Weight per Unit Angle") +
  scale_x_continuous( breaks=c(1.82, 1.9, 1.98))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


################################
plot_id <- "fade_last_bump_all"
################################
dens <- get_density(f[2.2 <= f$AHdist & f$AHdist <= 2.4, ], c("sample_source"))
ggplot(data=dens) +
  geom_vline(xintercept = 2.3, colour="darkgray") +
  plot_parts +
  ggtitle("Hydrogen Bond A-H Distance Last Fade Bump\nNormalized for equal Weight per Unit Angle") +
  scale_x_continuous( breaks=c(2.25, 2.1, 2.45))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)

#########################################
plot_id <- "fade_last_bump_don_chem_type"
#########################################
dens <- get_density(f[2.2 <= f$AHdist & f$AHdist <= 2.4, ], c("sample_source", "don_chem_type"))
ggplot(data=dens) +
  geom_vline(xintercept = 2.3, colour="darkgray") +
  plot_parts +
  facet_wrap(~don_chem_type) +
  ggtitle("Hydrogen Bond A-H Distance Last Fade Bump by Donor Chem Type\nNormalized for equal Weight per Unit Angle") +
  scale_x_continuous( breaks=c(2.25, 2.1, 2.45))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


#########################################
plot_id <- "fade_last_bump_acc_chem_type"
#########################################
dens <- get_density(f[2.2 <= f$AHdist & f$AHdist <= 2.4, ], c("sample_source", "acc_chem_type"))
ggplot(data=dens) +
  geom_vline(xintercept = 2.3, colour="darkgray") +
  plot_parts +
  facet_wrap(~acc_chem_type) +
  ggtitle("Hydrogen Bond A-H Distance Last Fade Bump by Acceptor Chem Type\nNormalized for equal Weight per Unit Angle") +
    scale_x_continuous( breaks=c(2.25, 2.1, 2.45))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


#####################################
plot_id <- "fade_last_bump_chem_type"
#####################################
dens <- get_density(f[2.2 <= f$AHdist & f$AHdist <= 2.4, ], c("sample_source", "don_chem_type", "acc_chem_type"))
ggplot(data=dens) +
  geom_vline(xintercept = 2.3, colour="darkgray") +
  plot_parts +
  facet_grid(don_chem_type ~ acc_chem_type) +
  ggtitle("Hydrogen Bond A-H Distance Last Fade Bump by Chem Type\nNormalized for equal Weight per Unit Angle") +
  scale_x_continuous( breaks=c(2.22, 2.38))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


################################
plot_id <- "fade_range_bump_all"
################################
dens <- get_density(f[1.85 <= f$AHdist & f$AHdist <= 2.35, ], c("sample_source"))
ggplot(data=dens) +
#  geom_vline(xintercept = 1.9, colour="darkgray") +
#  geom_vline(xintercept = 2.3, colour="darkgray") +
  plot_parts +
  ggtitle("Hydrogen Bond A-H Distance Fade Range\nNormalized for equal Weight per Unit Angle") +
  scale_x_continuous( breaks=c(1.9, 2.1, 2.3))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


##########################################
plot_id <- "fade_range_bump_don_chem_type"
##########################################
dens <- get_density(f[1.85 <= f$AHdist & f$AHdist <= 2.32, ], c("sample_source", "don_chem_type"))
ggplot(data=dens) +
  geom_vline(xintercept = 1.9, colour="darkgray") +
  geom_vline(xintercept = 2.3, colour="darkgray") +
  plot_parts +
  facet_wrap( ~ don_chem_type) +
  ggtitle("Hydrogen Bond A-H Distance Fade Range by Donor Chem Type\nNormalized for equal Weight per Unit Angle") +
  scale_x_continuous( breaks=c(1.9, 2.3))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


##########################################
plot_id <- "fade_range_bump_acc_chem_type"
##########################################
dens <- get_density(f[1.85 <= f$AHdist & f$AHdist <= 2.35, ], c("sample_source", "acc_chem_type"))
ggplot(data=dens) +
  geom_vline(xintercept = 1.9, colour="darkgray") +
  geom_vline(xintercept = 2.3, colour="darkgray") +
  plot_parts +
  facet_wrap( ~ acc_chem_type) +
  ggtitle("Hydrogen Bond A-H Distance Fade Range by Acceptor Chem Type\nNormalized for equal Weight per Unit Angle") +
  scale_x_continuous( breaks=c(1.9, 2.3))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


######################################
plot_id <- "fade_range_bump_chem_type"
######################################
dens <- get_density(f[1.85 <= f$AHdist & f$AHdist <= 2.35, ], c("sample_source", "don_chem_type", "acc_chem_type"))
ggplot(data=dens) +
  geom_vline(xintercept = 1.9, colour="darkgray") +
  geom_vline(xintercept = 2.3, colour="darkgray") +
  plot_parts +
  facet_grid(don_chem_type ~ acc_chem_type) +
  ggtitle("Hydrogen Bond A-H Distance Fade Range by Chem Type\nNormalized for equal Weight per Unit Angle") +
  scale_x_continuous( breaks=c(1.9, 2.3))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)



dens

})) # end FeaturesAnalysis
