# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

library(ggplot2)


library(plyr)

library("scales")

#check_setup()

feature_analyses <- c(feature_analyses, new("FeaturesAnalysis",
id = "nhbonds_per_site",
author = "Andrew Leaver-Fay",
brief_description = "Count the number of hydrogen bonds formed per site conditional on the donor and acceptor chemical types.",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

#source("/home/andrew/GIT/Rosetta/main6/tests/features/scripts/analysis/plots/hbonds/hbond_geo_dim_scales.R")

# source("../hbond_geo_dim_scales.R")

source( "/usr/local/lib/R/site-library/RosettaFeatures/scripts/analysis/plots/hbonds/hbond_geo_dim_scales.R" )

sele <-"
SELECT
  acc_site.HBChemType AS acc_chem_type,
  acc_env.num_hbonds AS num_hbonds
FROM
  hbond_sites AS acc_site,
	hbond_site_environment AS acc_env
WHERE
  acc_env.struct_id == acc_site.struct_id AND acc_env.site_id == acc_site.site_id AND
	acc_site.is_donor == 0
;"


f <- query_sample_sources(sample_sources, sele)

tf <- as.data.frame(table(f))
print(tf)

# tf2 <- cut(tf,tf.acc_chem_type,labels=T)
agg<-aggregate(tf$Freq, by=list(sample_source=tf$sample_source,acc_chem_type=tf$acc_chem_type), FUN=sum )
print(agg)

count = 0
tf2 <- ddply( tf, c("sample_source","acc_chem_type","num_hbonds"), function(subf) {
		data.frame(perc=subf$Freq / agg[ agg$sample_source == subf$sample_source & agg$acc_chem_type == subf$acc_chem_type, "x" ], count=agg[ agg$sample_source == subf$sample_source & agg$acc_chem_type == subf$acc_chem_type, "x" ] )
})
cat(paste("count",count))
print(tf2)



tf2$acc_chem_type_name <- acc_chem_type_name_linear(tf2$acc_chem_type)

plot_id <- "nhbonds_per_site_by_acc_chem_type"
ggplot(tf2, aes(x=num_hbonds,group=sample_source,colour=sample_source,fill=sample_source)) + theme_bw() +
	#geom_indicator(aes(indicator=count)) +
	geom_bar( aes(y=perc), stat="identity", position="dodge" ) +

	facet_grid( ~ acc_chem_type_name ) +
  ggtitle("Number of HBonds Per Site by Acceptor Type") +
  labs(x="Acceptor Type", y="counts") +
	#scale_y_continuous( labels = percent_format() )

save_plots(self, plot_id, sample_sources, output_dir, output_formats)

# ................ Now the donors

sele <-"
SELECT
  don_site.HBChemType AS don_chem_type,
  don_env.num_hbonds AS num_hbonds
FROM
  hbond_sites AS don_site,
	hbond_site_environment AS don_env
WHERE
  don_env.struct_id == don_site.struct_id AND don_env.site_id == don_site.site_id AND
	don_site.is_donor == 1
;"


f <- query_sample_sources(sample_sources, sele)

tf <- as.data.frame(table(f))
print(tf)

# tf2 <- cut(tf,tf.don_chem_type,labels=T)
agg<-aggregate(tf$Freq, by=list(sample_source=tf$sample_source,don_chem_type=tf$don_chem_type), FUN=sum )
print(agg)

count = 0
tf2 <- ddply( tf, c("sample_source","don_chem_type","num_hbonds"), function(subf) {
		data.frame(perc=subf$Freq / agg[ agg$sample_source == subf$sample_source & agg$don_chem_type == subf$don_chem_type, "x" ], count=agg[ agg$sample_source == subf$sample_source & agg$don_chem_type == subf$don_chem_type, "x" ] )
})
cat(paste("count",count))
print(tf2)



tf2$don_chem_type_name <- don_chem_type_name_linear(tf2$don_chem_type)

plot_id <- "nhbonds_per_site_by_don_chem_type"
ggplot(tf2, aes(x=num_hbonds,group=sample_source,colour=sample_source,fill=sample_source)) + theme_bw() +
	#geom_indicator(aes(indicator=count)) +
	geom_bar( aes(y=perc), stat="identity", position="dodge" ) +

	facet_grid( ~ don_chem_type_name ) +
  ggtitle("Number of HBonds Per Site by Donor Type") +
  labs(x="Donor Type", y="counts") +
	#scale_y_continuous( labels = percent_format() )

save_plots(self, plot_id, sample_sources, output_dir, output_formats)



# TEMP ! sele <-"
# TEMP ! SELECT
# TEMP !   don_site.HBChemType AS don_chem_type,
# TEMP !   don_env.num_hbonds AS num_hbonds
# TEMP ! FROM
# TEMP !   hbond_sites AS don_site,
# TEMP ! 	hbond_site_environment AS don_env
# TEMP ! WHERE
# TEMP !   don_env.struct_id == don_site.struct_id AND don_env.site_id == don_site.site_id AND
# TEMP ! 	don_site.is_donor == 1
# TEMP ! ;"
# TEMP !
# TEMP !
# TEMP !
# TEMP ! f2 <- query_sample_sources(sample_sources, sele)
# TEMP !
# TEMP ! f2$don_chem_type_name <- don_chem_type_name_linear(f2$don_chem_type)
# TEMP ! f2 <- na.omit(f2, method="r")
# TEMP !
# TEMP ! ddply(f2, .(sample_source), function(sub_f){
# TEMP ! 	ss_id <- sub_f$sample_source[1]
# TEMP ! 	ss <- sample_sources[sample_sources$sample_source == ss_id,]
# TEMP !
# TEMP ! 	#print("testing1")
# TEMP ! 	#print(ss)
# TEMP ! 	#print("testing2")
# TEMP !   plot_id <- "nhbonds_per_site_by_don_chem_type"
# TEMP !   ggplot(sub_f, aes( x=num_hbonds )) + theme_bw() +
# TEMP !   	#geom_histogram( aes( x=num_hbonds, y=..density..) ) +
# TEMP !   	geom_histogram(aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) +
# TEMP !   	#geom_histogram() +
# TEMP !   	facet_grid( ~ don_chem_type_name ) +
# TEMP !     ggtitle(paste("Number of HBonds Per Site by Donor Type:",ss$sample_source)) +
# TEMP !     labs(x="Donor Type", y="counts") +
# TEMP !   	scale_y_continuous( labels = percent_format() )
# TEMP !
# TEMP !   save_plots(self, plot_id, ss, output_dir, output_formats)
# TEMP ! })


# t <- cast(f, sample_source ~ acc_chem_type_name, value="acc_chem_type_count")
# save_tables(self, t, plot_id, sample_sources, output_dir, output_formats, caption="Number of HBonds by Acceptor Type", caption.placement="top")
#
# sele <-"
# SELECT
#   don_site.HBChemType AS don_chem_type,
#   count(don_site.HBChemType) AS don_chem_type_count
# FROM
#   hbonds AS hbond,
#   hbond_sites AS don_site
# WHERE
#   hbond.struct_id == don_site.struct_id AND hbond.don_id == don_site.site_id
# GROUP BY
#   don_site.HBChemType;"
#
# f <- query_sample_sources(sample_sources, sele)
#
# f$don_chem_type_name <- don_chem_type_name_linear(f$don_chem_type)
# f <- na.omit(f, method="r")
#
# plot_id <- "hbonds_num_by_don_chem_type"
# ggplot(f, aes(don_chem_type_name, log(don_chem_type_count))) + theme_bw() +
#   geom_line(aes(colour=sample_source, group=sample_source), size=.2) +
# 	geom_point(aes(colour=sample_source), size=3) +
#   ggtitle("Number of HBonds by Donor Type") +
#   labs(x="Donor Type", y="log(Counts)")
# save_plots(self, plot_id, sample_sources, output_dir, output_formats)
#
# t <- cast(f, sample_source ~ don_chem_type_name, value="don_chem_type_count")
# save_tables(self, t, plot_id, sample_sources, output_dir, output_formats, caption="Number of HBonds by Donor Type", caption.placement="top")
#
#
# sele <-"
# SELECT
# 	don_counts.num_don_sites,
# 	acc_counts.num_acc_sites,
#   bond_counts.don_chem_type,
#   bond_counts.acc_chem_type,
# 	bond_counts.bond_counts
# FROM
# 	(SELECT
# 	  don_site.HBChemType AS don_chem_type,
# 	  acc_site.HBChemType AS acc_chem_type,
# 	  count(*) AS bond_counts
# 	FROM
# 	  hbonds AS hbond,
# 	  hbond_sites AS don_site,
# 	  hbond_sites AS acc_site
# 	WHERE
# 	  hbond.struct_id == don_site.struct_id AND hbond.don_id == don_site.site_id AND
# 	  hbond.struct_id == acc_site.struct_id AND hbond.acc_id == acc_site.site_id
# 	GROUP BY
# 	  don_site.HBChemType,
# 	  acc_site.HBChemType) AS bond_counts,
# 	(SELECT
# 		don.HBChemType AS don_chem_type,
# 		count(*) AS num_don_sites
# 	FROM
# 		hbond_sites AS don
# 	GROUP BY
# 		don_chem_type) AS don_counts,
# 	(SELECT
# 		acc.HBChemType AS acc_chem_type,
# 		count(*) AS num_acc_sites
# 	FROM
# 		hbond_sites AS acc
# 	GROUP BY
# 		acc_chem_type) AS acc_counts
# WHERE
# 	bond_counts.don_chem_type = don_counts.don_chem_type AND
# 	bond_counts.acc_chem_type = acc_counts.acc_chem_type;"
#
# f <- query_sample_sources(sample_sources, sele)
#
# f$don_chem_type_name <- don_chem_type_name_linear(f$don_chem_type)
# f$acc_chem_type_name <- acc_chem_type_name_linear(f$acc_chem_type)
# f <- na.omit(f, method="r")
#
# plot_id <- "hbonds_num_don_acc_by_chem_type"
# ggplot(f, aes(don_chem_type_name, log(bond_counts))) + theme_bw() +
#   geom_line(aes(colour=sample_source, group=sample_source), size=.2) +
# 	geom_point(aes(colour=sample_source), size=2) +
# 	facet_wrap( ~ acc_chem_type_name) +
#   ggtitle("Number of HBonds by Donor and Acceptor Type") +
#   labs(x="Donor Type", y="log(Counts)") +
# 	coord_flip()
#
# save_plots(self, plot_id, sample_sources, output_dir, output_formats)
#
# plot_id <- "hbonds_num_acc_don_by_chem_type"
# ggplot(f, aes(acc_chem_type_name, log(bond_counts))) + theme_bw() +
#   geom_line(aes(colour=sample_source, group=sample_source), size=.2) +
# 	geom_point(aes(colour=sample_source), size=2) +
# 	facet_wrap( ~ don_chem_type_name) +
#   ggtitle("Number of HBonds by Acceptor and Donor Type") +
#   labs(x="Acceptor Type", y="log(Counts)") +
# 	coord_flip() +
#
# save_plots(self, plot_id, sample_sources, output_dir, output_formats)
#
# table_id <- "hbond_num_bonds_by_don_acc_chem_type"
# t <- cast(f,
# 	sample_source +
# 	acc_chem_type_name ~ don_chem_type_name,
# 	value="bond_counts")
# save_tables(self, t, plot_id, sample_sources, output_dir, output_formats, caption="Number of HBonds by Acceptor and Donor Type", caption.placement="top")
#

})) # end FeaturesAnalysis
