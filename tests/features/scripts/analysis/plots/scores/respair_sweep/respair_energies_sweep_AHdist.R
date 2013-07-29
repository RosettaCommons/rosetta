# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

check_setup()
feature_analyses <- c(feature_analyses, new("FeaturesAnalysis",
id = "respair_energies_sweep_AHdist",
author = "Andrew Leaver-Fay",
brief_description = "This is a script to plot the inter-residue energies achieved by scanning through distances separating an acceptor on one residue and a donor hydrogen on another residue.  The schema used for the database in this defined in the pilot app sweep_respair_energies",
#feature_reporter_dependencies = c(""),
run=function(self, sample_sources, output_dir, output_formats){

sele <-"
SELECT
  respair.res1name AS acc_rsd,
  respair.res2name AS don_rsd,
	geom.geom_value AS AHdist,
  score.score_name,
	score.score_value
FROM
  residue_pairs as respair,
  conformations as conf,
  conf_scores as score,
  geometries as geom
WHERE
	respair.focused_geom_param = 'sweep_AHdist' AND
  respair.respair_id == conf.respair_id AND
  score.conf_id == conf.conf_id AND
  geom.conf_id == conf.conf_id AND
	geom.geom_name == 'AHdist'
ORDER BY
  geom.geom_value;"

f <- query_sample_sources(sample_sources, sele)
f <- na.omit(f, method="r")

# take the row with the lowest energy
ftot <- subset( f, subset=score_name=="total" )
minvals <- ddply( ftot, acc_rsd ~ don_rsd, function(dat)dat[order(dat$score_value,decreasing=FALSE)[1],] )
print(summary(minvals))

ylimits <- c( min(minvals$score_value), 1 )

plot_id <- "sweep_AHdist_respair_scores_by_component"

ss_id <- as.character( sample_sources[1,"sample_source"] )
print(paste("Testing! ",ss_id))

p <- ggplot(data=f) + theme_bw() +
	#geom_path(aes(x=AHdist, y=score_value, colour=score_name)) +
	geom_path(data=subset(f,subset=score_name != "total"), aes(x=AHdist, y=score_value, colour=score_name)) +
	geom_path(data=subset(f,subset=score_name == "total"), aes(x=AHdist, y=score_value, colour=I("total")),size=1.2) +
	#geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
  geom_indicator(data=minvals, aes(indicator=AHdist), xpos="right", ypos="bottom",group=1) +
	facet_grid(acc_rsd ~ don_rsd) +
	opts(title = paste("Residue pair energy as a function of Acceptor-Hydrogen distance\nsample source:", ss_id, " ")) +
	scale_y_continuous("Weighted Energies", limits=ylimits, breaks=c(-3,-2,-1,0,1)) +
	scale_x_continuous(expression(paste('Acceptor -- Proton Distance (', ring(A), ')')), limits=c(1.5,3.5), breaks=c(2.0,2.5,3.0))

if(nrow(sample_sources) <= 3){
	p <- p + opts(legend.position="bottom", legend.direction="horizontal")
}
save_plots(self, plot_id, sample_sources, output_dir, output_formats)


# now focus in on the range between 1.5 and 2.1 to look closely at the energies in the hbonds
f2 <- subset( f, subset=(AHdist >= 1.5 & AHdist <= 2.1 ))
f2 <- subset( f2, subset=(score_name == "hbond_lr_bb" | score_name == "hbond_sr_bb" | score_name == "hbond_bb_sc" | score_name == "hbond_sc" | score_name == "hack_elec" | score_name == "total" | score_name == "fa_rep" ))
ftot <- subset( f2, subset=score_name=="total" )
minvals <- ddply( ftot, acc_rsd ~ don_rsd, function(dat)dat[order(dat$score_value,decreasing=FALSE)[1],] )

ylimits <- c( min(minvals$score_value), 1 )

plot_id <- "sweep_AHdist_respair_scores_1p5_to_2p1A"
p <- ggplot(data=f2) + theme_bw() +
	geom_path(data=subset(f2,subset=score_name != "total"), aes(x=AHdist, y=score_value, colour=score_name)) +
	geom_path(data=subset(f2,subset=score_name == "total"), aes(x=AHdist, y=score_value, colour=I("total")),size=1.2) +
	#geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
  geom_indicator(data=minvals, aes(indicator=AHdist), xpos="right", ypos="bottom",group=1) +
	facet_grid(acc_rsd ~ don_rsd) +
	opts(title = paste("Residue pair energy as a function of Acceptor-Hydrogen distance\nsample source:", ss_id, " ")) +
	scale_y_continuous("Weighted Energies", limits=ylimits, breaks=c(-3,-2,-1,0,1)) +
	scale_x_continuous(expression(paste('Acceptor -- Proton Distance (', ring(A), ')')), limits=c(1.5,2.1), breaks=c(1.6,1.7,1.8,1.9,2.0))

if(nrow(sample_sources) <= 3){
	p <- p + opts(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)

# now zoom out to the range between 0 and 4
f3 <- subset( f, subset=(AHdist >= 0 & AHdist <= 4 ))
f3 <- subset( f3, subset=(score_name == "hbond_lr_bb" | score_name == "hbond_sr_bb" | score_name == "hbond_bb_sc" | score_name == "hbond_sc" | score_name == "hack_elec" | score_name == "total" | score_name == "fa_rep" ))
ftot <- subset( f3, subset=score_name=="total" )
minvals <- ddply( ftot, acc_rsd ~ don_rsd, function(dat)dat[order(dat$score_value,decreasing=FALSE)[1],] )

ylimits <- c( min(minvals$score_value), 1 )

plot_id <- "sweep_AHdist_respair_scores_0_to_4"
p <- ggplot(data=f3) + theme_bw() +
	geom_path(data=subset(f3,subset=score_name != "total"), aes(x=AHdist, y=score_value, colour=score_name)) +
	geom_path(data=subset(f3,subset=score_name == "total"), aes(x=AHdist, y=score_value, colour=I("total")),size=1.2) +
	#geom_indicator(aes(indicator=counts, colour=sample_source, group=sample_source)) +
  geom_indicator(data=minvals, aes(indicator=AHdist), xpos="right", ypos="bottom",group=1) +
	facet_grid(acc_rsd ~ don_rsd) +
	opts(title = paste("Residue pair energy as a function of Acceptor-Hydrogen distance\nsample source:", ss_id, " ")) +
	scale_y_continuous("Weighted Energies", limits=ylimits, breaks=c(-3,-2,-1,0,1)) +
	scale_x_continuous(expression(paste('Acceptor -- Proton Distance (', ring(A), ')')), limits=c(0,4), breaks=c(1,2,3))

if(nrow(sample_sources) <= 3){
	p <- p + opts(legend.position="bottom", legend.direction="horizontal")
}

save_plots(self, plot_id, sample_sources, output_dir, output_formats)


})) # end FeaturesAnalysis
