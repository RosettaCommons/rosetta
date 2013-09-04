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
id = "tau_angles",
filename = "scripts/analysis/plots/loops/tau_angles.R",
author = "Brian D. Weitzner",
brief_description = "",
long_description = "
This features analysis scripts performs a expectation maximization (EM) mixture
model of two Gaussians to describe the tau angle in CDR H3 loops.
",

feature_reporter_dependencies = c("loop_anchor_features"),
run=function(self, sample_sources, output_dir, output_formats){

library(ggthemes)
library(mixtools)

# helper function to format numbers nicely for plot annotations
disp <- function(x) {
  formatC(signif(x, digits=3), digits=3, format="fg", flag="#")
}

# Loop anchor transfrom parameters for all candidate loops of a given length
sele <-paste("SELECT omega FROM loop_anchor_transforms_three_res;")

f <- query_sample_sources(sample_sources, sele)

# because we are querying a db that uses "omega" in place of "tau"
# the variable is named omega, but it refers to tau.  I know, I know...
expectation_maximization.model <- normalmixEM(x=f$omega)
parameters <- expectation_maximization.model[c("lambda", "mu", "sigma")]

# add Gaussians and mixture model to data.frame
f$curve.1 <- parameters$lambda[1] * dnorm(f$omega, parameters$mu[1], 
                                          parameters$sigma[1])
f$curve.2 <- parameters$lambda[2] * dnorm(f$omega, parameters$mu[2], 
                                          parameters$sigma[2])
f$combined <- f$curve.1 + f$curve.2

p <- ggplot(f, aes(x=omega)) +
  geom_histogram(aes(y=..density..), fill="lightgrey", colour="white") + 
  geom_hline(yintercept=1, colour="white") + 
  geom_hline(yintercept=2, colour="white") + 
  geom_hline(yintercept=3, colour="white") + 
  geom_line(aes(y=curve.1), size=1, colour="darkgrey") + 
  geom_line(aes(y=curve.2), size=1, colour="darkgrey") + 
  geom_line(aes(y=combined), size=1) + 
  ggtitle(expression(paste(tau, " in antibodies"))) +
  scale_x_continuous(expression(paste(tau, " (radians)")), expand=c(0, 0)) +
  scale_y_continuous("Density") +
  annotate("text", family="serif", parse=TRUE, x=2.2, y=0.9, hjust=0,
           label=paste(expression(lambda[1]),": ", 
                       disp(parameters$lambda[1]))) +
  annotate("text", family="serif", parse=TRUE, x=2.2, y=0.8, hjust=0,
           label=paste(expression(mu[1]),": ", disp(parameters$mu[1]))) +
  annotate("text", family="serif", parse=TRUE, x=2.2, y=0.7, hjust=0,
           label=paste(expression(sigma[1]),": ", 
                       disp(parameters$sigma[1]))) +
  annotate("text", family="serif", parse=TRUE, x=2.35, y=0.9, hjust=0,
           label=paste(expression(lambda[2]),": ", 
                       disp(parameters$lambda[2]))) +
  annotate("text", family="serif", parse=TRUE, x=2.35, y=0.8, hjust=0,
           label=paste(expression(mu[2]),": ", disp(parameters$mu[2]))) +
  annotate("text", family="serif", parse=TRUE, x=2.35, y=0.7, hjust=0,
           label=paste(expression(sigma[2]),": ", 
                       disp(parameters$sigma[2]))) +
  annotate("text", family="serif", parse=TRUE, x=2.2, y=1.2, hjust=0,
           label=paste(expression(f(tau) == sum(lambda[i] * frac(1, sigma[i] * sqrt(2*pi)) ~~ e^{- ~~ frac((tau-mu[i])^2, 2*sigma[i]^2)}, i==1, n)))) +
  theme_tufte()

save_plots(self, paste("tau_density_estimate", sep = "_"), sample_sources, 
           output_dir, output_formats, fonts=c("serif"))

})) # end FeaturesAnalysis