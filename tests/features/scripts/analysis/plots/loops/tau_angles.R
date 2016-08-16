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
library(movMF)

# helper function to format numbers nicely for plot annotations
disp <- function(x) {
  formatC(signif(x, digits=3), digits=3, format="fg", flag="#")
}

model_tau101 <- function(df) {
  em.model <- normalmixEM(x=df$tau101)
  parameters <- em.model[c("lambda", "mu", "sigma")]
  
  # add Gaussians and mixture model to data.frame
  df$curve.1 <- parameters$lambda[1] * dnorm(df$tau101, parameters$mu[1], 
                                            parameters$sigma[1])
  df$curve.2 <- parameters$lambda[2] * dnorm(df$tau101, parameters$mu[2], 
                                            parameters$sigma[2])
  df$combined <- df$curve.1 + df$curve.2
  
  hline.data <- data.frame(z = c(0.02, 0.04, 0.06))
  p <- ggplot(df, aes(x=tau101)) +
    geom_histogram(aes(y=..density..), fill="lightgrey", colour="white") +
    geom_hline(aes(yintercept = z), hline.data, colour="white") +
    geom_line(aes(y=curve.1), size=1, colour="darkgrey") +
    geom_line(aes(y=curve.2), size=1, colour="darkgrey") + 
    geom_line(aes(y=combined), size=1) +
    scale_x_continuous(expression(paste(tau[101], " (degrees)")), expand=c(0, 0)) +
    scale_y_continuous("Density") +
    annotate("text", family="sans", parse=TRUE, x=120, y=0.023, hjust=0,
             label=paste(expression(lambda[1]),": ", 
                         disp(parameters$lambda[1]))) +
    annotate("text", family="sans", parse=TRUE, x=120, y=0.019, hjust=0,
             label=paste(expression(mu[1]),": ", disp(parameters$mu[1]))) +
    annotate("text", family="sans", parse=TRUE, x=120, y=0.015, hjust=0,
             label=paste(expression(sigma[1]),": ", 
                         disp(parameters$sigma[1]))) +
    annotate("text", family="sans", parse=TRUE, x=133, y=0.023, hjust=0,
             label=paste(expression(lambda[2]),": ", 
                         disp(parameters$lambda[2]))) +
    annotate("text", family="sans", parse=TRUE, x=133, y=0.019, hjust=0,
             label=paste(expression(mu[2]),": ", disp(parameters$mu[2]))) +
    annotate("text", family="sans", parse=TRUE, x=133, y=0.015, hjust=0,
             label=paste(expression(sigma[2]),": ", 
                         disp(parameters$sigma[2]))) +
    annotate("text", family="sans", parse=TRUE, x=120, y=0.031, hjust=0,
             label=paste(expression(italic(f)(tau[101]) == sum(lambda[italic(i)] 
                * frac(1, sigma[italic(i)] * sqrt(2*pi)) ~ italic(e)^{- ~ 
                frac((tau-mu[italic(i)])^2, 2*sigma[italic(i)]^2)}, 
                italic(i)==1, italic(n))))) +
    theme_tufte(base_family="sans")
  
  save_plots(self, paste("tau_density_estimate", sep = "_"), sample_sources, 
             output_dir, output_formats,
             fonts=c("sans"))
  
  # return parameters
  bigger_curve <- ifelse(parameters$lambda[1] > parameters$lambda[2], 1, 2)
  return(list(mu=parameters$mu[bigger_curve], 
              sigma=parameters$sigma[bigger_curve]))
}

model_alpha101 <- function(df) {
  pts_on_unit_circle <- cbind(cos(df$alpha101 * pi / 180), 
                              sin(df$alpha101 * pi / 180))
  
  em.model <- movMF(pts_on_unit_circle, 2)
  
  # compute mean angles (in degrees) of distributions
  mu <- atan2(em.model$theta[,2], em.model$theta[,1]) * 180 / pi
  
  # compute the concentration parameter
  kappa <- sqrt(rowSums(em.model$theta^2))
  
  
  # Because the data are laid out on a unit circle and we will be plotting
  # on a domain of 360, divide the densities by 360.
  df$curve.1 <- em.model$alpha[1] * dmovMF(pts_on_unit_circle, 
                                           em.model$theta[1,]) / 360
  df$curve.2 <- em.model$alpha[2] * dmovMF(pts_on_unit_circle, 
                                           em.model$theta[2,]) / 360
  df$combined <- df$curve.1 + df$curve.2
  
  hline.data <- data.frame(z = seq(0.005, 0.025, by=0.005))
  p <- ggplot(df, aes(x=alpha101)) +
    geom_histogram(aes(y=..density..), fill="lightgrey", colour="white") +
    geom_hline(aes(yintercept = z), hline.data, colour="white") +
    geom_line(aes(y=curve.1), size=1, colour="grey50") +
    geom_line(aes(y=curve.2), size=1, colour="darkgrey") + 
    scale_x_continuous(expression(paste(alpha[101], " (degrees)")), 
                       limits=c(-180, 180), breaks=seq(-180,180, by=45)) +
    scale_y_continuous("Density", expand=c(0, 0)) + 
    annotate("text", family="sans", parse=TRUE, x=-150, y=0.013, hjust=0,
             label=paste(expression(lambda[1]),": ", 
                         disp(em.model$alpha[1]))) +
    annotate("text", family="sans", parse=TRUE, x=-150, y=0.011, hjust=0,
             label=paste(expression(mu[1]),": ", disp(mu[1]))) +
    annotate("text", family="sans", parse=TRUE, x=-150, y=0.009, hjust=0,
             label=paste(expression(kappa[1]),": ", 
                         disp(kappa[1]))) +
    annotate("text", family="sans", parse=TRUE, x=-60, y=0.013, hjust=0,
             label=paste(expression(lambda[2]),": ", 
                         disp(em.model$alpha[2]))) +
    annotate("text", family="sans", parse=TRUE, x=-60, y=0.011, hjust=0,
             label=paste(expression(mu[2]),": ", disp(mu[2]))) +
    annotate("text", family="sans", parse=TRUE, x=-60, y=0.009, hjust=0,
             label=paste(expression(kappa[2]),": ", 
                         disp(kappa[2]))) +
    annotate("text", family="sans", parse=TRUE, x=-150, y=0.017, hjust=0,
             label=paste(expression(italic(f)(alpha[101]) == 
                         sum(lambda[italic(i)] ~ frac(italic(e)^{kappa~cos(frac(
                         (alpha[101]-mu[italic(i)])~pi, 180))}, 2*pi ~ 
                         italic(I)[0](kappa[italic(i)])), italic(i)==1, 
                         italic(n))))) +
    theme_tufte(base_family="sans")
  
  save_plots(self, paste("alpha_density_estimate", sep = "_"), sample_sources, 
             output_dir, output_formats, fonts=c("sans"))
  
  # return parameters
  bigger_curve <- ifelse(em.model$alpha[1] > em.model$alpha[2], 1, 2)
  return(list(mu=mu[bigger_curve], kappa=kappa[bigger_curve]))
}

tau_vs_alpha <- function(df, tau, alpha, n_std_dev=3.0) {
  # kappa closely approximates sigma (in radians): (kappa)^-1/2 ~ sigma
  alpha$sigma <- (alpha$kappa)^-(1/2) * 180 / pi
  
  # horizontal box value, vertical box value
  ymin <- rbind(tau$mu - (n_std_dev * tau$sigma), -Inf)
  ymax <- rbind(tau$mu + (n_std_dev * tau$sigma), Inf)
  xmax <- rbind(Inf, alpha$mu + (n_std_dev * alpha$sigma))
  xmin <- rbind(-Inf, alpha$mu - (n_std_dev * alpha$sigma))
  boxes <- data.frame(xmin=as.numeric(xmin), xmax=as.numeric(xmax), 
                      ymin=as.numeric(ymin), ymax=as.numeric(ymax))
  
  hline.data <- data.frame(z = c(boxes$ymin[1], boxes$ymax[1]))
  vline.data <- data.frame(z = c(boxes$xmin[2], boxes$xmax[2]))
  p <- ggplot(df, aes(x=alpha101, y=tau101, colour=num.hbonds)) + 
    geom_rect(aes(NULL, NULL, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), 
              data=boxes, fill="lightgrey", colour=NA) +
    geom_hline(aes(yintercept = z), hline.data, colour="white", size=1) +
    geom_vline(aes(xintercept = z), vline.data, colour="white", size=1) +
    geom_point() +
    scale_x_continuous(expression(paste(alpha[101], " (Degrees)"))) +
    scale_y_continuous(expression(paste(tau[101], " (Degrees)"))) +
    scale_colour_discrete(name = "Number of stabilizing Hydrogen bonds") +
    theme_tufte(base_family="sans") +
    theme(legend.position="bottom") +
    geom_rangeframe()
  
  # use cairo to preserve transparency in EPS format
  # NOTE: cairo fucks up Greek glyphs, so they will need to be added back in
  # Illustrator or some other post-processing application.
#   save_plots(self, paste("alpha_vs_tau", sep = "_"), 
#              sample_sources, output_dir, 
#              output_formats[output_formats$extension == ".eps",], 
#              device=cairo_ps, family=c("sans"))
  
  # don't use cairo for non-EPS
  save_plots(self, paste("alpha_vs_tau", sep = "_"), 
             sample_sources, output_dir, 
             output_formats, 
             fonts=c("sans"))
}
# Loop anchor transfrom parameters for all candidate loops of a given length
sele <-paste("SELECT struct_id, tau101, alpha101 FROM loop_anchor_transforms")

f <- query_sample_sources(sample_sources, sele)
tau.params <- model_tau101(f)
alpha.params <- model_alpha101(f)

# There is a seemingly strange off-by-one assertion in the struct_id for the lat
# and hbond tables. This is due to the fact that for the set of known 
# antibodies the FeatureReporter is restricted to the H3 loop.  This makes in-
# cluding Hbonds with residues not in the H3 loop impossible.  To get around
# this, I add a second FeatureReporter mover that can compute Hbonds for the
# entire structure in my XML file, which has the effect of struct_id 1 corr-
# esponding to the restricted, H3-only calculation and struct_id 2 corresponding
# to the whole structure calculation.  The off-by-one assertion in the WHERE
# clause allows the H3-only information (that LAT parameters) to be joined with
# the whole strucutre information (all Hbonds).
hbond_query <- "
SELECT
  lat.struct_id, lat.residue_begin, lat.residue_end,
(lat.residue_end - lat.residue_begin + 1) AS length,
don.resNum as don, acc.resNum as acc,
(don.resNum - lat.residue_begin + 1) AS hbond_don,
(acc.resNum - lat.residue_begin + 1) AS hbond_acc,
CASE
WHEN don.HBChemType == 'hbdon_PBA' AND acc.HBChemType == 'hbacc_PBA' THEN 'BB-BB'
WHEN don.HBChemType != 'hbdon_PBA' AND acc.HBChemType != 'hbacc_PBA' THEN 'SC-SC'
ELSE 'BB-SC' END AS hbType,
structures.tag

FROM
loop_anchor_transforms AS lat, hbonds AS hb,
hbond_sites AS don, hbond_sites AS acc, structures
WHERE
lat.struct_id = hb.struct_id - 1 AND lat.struct_id = structures.struct_id AND
don.struct_id = hb.struct_id AND don.site_id =hb.don_id AND
acc.struct_id = hb.struct_id AND acc.site_id =hb.acc_id AND
((hbond_don = 4 AND hbond_acc = length - 3 AND don.HBChemType = 'hbdon_PBA' AND acc.HBChemType = 'hbacc_PBA')
OR
(hbond_don = length - 1 AND hbond_acc = 2 AND don.HBChemType = 'hbdon_PBA' AND acc.HBChemType = 'hbacc_PBA')
OR
(hbond_don = length + 1 AND hbond_acc = length - 2 AND don.HBChemType != 'hbdon_PBA' AND acc.HBChemType = 'hbacc_PBA')	OR
(hbond_don = 2 AND hbond_acc = length - 1 AND don.HBChemType != 'hbdon_PBA' AND acc.HBChemType != 'hbacc_PBA'))
GROUP BY lat.struct_id, don, acc, hbType
"

hbq <- query_sample_sources(sample_sources, hbond_query)

hbdata <- ddply(hbq, .(struct_id, residue_begin, residue_end), function(df) {
  data.frame(struct_id=unique(df$struct_id), 
             num.hbonds=nrow(df))
})

f$num.hbonds <- 0
for (struct_id in f$struct_id) { 
  if (struct_id %in% hbdata$struct_id) {
    f[f$struct_id == struct_id,]$num.hbonds <- hbdata[hbdata$struct_id == struct_id,]$num.hbonds
  }
}
f$num.hbonds <- as.factor(f$num.hbonds)

tau_vs_alpha(f, tau.params, alpha.params)


})) # end FeaturesAnalysis