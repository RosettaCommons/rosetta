# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file test/scientific/rotamer_recovery_analysis.R


load_it <- function(pkg_name){
  loaded_ok <- require(pkg_name, character.only=TRUE)
  if(!loaded_ok){
    stop(paste("Package '", pkg_name,"' is required but not found. To install R packages, open R and run:\n    install.packages('", pkg_name, "')", sep=""))
  }
}

sapply(c("plyr", "reshape", "ggplot2", "RSQLite"), load_it)


########## I/O files and Parameters #############
rotamer_recovery_results_db_fname="outputs/rotamer_recovery_results.db3"

#################################################

uniform_weight <- function(x){ rep(1/length(x),length(x)) }
estimate_density_1d <-function(df, ids, variable, weight_fun=uniform_weight, min_count=20){
	n_pts <- 200
	xlim <- range(df[,variable])
	compute_density <- function(factor_df){
		if (nrow(factor_df) < min_count){
			return( data.frame(x=seq(xlim[1], xlim[2], n_pts), y=0))
		} else {
			weights <- weight_fun(factor_df[,variable])
			d <- density(x=factor_df[,variable], from=xlim[1], to=xlim[2], n=n_pts,
									 weights=weights)
			return(data.frame(x=d$x, y=d$y))
		}
	}
	ddply(df, ids, compute_density)
}


engine <-SQLite()
con <- dbConnect(engine, rotamer_recovery_results_db_fname)
z <- dbSendQuery(con, "create index if not exists rotamer_recovery_index on rotamer_recovery ( struct1_name, chain1, res1 );")

sele <- "
SELECT score as score,
       aa as aa
FROM   rotamer_recovery"

recovery_data <- dbGetQuery(con, sele)

dens <- estimate_density_1d(
	recovery_data,
	ids=c("aa"),
	variable=c("score"))

p <- ggplot(dens, aes(x=x, y=log(y+1)))
p <- p + geom_line()
p <- p + facet_wrap(~aa, ncol=5)
p <- p + opts(title="Rotamer Recovery by Amino Acid Type")
p <- p + labs(x="Automorphic RMSD (lower is better)")
p <- p + labs(y="log(Density + 1)")

if( capabilities()['png'] ){
	png("outputs/rotamer_recovery_by_amino_acid.png")
	print(p)
	dev.off()
} else {
	print("Unable to generate png output because the device has not been installed into R")
}

ggsave(
	filename="outputs/rotamer_recovery_by_amino_acid.pdf",
	plot=p)

