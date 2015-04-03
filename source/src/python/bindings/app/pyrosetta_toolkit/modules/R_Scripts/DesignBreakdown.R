#!/usr/bin/Rscript

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/modules/R_Scripts/BreakdownDesign.R
## @brief  R Script for plotting design breakdown for modules/DesignBreakdown.py
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Args are: Database_In, Output_directory

#Database Columns Required: 
 #region, type (design or reference), rosetta_position, pdb_position(ex: 24L), aa, prob

#Requirements:
require(RSQLite)
require(rgl)


args=(commandArgs(TRUE))
if(length(args)==0){
    print("No arguments supplied.")}

if (length(args)!=2){
	print("Output Directory Required")}


drv <- dbDriver("SQLite")

con <- dbConnect(drv, dbname = toString(args[1]))
outdir = toString(args[2])

sql <- "SELECT DISTINCT aa FROM design_data WHERE type='design'"
aminos = dbGetQuery(con, sql)

sql <- "SELECT DISTINCT region FROM design_data"
regions = dbGetQuery(con, sql)

#sink(paste(sep="", outdir, "/", "LOG.txt"))
for (r in regions$region){
	f = paste(sep="", outdir,"/", r, "_PLOTS.pdf")
	if (regexpr(':', r)){
		sp = strsplit(r, ':')
		new = paste(sep="_", sp[[1]][1], sp[[1]][2], sp[[1]][3])
		f = paste(sep="", outdir,"/", new, "_PLOTS.pdf")
	}
	
	pdf(file=f, width=7*5/2, height=7*4/2)
	par(mfrow=c(4,3))


	sql <-"SELECT DISTINCT pdb_position FROM design_data WHERE region=?"
	pdb_positions = dbGetQuery(con, sql, r)

	for (position in pdb_positions$pdb_position){

		sql<-"SELECT prob FROM design_data WHERE type='design' and region=? and pdb_position=?"
		df = data.frame(region=r, pdb_position=position)
		prob = dbGetQuery(con, sql, df)

		sql<-"SELECT aa FROM design_data WHERE type='reference' and region=? and pdb_position=?"
		reference_aa = dbGetQuery(con, sql, df)

		aa = reference_aa$aa[1]
		
		#I think blue looks great.  If you want something else, change it here.
		barplot(prob$prob, names.arg=aminos$aa, density=prob$prob*100, col='blue')
		title(main=position, xlab=paste('AA (reference: ',aa, ')'), ylab='Percent Observed', font.main=2, cex.main=1.6, cex.lab=1.2)
		#legend('topright', legend = c(paste('reference:', aa)), box.col='blue', text.col='blue')
		
		

	}
	dev.off()
}
warnings()
#sink()
#close(con)


