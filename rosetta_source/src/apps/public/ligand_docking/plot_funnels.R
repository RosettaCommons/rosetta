# Just pull the data in from file
dock.load <- function(filename) {
    if(!file.exists(filename)) return(NULL)
	x <- read.table(filename, header=TRUE)
	x$score_rank <- rank(x$total_score, ties='min') / nrow(x)
	if(!is.null(x$native_ligand_auto_rms_no_super)) { x$ligand_auto_rms_no_super <- x$native_ligand_auto_rms_no_super }
	if(!is.null(x$native_ligand_auto_rms_with_super)) { x$ligand_auto_rms_with_super <- x$native_ligand_auto_rms_with_super }
	invisible(x)
}

# Once data is loaded, select the best runs and sort them.
# This is the original filter by total score, sort by interface delta.
dock.filter <- function(x, fraction=0.05) {
	if(is.null(x)) return(invisible(x))
	x <- x[order(-x$ligand_is_touching, x$total_score),][1:round(fraction*dim(x)[1]),]
    x <- x[order(x$interface_delta),]
    invisible(x)
}

pdf(file="funnel_plots.pdf", width=11, height=8.5)
for(f in list.files(pattern=glob2rx("*.tab"))) {
	basename <- substr(f, 1, nchar(f)-4)
	print(basename)
	#pdf(file=sprintf("%s.pdf", basename), width=11, height=8.5)
	x <- dock.load(f)
	x1 <- dock.filter(x)
	x2 <- dock.filter(x, 1.0) # just to sort
	max.ifd <- max(x1$interface_delta)
	x2 <- x2[ x2$interface_delta <= max.ifd, ]
	z <- dock.filter(dock.load(sprintf("minnat/%s", f)), 1.0) # may be null, won't hurt anything
	z <- z[ z$interface_delta <= max.ifd, ]
	rx <- range(0, x2$ligand_auto_rms_no_super)
	ry <- range(x2$interface_delta, z$interface_delta)
	plot(x2$ligand_auto_rms_no_super, x2$interface_delta, xlim=rx, ylim=ry, main=basename, xlab='RMSD to native (Å)', ylab='interaction energy', col='#ccccff')
	abline(v=2, col='#cccccc')
	if(!is.null(z)) points(z$ligand_auto_rms_no_super, z$interface_delta, col='#33cc33', pch=16)
	points(x1$ligand_auto_rms_no_super, x1$interface_delta, col='#3333cc')
	#dev.off()
}
dev.off()