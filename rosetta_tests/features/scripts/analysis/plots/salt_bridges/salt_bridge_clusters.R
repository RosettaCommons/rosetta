
# classes and methods for spatial data
library(sp)

# label points that fall into clusters of salt bridge interactions
# input: data.frame with columns psi_degrees and rho
# output: input data.frame with extra column 
add_salt_bridge_clusters <- function(f){
  	if(!("psi_degrees" %in% names(f)) || !("rho" %in% names(f))){
          stop("add_salt_bridge_clusters requires a data.frame with the columns \"psi_degrees\" and \"rho\"")
        }
	clusters <-rbind(
		data.frame(
			x=c(-158, -118, -118, -138, -158, -158),
			y=c(4, 4, 3.58, 3.25, 3.58, 4),
			cluster_label=factor("bidentite_a1")),
		data.frame(
			x=c(-118, -78, -78, -98, -118, -118),
			y=c(4, 4, 3.58, 3.25, 3.58, 4),
			cluster_label=factor("bidentite_a2")),
		data.frame(
			x=c(-118, -98, -98, -138, -138, -118),
			y=c(3.58, 3.25, 3.1, 3.1, 3.25, 3.58),
			cluster_label=factor("bifurcated_a")),
	        data.frame(
			x=c(-40, 0, 0, -20, -40, -40),
			y=c(4, 4, 3.58, 3.25, 3.58, 4),
			cluster_label=factor("bidentite_b1")),
		data.frame(
			x=c(0, 40, 40, 20, 0, 0),
			y=c(4, 4, 3.58, 3.25, 3.58, 4),
			cluster_label=factor("bidentite_b2")),
		data.frame(
			x=c(0, 20, 20, -20, -20, 0),
			y=c(3.58, 3.25, 3.1, 3.1, 3.25, 3.58),
			cluster_label=factor("bifurcated_b")),
		data.frame(
			x=c(70, 120, 120, 110, 80, 70, 70),
			y=c(4.1, 4.1, 3.6, 3.45, 3.6, 3.8, 4.1),
			cluster_label=factor("GDH")))
	
	# spatial points
	cl <- rep(NA, nrow(f))
	points <- SpatialPoints(f[,c("psi_degrees", "rho")])
	d_ply(clusters, .(cluster_label), function(cluster){
		cluster_label <- as.character(cluster[1, "cluster_label"])
		print(cluster_label)
		polygon <- SpatialPolygons(list(Polygons(list(Polygon(cluster[,c("x", "y")])), ID=cluster_label)))
		cl[!is.na(over(points, polygon))] <<- cluster_label
	})
	factor(cl)
}
