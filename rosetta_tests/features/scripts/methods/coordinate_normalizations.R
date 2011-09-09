

radial_3d_normalization <- function(x){ 1/(x^2*sum(1/x^2))}


uniform_normalization <- function(x){ rep(1/length(x),length(x)) }


no_normalization <- function(x) rep(1,length(x))


# Equal Area Coordinate Grids 
major_long_coords <- transform(
	expand.grid(long=seq(0, 2*pi, length.out=200), lat=c(pi/6, pi/2)),
	capx = 2*sin(lat/2)*cos(long),
	capy = 2*sin(lat/2)*sin(long))

minor_long_coords <- transform(
	expand.grid(long=seq(0, 2*pi, length.out=200), lat=c(pi/3, pi*2/3)),
	capx = 2*sin(lat/2)*cos(long),
	capy = 2*sin(lat/2)*sin(long))

major_lat_coords <- transform(
	expand.grid(long=seq(0,2*pi, length.out=5), lat=c(pi/3, pi*2/3, length.out=200)),
	capx = 2*sin(lat/2)*cos(long),
	capy = 2*sin(lat/2)*sin(long))

minor_lat_coords <- transform(
	expand.grid(long=seq(0,2*pi, length.out=5) + pi/4, lat=c(pi/3, pi*2/3, length.out=200)),
	capx = 2*sin(lat/2)*cos(long),
	capy = 2*sin(lat/2)*sin(long))

polar_equal_area_grids_bw <- list(
	geom_path(data=minor_long_coords, aes(x=capx, y=capy, group=lat), size=.5, colour="grey98"),
	geom_path(data=minor_lat_coords, aes(x=capx, y=capy, group=long), size=.5, colour="grey98"),
	geom_path(data=major_long_coords, aes(x=capx, y=capy, group=lat), size=.2, colour="grey90"),
	geom_path(data=major_lat_coords, aes(x=capx, y=capy, group=long), size=.2, colour="grey90"),
	opts(panel.grid.major = theme_blank(), panel.grid.minor = theme_blank()))
