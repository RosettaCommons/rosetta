

scale_x_AHdist <- scale_x_continuous(
	expression(paste('Acceptor -- Hydrogen Distance (', ring(A), ')')),
	limits=c(1.4, 3), breaks=c(1.4, 1.8, 2.2, 2.6, 3))

scale_y_AHdist <- scale_y_continuous(
	expression(paste('Acceptor -- Hydrogen Distance (', ring(A), ')')),
	limits=c(1.4, 3), breaks=c(1.4, 1.8, 2.2, 2.6, 3))

scale_x_ADdist <- scale_x_continuous(
	expression(paste('Acceptor -- Donor Distance (', ring(A), ')')),
	limits=c(2.4, 3.3), breaks=c(2.6, 2.9, 3.2))

scale_y_ADdist <- scale_y_continuous(
	expression(paste('Acceptor -- Donor Distance (', ring(A), ')')),
	breaks=c(2.3, 2.8, 3.3))

scale_x_cosBAH <- scale_x_continuous(
	"cos(Base -- Acceptor -- Hydrogen)",
	limit=c(-.3,1), breaks=c(-.2, .2, .6, 1))

scale_y_cosBAH <- scale_y_continuous(
	"cos(Base -- Acceptor -- Hydrogen)",
	limit=c(-.3,1), breaks=c(-.2, .2, .6, 1))

scale_x_cosAHD <- scale_x_continuous(
	"cos(Acceptor -- Hydrogen -- Donor)",
	limit=c(0,1), breaks=c(.2, .4, .6, .8, 1))

scale_y_cosAHD <- scale_y_continuous(
	"cos(Acceptor -- Hydrogen -- Donor)",
	limit=c(0,1), breaks=c(.2, .4, .6, .8, 1))

scale_x_chi <- scale_x_continuous(
	"Base -- Acceptor Torsion (Radians)",
	limit=c(-pi,pi), breaks=c(-pi, -pi*2/3, -pi/3, 0, pi/3, pi*2/3, pi))

scale_y_chi <- scale_y_continuous(
	"Base -- Acceptor Torsion (Radians)",
	limit=c(-pi,pi), breaks=c(-pi, -pi*2/3, -pi/3, 0, pi/3, pi*2/3, pi))

scale_x_chi_degrees <- scale_x_continuous(
	"Base -- Acceptor Torsion (Degrees)",
	limit=c(-180,180), breaks=c(-180, -90, 0, 90, 180))

scale_y_chi_degrees <- scale_y_continuous(
	"Base -- Acceptor Torsion (Degrees)",
	limit=c(-180,180), breaks=c(-180, -90, 0, 90, 180))
