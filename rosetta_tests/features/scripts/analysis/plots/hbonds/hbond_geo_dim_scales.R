

scale_x_AHdist <- scale_x_continuous(
	expression(paste('Acceptor -- Hydrogen Distance (', ring(A), ')')),
	limits=c(1.4, 3), breaks=c(1.4, 1.8, 2.2, 2.6, 3))

scale_y_AHdist <- scale_y_continuous(
	expression(paste('Acceptor -- Hydrogen Distance (', ring(A), ')')),
	limits=c(1.4, 3), breaks=c(1.4, 1.8, 2.2, 2.6, 3))

scale_x_ADdist <- scale_x_continuous(
	expression(paste('Acceptor -- Donor Distance (', ring(A), ')')),
	breaks=c(2.3, 2.8, 3.3))

scale_y_ADdist <- scale_y_continuous(
	expression(paste('Acceptor -- Donor Distance (', ring(A), ')')),
	breaks=c(2.3, 2.8, 3.3))

scale_x_cosBAH <- scale_x_continuous(
	"cos(Base -- Acceptor -- Hydrogen)",
	limit=c(-.4,1), breaks=c(-.4 -.2, 0, .2, .4, .6, .8, 1))

scale_y_cosBAH <- scale_y_continuous(
	"cos(Base -- Acceptor -- Hydrogen)",
	limit=c(-.4,1), breaks=c(-.4, -.2, 0, .2, .4, .6, .8, 1))

scale_x_cosAHD <- scale_x_continuous(
	"cos(Acceptor -- Hydrogen -- Donor)",
	limit=c(0,1), breaks=c(0, .2, .4, .6, .8, 1))

scale_y_cosAHD <- scale_y_continuous(
	"cos(Acceptor -- Hydrogen -- Donor)",
	limit=c(0,1), breaks=c(0, .2, .4, .6, .8, 1))

scale_x_chi <- scale_x_continuous(
	"Base -- Acceptor Torsion (Radians)",
	limit=c(0,2*pi), breaks=c(0, pi/3, pi*2/3, pi, pi*4/3, pi*5/3, 2*pi))

scale_y_chi <- scale_y_continuous(
	"Base -- Acceptor Torsion (Radians)",
	limit=c(0,2*pi), breaks=c(0, pi/3, pi*2/3, pi, pi*4/3, pi*5/3, 2*pi))
