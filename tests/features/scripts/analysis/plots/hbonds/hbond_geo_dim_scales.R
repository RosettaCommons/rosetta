# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


scale_x_AHdist <- scale_x_continuous(
	expression(paste('Acceptor -- Hydrogen Distance (', ring(A), ')')),
	limits=c(1.4, 3), breaks=c(1.4, 1.8, 2.2, 2.6, 3))

scale_y_AHdist <- scale_y_continuous(
	expression(paste('Acceptor -- Hydrogen Distance (', ring(A), ')')),
	limits=c(1.4, 3), breaks=c(1.4, 1.8, 2.2, 2.6, 3))

scale_x_AHdist3 <- scale_x_continuous(
	expression(paste('(Acceptor -- Hydrogen Distance)^3')),
	limits=c(1.4^3, 3.2^3), breaks=c(1.5^3, 2^3, 3^3), labels=c(1.5, 2, 3))

scale_y_AHdist3 <- scale_y_continuous(
	expression(paste('(Acceptor -- Hydrogen Distance)^3')),
	limits=c(1.4^3, 3.2^3), breaks=c(1.4^3, 2^3, 3^3), labels=c(1.4, 2, 3))

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


hbond_scale_x_by_name <- function(dim){
	if(as.character(dim)=="AHdist") scale_x_AHdist
	else if(as.character(dim)=="AHdist") scale_x_AHdist
	else if(as.character(dim)=="cosBAH") scale_x_cosBAH
	else if(as.character(dim)=="cosAHD") scale_x_cosAHD
	else if(as.character(dim)=="chi") scale_x_chi_degrees
}

hbond_scale_y_by_name <- function(dim){
	if(as.character(dim)=="AHdist") scale_y_AHdist
	else if(as.character(dim)=="AHdist") scale_y_AHdist
	else if(as.character(dim)=="cosBAH") scale_y_cosBAH
	else if(as.character(dim)=="cosAHD") scale_y_cosAHD
	else if(as.character(dim)=="chi") scale_y_chi_degrees
}


don_chem_type_name_linear <- function(don_chem_type){
	factor(don_chem_type,
		levels = c("hbdon_IMD", "hbdon_IME", "hbdon_GDE", "hbdon_GDH",
			"hbdon_AHX", "hbdon_HXL", "hbdon_IND", "hbdon_AMO", "hbdon_CXA", "hbdon_PBA"),
		labels = c("dIMD: h", "dIME: h", "dGDE: r", "dGDH: r",
			"dAHX: y", "dHXL: s,t", "dIND: w", "dAMO: k", "dCXA: n,q", "dPBA: bb"))
}

don_chem_type_name_wrap <- function(don_chem_type){
	factor(don_chem_type,
		levels = c("hbdon_HXL", "hbdon_IMD", "hbdon_GDE", "hbdon_AHX",
			"hbdon_IME", "hbdon_GDH", "hbdon_CXA", "hbdon_AMO", "hbdon_IND", "hbdon_PBA"),
		labels = c("dHXL: s,t", "dIMD: h", "dGDE: r", "dAHX: y", "dIME: h", "dGDH: r",
			"dCXA: n,q", "dAMO: k", "dIND: w", "dPBA: bb"))
}


acc_chem_type_name_linear <- function(acc_chem_type){
	factor(acc_chem_type,
		levels = c("hbacc_IMD", "hbacc_IME", "hbacc_AHX", "hbacc_HXL",
			"hbacc_CXA", "hbacc_CXL", "hbacc_PBA"),
		labels = c("aIMD: h", "aIME: h", "aAHX: y", "aHXL: s,t",
			"aCXA: n,q", "aCXL: d,e", "aPBA: bb"))
}

acc_chem_type_name_wrap <- function(acc_chem_type){
	factor(acc_chem_type,
		levels = c("hbacc_HXL", "hbacc_CXL", "hbacc_IMD", "hbacc_AHX",
			"hbacc_CXA", "hbacc_IME",  "hbacc_PBA"),
		labels = c("aHXL: s,t", "aCXL: d,e", "aIMD: h",   "aAHX: y",
			"aCXA: n,q", "aIME: h",    "aPBA: bb"))
}

chem_type_name_linear <- function(chem_type){
	factor(chem_type,
		levels = c("hbdon_IMD", "hbdon_IME", "hbdon_GDE", "hbdon_GDH",
			"hbdon_AHX", "hbdon_HXL", "hbdon_IND", "hbdon_AMO", "hbdon_CXA", "hbdon_PBA",
			"hbacc_IMD", "hbacc_IME", "hbacc_AHX", "hbacc_HXL",
			"hbacc_CXA", "hbacc_CXL", "hbacc_PBA"),
		labels = c("dIMD: h", "dIME: h", "dGDE: r", "dGDH: r",
			"dAHX: y", "dHXL: s,t", "dIND: w", "dAMO: k", "dCXA: n,q", "dPBA: bb",
			"aIMD: h", "aIME: h", "aAHX: y", "aHXL: s,t",
			"aCXA: n,q", "aCXL: d,e", "aPBA: bb"))
}

