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
id = "chi_sinBAH_polar_scatter_bb_by_ss",
author = "Matthew O'Meara",
brief_description = "",
feature_reporter_dependencies = c("HBondFeatures"),
run=function(self, sample_sources, output_dir, output_formats){

d_ply(sample_sources, .variables=("sample_source"), function(sample_source){
	ss <- sample_source[1,"sample_source"]

	sele <-"
	SELECT
		geom.cosBAH,
		geom.chi,
		acc_atoms.base2_x AS p1x, acc_atoms.base2_y AS p1y, acc_atoms.base2_z AS p1z,
		acc_atoms.base_x  AS p2x, acc_atoms.base_y  AS p2y, acc_atoms.base_z  AS p2z,
		acc_atoms.atm_x   AS p3x, acc_atoms.atm_y   AS p3y, acc_atoms.atm_z   AS p3z,
		don_atoms.atm_x   AS p4x, don_atoms.atm_y   AS p4y, don_atoms.atm_z   AS p4z,
		don_ss.dssp AS don_ss, acc_ss.dssp AS acc_ss
	FROM
		hbond_geom_coords AS geom,
		hbonds AS hbond,
		hbond_sites AS don_site,
		hbond_sites AS acc_site,
		hbond_site_atoms AS don_atoms,
		hbond_site_atoms AS acc_atoms,
		residue_secondary_structure AS don_ss,
		residue_secondary_structure AS acc_ss
	WHERE
		hbond.struct_id = geom.struct_id AND hbond.hbond_id = geom.hbond_id AND
		hbond.struct_id = don_site.struct_id AND hbond.don_id = don_site.site_id AND
		hbond.struct_id = acc_site.struct_id AND hbond.acc_id = acc_site.site_id AND
		don_site.HBChemType = 'hbdon_PBA' AND acc_site.HBChemType = 'hbacc_PBA' AND
		don_atoms.struct_id = hbond.struct_id AND don_atoms.site_id = hbond.don_id AND
		acc_atoms.struct_id = hbond.struct_id AND acc_atoms.site_id = hbond.acc_id AND
		don_ss.struct_id = don_site.struct_id AND don_ss.resNum = don_site.resNum AND
		acc_ss.struct_id = acc_site.struct_id AND acc_ss.resNum = acc_site.resNum;"
	f <- query_sample_sources(sample_source, sele)

	normalize <- function(xin,yin,zin) {
		len = sqrt(xin*xin+yin*yin+zin*zin)
		data.frame(x=xin/len, y=yin/len, z=zin/len)
	}

	vlen <- function(xin,yin,zin) {
		sqrt(xin*xin+yin*yin+zin*zin)
	}

	cross <- function(x1,y1,z1,x2,y2,z2){
		data.frame(x=(y1 * z2 ) - (z1 * y2 ),y=(z1 * x2 ) - (x1 * z2 ),z=(x1 * y2 ) - (y1 * x2 ))
	}

	dot_prod <- function(p1,p2){
		p1$x*p2$x + p1$y*p2$y + p1$z*p2$z
	}

	dihedral_angle <- function(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4){
		v1 <- normalize(x2-x1, y2-y1, z2-z1)
		v2 <- normalize(x3-x2, y3-y2, z3-z2)
		v3 <- normalize(x4-x3, y4-y3, z4-z3)
		x <- dot_prod(v1,v2) * dot_prod(v2,v3) - dot_prod(v1,v3)
		y <- dot_prod(v1,cross(v2$x, v2$y, v2$z, v3$x, v3$y, v3$z))
		atan2(y,x)
	}

	alt_chi_dihedral_angle <- function(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4){
		vbx <- (x2+x1)/2; vby <- (y2+y1)/2; vbz <- (z2+z1)/2  # virtual base
		v1x <- x2-x1;     v1y <- y2-y1;     v1z <- z2-z1      # b2 -> b
		v2x <- x3-x2;     v2y <- y3-y2;     v2z <- z3-z2      # b2 -> b
		vb2 <- cross(v1x, v1y, v1z, v2x, v2y, v2z)            # virtual abase2
#		print(paste("vb2:", vb2$x, vb2$y, vb2$z))
#		print(paste("vb: ", vbx, vby, vbz))
#		print(paste("p3: ", x3, y3, z3))
#		print(paste("p4: ", x4, y4, z4))
		dihedral_angle(vb2$x, vb2$y, vb2$z, vbx, vby, vbz, x3, y3, z4, x4, y4, z4)
	}

#	print(dihedral_angle(1,0,0, 0,1,0, 0,0,0, -1,-1,-1))
#	print(alt_chi_dihedral_angle(1,0,0, 0,1,0, 0,0,0, -1,-1,1))
#	print(alt_chi_dihedral_angle(c(1,1),c(0,0),c(0,0), c(0,0),c(1,1),c(0,0), c(0,0),c(0,0),c(0,0), c(-1,-1),c(-1,-1), c(1,-1)))

	f[f$hybrid == "sp3" | f$hybrid == "ring", "chi"] <- with(f[f$hybrid == "sp3" | f$hybrid == "ring",],
		alt_chi_dihedral_angle(p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z, p4x, p4y, p4z))

	alt_sp3_cosBAH <- function(x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4){
		vbx <- (x2+x1)/2; vby <- (y2+y1)/2; vbz <- (z2+z1)/2  # virtual base
		v1 <- normalize(x3-vbx, y3-vby, z3-vbz)
		v2 <- normalize(x4-x3, y4-y3, z4-z3)
		dot_prod(v1, v2)
	}
	f[f$hybrid == "sp3", "cosBAH"] <- with(f[f$hybrid == "sp3",],
		alt_sp3_cosBAH(p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z, p4x, p4y, p4z))

	#equal area projection
	f <- transform(f,
		capx = 2*sin(acos(cosBAH)/2)*cos(chi),
		capy = 2*sin(acos(cosBAH)/2)*sin(chi))

	sub_f <- ddply(f, .variables=c("don_ss", "don_ss"),
		function(df){sample_rows(df, 10000)})

	plot_id = "hbond_sinBAH_eq_polar_scatter_bb_by_ss"
	ggplot(data=sub_f) + theme_bw() +
		polar_equal_area_grids_bw() +
		geom_point(aes(x=capx, y=capy), size=.4, alpha=.5) +
		facet_grid(don_ss ~ acc_ss) +
		ggtitle(paste("Backbone-Backbone Hydrogen Bonds chi vs sinBAH Angles by Secondary Structure\nEqual Coordinate Projection   Sample Source: ", ss, sep="")) +
		scale_x_continuous('2*sin(BAH/2) * cos(CHI)', breaks=c(-1, 0, 1)) +
		scale_y_continuous('2*sin(BAH/2) * sin(CHI)', breaks=c(-1, 0, 1))
	save_plots(self, plot_id, sample_source, output_dir, output_formats)

	#orthographic projection
	f <- transform(f,
		capx = sin(acos(cosBAH))*cos(chi),
		capy = sin(acos(cosBAH))*sin(chi))

	sub_f <- ddply(f, .variables=c("don_ss", "don_ss"),
		function(df){sample_rows(df, 10000)})

	plot_id = "hbond_sinBAH_ortho_polar_scatter_bb_by_ss"
	ggplot(data=sub_f) + theme_bw() +
		geom_point(aes(x=capx, y=capy), size=.5, alpha=.4) +
		facet_grid(don_ss ~ acc_ss) +
		ggtitle(paste("Backbone-Backbone Hydrogen Bonds chi vs sinBAH Angles by Secondary Structure\nOrthographic Projection   Sample Source: ", ss, sep="")) +
		scale_x_continuous('sin(BAH) * cos(CHI)', breaks=c(-1, 0, 1)) +
		scale_y_continuous('sin(BAH) * sin(CHI)', breaks=c(-1, 0, 1))
	save_plots(self, plot_id, sample_source, output_dir, output_formats)
})


})) # end FeaturesAnalysis
