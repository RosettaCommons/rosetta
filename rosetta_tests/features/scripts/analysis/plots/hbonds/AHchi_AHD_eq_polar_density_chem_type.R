# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

check_setup()

sele <-"
SELECT
	geom.cosAHD, geom.chi,
	don_site.HBChemType AS don_chem_type,
	acc_site.HBChemType AS acc_chem_type,
 	acc_atoms.base_x AS abx, acc_atoms.base_y AS aby, acc_atoms.base_z AS abz, -- acceptor base atom
	acc_atoms.atm_x  AS ax,  acc_atoms.atm_y  AS ay,  acc_atoms.atm_z  AS az,  -- acceptor atom
	don_atoms.atm_x  AS hx,  don_atoms.atm_y  AS hy,  don_atoms.atm_z  AS hz,  -- hydrogen atom
	don_atoms.base_x AS dx,  don_atoms.base_y AS dy,  don_atoms.base_z AS dz,  -- donor atom
	hbond.accRank AS acc_rank, don_site.resNum - acc_site.resNum AS seq_sep,
	struct.tag,
	don_pdb.chain AS don_chain,
	don_pdb.resNum AS don_resNum,
	don_pdb.iCode AS don_iCode,
	acc_pdb.chain AS acc_chain,
	acc_pdb.resNum AS acc_resNum,
	acc_pdb.iCode AS acc_iCode
FROM
	structures AS struct,
	hbond_geom_coords AS geom,
	hbonds AS hbond,
	hbond_sites AS don_site,
	hbond_sites AS acc_site,
	hbond_site_atoms AS don_atoms, hbond_site_atoms AS acc_atoms,
	hbond_sites_pdb AS don_pdb, hbond_sites_pdb AS acc_pdb
WHERE
	struct.struct_id = hbond.struct_id AND
	don_site.HBChemType != 'hbdon_PBA' AND acc_site.HBChemType != 'hbacc_PBA' AND
	hbond.struct_id = geom.struct_id AND hbond.hbond_id = geom.hbond_id AND
	hbond.struct_id = don_site.struct_id AND hbond.don_id = don_site.site_id AND
	hbond.struct_id = acc_site.struct_id AND hbond.acc_id = acc_site.site_id AND
	don_atoms.struct_id = hbond.struct_id AND don_atoms.site_id = hbond.don_id AND
	acc_atoms.struct_id = hbond.struct_id AND acc_atoms.site_id = hbond.acc_id AND
  don_pdb.struct_id = hbond.struct_id AND don_pdb.site_id = don_site.site_id AND
  acc_pdb.struct_id = hbond.struct_id AND acc_pdb.site_id = acc_site.site_id AND
	don_pdb.heavy_atom_temperature < 20 AND acc_pdb.heavy_atom_temperature < 20;"
f <- query_sample_sources(sample_sources, sele)

f$don_chem_type <- factor(f$don_chem_type,
	levels = c("hbdon_HXL", "hbdon_IMD", "hbdon_GDE", "hbdon_AHX",
		"hbdon_IME", "hbdon_GDH", "hbdon_CXA", "hbdon_AMO", "hbdon_IND", "hbdon_PBA"),
	labels = c("dHXL: s,t", "dIMD: h", "dGDE: r", "dAHX: y", "dIME: h", "dGDH: r",
		"dCXA: n,q", "dAMO: k", "dIND: w", "dPBA: bb"))

f$acc_chem_type <- factor(f$acc_chem_type,
	levels = c("hbacc_HXL", "hbacc_CXL", "hbacc_IMD", "hbacc_AHX", "hbacc_CXA", "hbacc_IME",  "hbacc_PBA"),
	labels = c("aHXL: s,t", "aCXL: d,e", "aIMD: h",   "aAHX: y",   "aCXA: n,q", "aIME: h",    "aPBA: bb"))

f$acc_rank <- factor(f$acc_rank)

f$orbital <- factor(ifelse(f$chi < pi/2 & f$chi > -pi/2, "anti", "syn"))
f[ abs(f$seq_sep) > 5, "seq_sep"] = factor("long_range")

f <- transform(f,
	AHchi = vector_dihedral(
		cbind(abx, aby, abz), cbind(ax, ay, az),
		cbind(hx, hy, hz), cbind(dx, dy, dz)))

#write.table(
#	f[acos(f$cosAHD)*180/pi > 20 & f$AHchi > pi/2 & f$AHchi < pi*3/4,
#		c("tag","don_chain","don_resNum", "don_iCode",
#			"acc_chain", "acc_resNum", "acc_iCode")],
#	file="/tmp/instance.csv", quote=FALSE, sep=",")

write.table( with(f,
	f[orbital=="anti" & seq_sep==4 & don_chem_type == "dGDH: r" &
		acc_chem_type == "aCXL: d,e",
		c("tag","don_chain","don_resNum", "don_iCode",
			"acc_chain", "acc_resNum", "acc_iCode") ]),
	file="/tmp/instances.csv", quote=FALSE, sep=",")

#equal area projection
f <- transform(f,
	capx = 2*sin(acos(cosAHD)/2)*cos(AHchi),
	capy = 2*sin(acos(cosAHD)/2)*sin(AHchi))

capx_limits <- c(-1.5,1.5)
capy_limits <- capx_limits

plot_parts <- list(
	theme_bw(),
	stat_density2d(
		aes(x=capx,y=capy, fill=log(..density..+1)), geom="tile", contour=FALSE),
	polar_equal_area_grids_bw(bgcolor="#00007F"),
	geom_indicator(aes(indicator=counts), color="white"),
	scale_x_continuous('2*sin(AHD/2) * cos(AHchi)', limits=capx_limits, breaks=c(-1, 0, 1)),
	scale_y_continuous('2*sin(AHD/2) * sin(AHchi)', limits=capy_limits, breaks=c(-1, 0, 1)),
	coord_fixed(ratio = 1),
	scale_fill_gradientn('log(Density+1)', colour=jet.colors(10)))

narrow_output_formats <- transform(output_formats, width=height)

d_ply(f, .(sample_source), function(sub_f){
	ss_id <- sub_f$sample_source[1]
	ss <- sample_sources[sample_sources$sample_source == ss_id,]

	##################
	plot_id = paste("AHchi_AHD_eq_polar_density_by_don", ss_id, sep="_")
	sub_f <- ddply(sub_f, .(don_chem_type), transform, counts = length(sample_source))
	ggplot(data=sub_f) + plot_parts +
		opts(title = paste("Hydrogen Bonds AHchi vs AHD Angles Sidechain-Sidechain\nEqual Coordinate Projection   Sample Source: ", ss_id, sep="")) +
		facet_wrap( ~ don_chem_type)
	save_plots(plot_id, ss, output_dir, narrow_output_formats)

	#################
	plot_id = paste("AHchi_AHD_eq_polar_density_by_acc", ss_id, sep="_")
	sub_f <- ddply(sub_f, .(acc_chem_type), transform, counts = length(sample_source))
	ggplot(data=sub_f) + plot_parts +
		opts(title = paste("Hydrogen Bonds AHchi vs AHD Angles Sidechain-Sidechain\nEqual Coordinate Projection   Sample Source: ", ss_id, sep="")) +
		facet_wrap( ~ acc_chem_type)
	save_plots(plot_id, ss,	output_dir, output_formats)

	#################
	plot_id = paste("AHchi_AHD_eq_polar_density", ss_id, sep="_")
	sub_f$counts <- nrow(sub_f)
	ggplot(data=sub_f) + plot_parts +
		opts(title = paste("Hydrogen Bonds AHchi vs AHD Angles Sidechain-Sidechain\nEqual Coordinate Projection   Sample Source: ", ss_id, sep=""))
	save_plots(plot_id, ss, output_dir, output_formats)

	#################
	sub_f <- ddply(sub_f, .(don_chem_type, acc_chem_type),
		transform, counts = length(sample_source))

	plot_id = paste("AHchi_AHD_eq_polar_density_dGDE", ss_id, sep="_")
	ggplot(data=subset(sub_f, don_chem_type == "dGDE: r")) + plot_parts +
		opts(title = paste("Hydrogen Bonds AHchi vs AHD Angles GDE Donor\nEqual Coordinate Projection   Sample Source: ", ss_id, sep=""))
		facet_wrap( ~ acc_chem_type)
	save_plots(plot_id, ss, output_dir, output_formats)

	plot_id = paste("AHchi_AHD_eq_polar_density_dGDH", ss_id, sep="_")
	ggplot(data=subset(sub_f, don_chem_type == "dGDH: r")) + plot_parts +
		opts(title = paste("Hydrogen Bonds AHchi vs AHD Angles GDH Donor\nEqual Coordinate Projection   Sample Source: ", ss_id, sep="")) +
		facet_wrap( ~ acc_chem_type)
	save_plots(plot_id, ss, output_dir, output_formats)

	################
	sub_f <- ddply(sub_f, .(don_chem_type, acc_chem_type, acc_rank),
		transform, counts = length(sample_source))

	plot_id = paste("AHchi_AHD_eq_polar_density_dGDE_aCXL_by_acc_rank", ss_id, sep="_")
	ggplot(data=subset(sub_f, don_chem_type == "dGDE: r" & acc_chem_type == "aCXL: d,e")) + plot_parts +
		opts(title = paste("Hydrogen Bonds AHchi vs AHD Angles aCXL-dGDE by Acc Rank\nEqual Coordinate Projection   Sample Source: ", ss_id, sep="")) +
		facet_wrap( ~ acc_rank)
	save_plots(plot_id, ss, output_dir, output_formats)

	plot_id = paste("AHchi_AHD_eq_polar_density_dGDH_aCXL_by_acc_rank", ss_id, sep="_")
	ggplot(data=subset(sub_f, don_chem_type == "dGDH: r" & acc_chem_type == "aCXL: d,e")) + plot_parts +
		opts(title = paste("Hydrogen Bonds AHchi vs AHD Angles aCXL-dGDH by Acc Rank\nEqual Coordinate Projection   Sample Source: ", ss_id, sep="")) +
		facet_wrap( ~ acc_rank)
	save_plots(plot_id, ss, output_dir, output_formats)

	################
	sub_f <- ddply(sub_f, .(don_chem_type, orbital),
		transform, counts = length(sample_source))

	plot_id = paste("AHchi_AHD_eq_polar_density_dGDE_aCXL_by_acc_orbital", ss_id, sep="_")
	ggplot(data=subset(sub_f, don_chem_type == "dGDE: r" & acc_chem_type == "aCXL: d,e")) + plot_parts +
		opts(title = paste("Hydrogen Bonds AHchi vs AHD Angles aCXL-dGDE by Acc Orbital\nEqual Coordinate Projection   Sample Source: ", ss_id, sep="")) +
		facet_wrap( ~ orbital)
	save_plots(plot_id, ss, output_dir, output_formats)

	plot_id = paste("AHchi_AHD_eq_polar_density_dGDH_aCXL_by_acc_orbital", ss_id, sep="_")
	ggplot(data=subset(sub_f, don_chem_type == "dGDH: r" & acc_chem_type == "aCXL: d,e")) + plot_parts +
		opts(title = paste("Hydrogen Bonds AHchi vs AHD Angles aCXL-dGDH by Acc Orbital\nEqual Coordinate Projection   Sample Source: ", ss_id, sep="")) +
		facet_wrap( ~ orbital)
	save_plots(plot_id, ss, output_dir, output_formats)

	################
	sub_sub_f <- subset(sub_f, orbital == "anti" & don_chem_type == "dGDH: r" & acc_chem_type == "aCXL: d,e")
	sub_sub_f <- ddply(sub_sub_f, .(seq_sep), transform, counts = length(sample_source))

	plot_id = paste("AHchi_AHD_eq_polar_density_dGDH_aCXL_acc_anti_orbital_by_seq_sep", ss_id, sep="_")
	ggplot(data=sub_sub_f) + plot_parts +
		opts(title = paste("Hydrogen Bonds AHchi vs AHD Angles aCXL-dGDH Acc Anti Orbital by don.resNum - acc.resNum\nEq Coord Proj   Sample Source: ", ss_id, sep="")) +
		facet_wrap( ~ seq_sep)
	save_plots(plot_id, ss, output_dir, output_formats)
})
