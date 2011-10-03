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
	geom.cosAHD,
	don_site.HBChemType as don_chem_type,
	acc_site.HBChemType as acc_chem_type,
 	acc_atoms.base_x AS abx, acc_atoms.base_y AS aby, acc_atoms.base_z AS abz, -- acceptor base atom
	acc_atoms.atm_x  AS ax,  acc_atoms.atm_y  AS ay,  acc_atoms.atm_z  AS az,  -- acceptor atom
	don_atoms.atm_x  AS hx,  don_atoms.atm_y  AS hy,  don_atoms.atm_z  AS hz,  -- hydrogen atom
	don_atoms.base_x AS dx,  don_atoms.base_y AS dy,  don_atoms.base_z AS dz   -- donor atom
FROM
	hbond_geom_coords AS geom,
	hbonds AS hbond,
	hbond_sites AS don_site,
	hbond_sites AS acc_site,
	hbond_site_atoms AS don_atoms, hbond_site_atoms AS acc_atoms
WHERE
	don_site.HBChemType != 'hbdon_PBA' AND acc_site.HBChemType != 'hbacc_PBA' AND
	hbond.struct_id = geom.struct_id AND hbond.hbond_id = geom.hbond_id AND
	hbond.struct_id = don_site.struct_id AND hbond.don_id = don_site.site_id AND
	hbond.struct_id = acc_site.struct_id AND hbond.acc_id = acc_site.site_id AND
	don_atoms.struct_id = hbond.struct_id AND don_atoms.site_id = hbond.don_id AND
	acc_atoms.struct_id = hbond.struct_id AND acc_atoms.site_id = hbond.acc_id;";
f <- query_sample_sources(sample_sources, sele)

f$don_chem_type <- factor(f$don_chem_type,
	levels = c("hbdon_HXL", "hbdon_IMD", "hbdon_GDE", "hbdon_AHX",
		"hbdon_IME", "hbdon_GDH", "hbdon_CXA", "hbdon_AMO", "hbdon_IND", "hbdon_PBA"),
	labels = c("dHXL: s,t", "dIMD: h", "dGDE: r", "dAHX: y", "dIME: h", "dGDH: r",
		  "dCXA: n,q", "dAMO: k", "dIND: w", "dPBA: bb"))

f$acc_chem_type <- factor(f$acc_chem_type,
	levels = c("hbacc_HXL", "hbacc_CXL", "hbacc_IMD", "hbacc_AHX", "hbacc_CXA", "hbacc_IME",  "hbacc_PBA"),
	labels = c("aHXL: s,t", "aCXL: d,e", "aIMD: h",   "aAHX: y",   "aCXA: n,q", "aIME: h",    "aPBA: bb"))

f <- transform(f,
  AHchi = vector_dihedral(
    cbind(abx, aby, abz), cbind(ax, ay, az),
    cbind(hx, hy, hz), cbind(dx, dy, dz)))
               
#equal area projection
f <- transform(f,
	capx = 2*sin(acos(cosAHD)/2)*cos(AHchi),
	capy = 2*sin(acos(cosAHD)/2)*sin(AHchi))

capx_limits <- c(-1.5,1.5)
capy_limits <- capx_limits

plot_parts <- list(
	theme_bw(),
	stat_density2d(aes(x=capx,y=capy, fill=log(..density..+1)), geom="tile", contour=FALSE),
	polar_equal_area_grids_bw,
	geom_indicator(aes(indicator=counts)),
	scale_x_continuous('2*sin(AHD/2) * cos(AHchi)', limits=capx_limits, breaks=c(-1, 0, 1)),
	scale_y_continuous('2*sin(AHD/2) * sin(AHchi)', limits=capy_limits, breaks=c(-1, 0, 1)),
	coord_fixed(ratio = 1),
	scale_fill_gradientn('log(Density+1)', colour=jet.colors(10)))

narrow_output_formats <- transform(output_formats, width=height)

d_ply(f, .(sample_source), function(sub_f){
	ss_id <- sub_f$sample_source[1]
        
	plot_id = paste("AHchi_AHD_eq_polar_density_by_don", ss_id, sep="_")
        sub_f <- ddply(sub_f, .(don_chem_type), transform, counts = length(sample_source))
	ggplot(data=sub_f) + plot_parts +
		opts(title = paste("Hydrogen Bonds AHchi vs AHD Angles Sidechain-Sidechain\nEqual Coordinate Projection   Sample Source: ", ss_id, sep="")) +
                facet_wrap( ~ don_chem_type)
	save_plots(plot_id, sample_sources[sample_sources$sample_source == ss_id,],
		output_dir, narrow_output_formats)

        plot_id = paste("AHchi_AHD_eq_polar_density_by_acc", ss_id, sep="_")
        sub_f <- ddply(sub_f, .(acc_chem_type), transform, counts = length(sample_source))
        ggplot(data=sub_f) + plot_parts +
		opts(title = paste("Hydrogen Bonds AHchi vs AHD Angles Sidechain-Sidechain\nEqual Coordinate Projection   Sample Source: ", ss_id, sep="")) +
                facet_wrap( ~ acc_chem_type)
	save_plots(plot_id, sample_sources[sample_sources$sample_source == ss_id,],
		output_dir, output_formats)

        plot_id = paste("AHchi_AHD_eq_polar_density", ss_id, sep="_")
        sub_f$counts <- nrow(sub_f)
	ggplot(data=sub_f) + plot_parts +
		opts(title = paste("Hydrogen Bonds AHchi vs AHD Angles Sidechain-Sidechain\nEqual Coordinate Projection   Sample Source: ", ss_id, sep=""))
	save_plots(plot_id, sample_sources[sample_sources$sample_source == ss_id,],
		output_dir, output_formats)

        sub_f <- ddply(sub_f, .(don_chem_type, acc_chem_type), transform, counts = length(sample_source))

        plot_id = paste("AHchi_AHD_eq_polar_density_dGDE", ss_id, sep="_")
	ggplot(data=subset(sub_f, don_chem_type == "dGDE: r")) + plot_parts +
		opts(title = paste("Hydrogen Bonds AHchi vs AHD Angles GDE Donor\nEqual Coordinate Projection   Sample Source: ", ss_id, sep="")) +
               	facet_wrap(~acc_chem_type)
	save_plots(plot_id, sample_sources[sample_sources$sample_source == ss_id,],
		output_dir, output_formats)

        plot_id = paste("AHchi_AHD_eq_polar_density_dGDH", ss_id, sep="_")
	ggplot(data=subset(sub_f, don_chem_type == "dGDH: r")) + plot_parts +
		opts(title = paste("Hydrogen Bonds AHchi vs AHD Angles GDH Donor\nEqual Coordinate Projection   Sample Source: ", ss_id, sep="")) +
               	facet_wrap(~acc_chem_type)
	save_plots(plot_id, sample_sources[sample_sources$sample_source == ss_id,],
		output_dir, output_formats)



        
})
