// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

#include <protocols/ligand_docking/PlaceWaters.hh>
#include <protocols/moves/MoverFactory.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <core/pose/Pose.hh>
#include <core/grid/CartGrid.hh>

//c++ headers
#include <string>
#include <iostream>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <utility/vector1.hh>

class PlaceWaters : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void test_tag() {
		protocols::moves::MoverOP base_mover_op = protocols::moves::MoverFactory::get_instance()->newMover("PlaceWaters");
		protocols::ligand_docking::PlaceWatersOP pw_op = protocols::ligand_docking::PlaceWatersOP(utility::pointer::dynamic_pointer_cast<protocols::ligand_docking::PlaceWaters> (base_mover_op) );

		//make tag without input size
		std::string xml_file = "<PlaceWaters chain=\"X\" water=\"TP3.params\"/>";
		std::stringstream ss( xml_file );
		utility::tag::TagOP tag ( new utility::tag::Tag() );
		tag->read (ss);

		basic::datacache::DataMap data;

		pw_op->parse_my_tag( tag, data);
		TS_ASSERT_EQUALS ( pw_op->get_grid_size(), 0.0 );

		xml_file = "<PlaceWaters grid_size=\"10.0\"  chain=\"X\" water=\"TP3.params\" />";
		std::stringstream ss2( xml_file );
		tag->read (ss2);

		pw_op->parse_my_tag( tag, data);
		TS_ASSERT_EQUALS ( pw_op->get_grid_size(), 10.0 );


	}

	//clustering and merging
	void test_clustering() {
		protocols::moves::MoverOP base_mover_op = protocols::moves::MoverFactory::get_instance()->newMover("PlaceWaters");
		protocols::ligand_docking::PlaceWatersOP pw_op = protocols::ligand_docking::PlaceWatersOP(utility::pointer::dynamic_pointer_cast<protocols::ligand_docking::PlaceWaters> (base_mover_op) );

		utility::vector1<core::Vector> spaces;

		for ( core::Size i = 0; i < 9; i+=3 ) {
			spaces.push_back(core::Vector(i, i, i));
		}
		spaces.push_back(core::Vector(9.25, 9.25, 9.25));
		spaces.push_back(core::Vector(8.75, 8.75, 8.75));

		pw_op->cluster_groups(spaces, 2.0, 2.0);
		TS_ASSERT_EQUALS( spaces.size(), 4 );
		TS_ASSERT(spaces.has_value(core::Vector(9.0, 9.0, 9.0)));

		spaces.push_back(core::Vector(1, 1, 1));
		pw_op->merge_close(spaces, 2.0);

		TS_ASSERT_EQUALS( spaces.size(), 4 );
		TS_ASSERT(spaces.contains(core::Vector(0.5, 0.5, 0.5)));
	}

	//hbond acceptors and donors
	void test_hbond_pruning(){
		protocols::moves::MoverOP base_mover_op = protocols::moves::MoverFactory::get_instance()->newMover("PlaceWaters");
		protocols::ligand_docking::PlaceWatersOP pw_op = protocols::ligand_docking::PlaceWatersOP(utility::pointer::dynamic_pointer_cast<protocols::ligand_docking::PlaceWaters> (base_mover_op));

		utility::vector1<core::Vector> spaces;

		for ( core::Size i = 0; i < 4; i+=3 ) {
			spaces.push_back(core::Vector(i, i, i));
		}

		spaces.push_back(core::Vector(-4, -4, -4));

		core::grid::CartGrid<core::Real> grid;
		initialize_grid_copy(grid, core::Vector(0, 0, 0), 10.0, 0.15);

		core::Real sphere_radius = 0.15*0.5;
		pw_op->set_sphere(grid, core::Vector(0,3,1.5), sphere_radius,-1.0);
		pw_op->set_sphere(grid, core::Vector(3,0,1.5), sphere_radius, -4.0);

		utility::vector1<protocols::ligand_docking::water_info> bridging;
		utility::vector1<core::Vector> kept;
		for ( core::Vector v : spaces ) {
			protocols::ligand_docking::water_info wi;
			wi.centroid = v;
			pw_op->find_nearby_hbonds(grid,v,wi,1.5,3.5);
			if ( prot_lig_hbonds_copy(wi) ) {
				bridging.push_back(wi);
				kept.push_back(v);
			}
		}
		TS_ASSERT_EQUALS(bridging.size(), 2);
		TS_ASSERT(!kept.has_value(core::Vector(-4, -4, -4)));

	}

	void initialize_grid_copy(core::grid::CartGrid<core::Real> & grid, core::Vector const & center,
		core::Real width,core::Real resolution ){

		auto num_pts = static_cast<core::Size>(width/resolution);

		core::Real const grid_halfwidth = width / 2.0;


		grid.setBase(
			center.x() - grid_halfwidth,
			center.y() - grid_halfwidth,
			center.z() - grid_halfwidth
		);
		grid.setDimensions(num_pts,num_pts,num_pts,resolution,resolution,resolution);
		grid.setupZones();
		grid.zero();
	}

	bool prot_lig_hbonds_copy(protocols::ligand_docking::water_info & wi){

		bool ligand_contacts = (wi.ligand_donors.size() || wi.ligand_acceptors.size());
		bool protein_contacts = (wi.protein_donors.size() || wi.protein_acceptors.size());

		return (ligand_contacts && protein_contacts);
	}

};

