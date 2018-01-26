// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/conformation/symmetry/SymmDataLoader.cxxtest.hh
/// @brief test suite for core::conformation::symmetry::SymmDataLoader
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Project headers
#include <core/conformation/symmetry/SymmData.hh>
#include <core/conformation/symmetry/SymmDataLoader.hh>


// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/pointer/owning_ptr.hh>

// Basic headers
#include <basic/resource_manager/ResourceManager.hh>

// C++ headers
#include <sstream>

using namespace basic::resource_manager;

class SymmDataLoaderTests : public CxxTest::TestSuite {

public:

	void setUp() {
		protocols_init();
	}

	std::string
	example_symmdef() {
		return
			"symmetry_name stest1__3\n"
			"E = 3*VRT0_base + 3*(VRT0_base:VRT1_base)\n"
			"anchor_residue 18\n"
			"virtual_coordinates_start\n"
			"xyz VRT0  0.9872615,0.0594791,-0.1475702  -0.0601376,0.9981901,0.0000000  0.6047347,-0.2520896,1.8401966\n"
			"xyz VRT0_base  0.9872615,0.0594791,-0.1475702  -0.0601376,0.9981901,0.0000000  -5.1031250,-0.5959687,2.6933750\n"
			"xyz VRT1  -0.4383883,-0.8937378,0.0951230  0.8902436,-0.4463399,-0.0908131  0.6047347,-0.2520896,1.8401966\n"
			"xyz VRT1_base  -0.4383883,-0.8937378,0.0951230  0.8902436,-0.4463399,-0.0908131  3.1392801,4.9150622,1.2902424\n"
			"xyz VRT2  -0.5488732,0.8342587,0.0524473  -0.8192128,-0.5493288,0.1647067  0.6047347,-0.2520896,1.8401966\n"
			"xyz VRT2_base  -0.5488732,0.8342587,0.0524473  -0.8192128,-0.5493288,0.1647067  3.7780490,-5.0753624,1.5369723\n"
			"xyz VRT  1.0000000,0.0000000,0.0000000  0.0000000,1.0000000,0.0000000  1.6047347,-0.2520896,1.8401966\n"
			"virtual_coordinates_stop\n"
			"connect_virtual JUMP0_to_com VRT0 VRT0_base\n"
			"connect_virtual JUMP0_to_subunit VRT0_base SUBUNIT\n"
			"connect_virtual JUMP1_to_com VRT1 VRT1_base\n"
			"connect_virtual JUMP1_to_subunit VRT1_base SUBUNIT\n"
			"connect_virtual JUMP2_to_com VRT2 VRT2_base\n"
			"connect_virtual JUMP2_to_subunit VRT2_base SUBUNIT\n"
			"connect_virtual JUMP0 VRT VRT0\n"
			"connect_virtual JUMP1 VRT0 VRT1\n"
			"connect_virtual JUMP2 VRT0 VRT2\n"
			"set_dof JUMP0_to_com x(5.78150747337524)\n"
			"set_dof JUMP0_to_subunit angle_x angle_y angle_z\n"
			"set_jump_group JUMPGROUP2 JUMP0_to_com JUMP1_to_com JUMP2_to_com\n"
			"set_jump_group JUMPGROUP3 JUMP1_to_subunit JUMP0_to_subunit JUMP2_to_subunit\n";
	}

	void test_SymmDataLoader(){
		using namespace core::conformation::symmetry;
		std::istringstream example_symmdef_stream( example_symmdef() );

		basic::resource_manager::ResourceManager rm;
		SymmDataLoader loader;
		utility::tag::TagCOP tag( utility::tag::Tag::create( "<SymmData name=\"bogus\"/>" ) );

		basic::resource_manager::ResourceCOP my_resource = loader.create_resource( rm, tag, "unit test", example_symmdef_stream );
		core::conformation::symmetry::SymmDataCOP symm_data(
			utility::pointer::dynamic_pointer_cast< core::conformation::symmetry::SymmData const > ( my_resource ));

		// make sure we got back the right resource type
		TS_ASSERT( symm_data );

		// Setup the alternate symmetry data to compare against
		core::conformation::symmetry::SymmDataOP symm_data_alt( new core::conformation::symmetry::SymmData() );
		std::istringstream iss2( example_symmdef() );
		symm_data_alt->read_symmetry_data_from_stream( iss2 );

		TS_ASSERT_EQUALS(*symm_data, *symm_data_alt);
	}

};
