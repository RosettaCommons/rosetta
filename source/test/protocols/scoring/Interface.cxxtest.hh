// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/conformation/Interface.cxxtest.hh
/// @brief  test suite for core::conformation::Interface.cc
/// @author Monica Berrondo


// Test headers
#include <cxxtest/TestSuite.h>

#include <test/protocols/init_util.hh>


#include <core/chemical/ChemicalManager.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// Unit headers
#include <protocols/scoring/Interface.hh>

// Package headers

#include <test/UTracer.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


using basic::Error;
using basic::Warning;

// using declarations
using namespace core;
using namespace core::pose;
using namespace core::conformation;
using namespace protocols::scoring;

/// @name InterfaceTest
/// @brief: test the interface calculation between two proteins
/// @details use the docking protein from demo/rosetta/docking
///  Values for the interface have been checked against pymol

class InterfaceTest : public CxxTest::TestSuite
{
public:
	chemical::ResidueTypeSetCAP residue_set;

	PoseOP the_pose;
	InterfaceTest() {}

	void setUp() {
		protocols_init();

		// This accomplishes nothing...
		residue_set = chemical::ChemicalManager::get_instance()->residue_type_set( chemical::CENTROID );

		the_pose = utility::pointer::make_shared< Pose >();
		core::import_pose::centroid_pose_from_pdb( *the_pose, "protocols/scoring/dock_in.pdb" );
	}

	void tearDown() {
		the_pose.reset();
	}

	void test_InterfaceTest() {
		//test::UTracer UT("core/conformation/test_simple_conformation.u");

		// create it with jump 1 as the dock jump
		InterfaceOP iface( utility::pointer::make_shared< Interface >(1) );

		// the Interface object uses an energy graph to calculate the iface,
		// therefore, energies must be accumulated before any calculations
		scoring::ScoreFunctionOP scorefxn ( scoring::ScoreFunctionFactory::create_score_function( "interchain_cen" ) );
		(*scorefxn)(*the_pose);
		/// Now handled automatically.  scorefxn->accumulate_residue_total_energies( *the_pose );

		// use an 8 A cutoff for the iface calculation
		iface->distance( 8.0 );

		iface->calculate( *the_pose );

		test::UTracer ut("protocols/scoring/Interface.u");
		iface->show( ut, *the_pose );

		// This is the expected interface residue position vector
		utility::vector1 < utility::vector1_int > expected_paired_list({
			utility::vector1_int{39, 40, 41, 42, 55, 57, 58, 59, 95, 96, 97, 98, 99, 102, 141, 142, 143, 146, 149, 151, 152, 172, 175, 180, 189, 190, 191, 192, 193, 194, 195, 213, 214, 215, 216, 217, 218, 219, 220, 226, 227},
			utility::vector1_int{255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 275, 277, 278, 280, 281, 284, 288, 299, 300}
			});

		// Verify that the calculated pair_list matches the expected paired list
		TS_ASSERT_EQUALS(expected_paired_list, iface->pair_list());
	}

	void test_OneChainInterfaceTest() {
		// Generate a one chain pose to test
		PoseOP onechain_pose = utility::pointer::make_shared< Pose >();
		core::import_pose::centroid_pose_from_pdb( *onechain_pose, "core/pose/onechain.pdb" );

		// Create a new interface
		InterfaceOP iface( utility::pointer::make_shared<Interface>() );

		// Score the pose so that we can calculate the interface using the energy graph
		scoring::ScoreFunctionOP scorefxn ( scoring::ScoreFunctionFactory::create_score_function( "interchain_cen" ) );
		(*scorefxn)(*onechain_pose);

		// use an 8 A cutoff for the iface calculation
		iface->distance( 8.0 );

		// Calculate the interface (we have one chain, so there should be none)
		iface->calculate( *onechain_pose );

		// Make sure the interface residue lists are empty
		for ( utility::vector1_int res_list: iface->pair_list() ) {
			TS_ASSERT(res_list.empty());
		}
	}
};
