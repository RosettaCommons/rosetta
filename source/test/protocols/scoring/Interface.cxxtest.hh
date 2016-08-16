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
#include <core/types.hh>

// Unit headers
#include <protocols/scoring/Interface.hh>

// Package headers

#include <test/UTracer.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <utility/vector1.hh>


using basic::T;
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

		the_pose = PoseOP( new Pose );
		core::import_pose::centroid_pose_from_pdb( *the_pose, "protocols/scoring/dock_in.pdb" );
	}

	void tearDown() {
		the_pose.reset();
	}

	void test_InterfaceTest() {
		//test::UTracer UT("core/conformation/test_simple_conformation.u");

		// create it with jump 1 as the dock jump
		InterfaceOP iface( new Interface(1) );

		// the Interface object uses an energy graph to calculate the iface,
		// therefore, energies must be accumulated before any calculations
		scoring::ScoreFunctionOP scorefxn ( scoring::ScoreFunctionFactory::create_score_function( "interchain_cen" ) );
		(*scorefxn)(*the_pose);
		/// Now handled automatically.  scorefxn->accumulate_residue_total_energies( *the_pose );

		// use an 8 A cutoff for the iface calculation
		iface->distance( 8.0 );

		// monitor the output from the iface calculation
		basic::otstreamOP ut( new test::UTracer("protocols/scoring/Interface.u") );
		basic::Tracer::set_ios_hook(ut, "core.conformation.Interface");
		iface->calculate( *the_pose );
		iface->print( *the_pose );
		basic::Tracer::set_ios_hook(0, "");
	}
};

