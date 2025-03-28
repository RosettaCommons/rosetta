// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/scoring/InterchainPairEnergy.cxxtest.hh
/// @brief  test suite for protocols/scoring/InterchainPairEnergy.cc
/// @author Monica Berrondo


// Test headers
#include <cxxtest/TestSuite.h>

#include <test/protocols/init_util.hh>


#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>

// Unit headers

// Package headers


//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <protocols/scoring/Interface.fwd.hh>

#include <basic/Tracer.hh> // AUTO IWYU For Error, Warning


using basic::Error;
using basic::Warning;

// using declarations
using namespace core;
using namespace core::pose;
using namespace core::conformation;
using namespace protocols::scoring;

/// @name InterchainPairEnergy
/// @brief: test the Interchain Pair Energy between two proteins

class InterchainPairEnergy : public CxxTest::TestSuite
{
public:
	PoseOP the_pose;
	InterchainPairEnergy() {}

	void setUp() {
		protocols_init();

		the_pose = utility::pointer::make_shared< Pose >();
		core::import_pose::centroid_pose_from_pdb( *the_pose, "protocols/scoring/dock_in.pdb" );
	}

	void tearDown() {
		the_pose.reset();
	}

	void test_InterchainPairEnergyTest() {
		using namespace core::pose;
		using namespace core::scoring;
		ScoreFunction sfxn;
		sfxn.set_weight( interchain_pair, 1.0 );
		TS_ASSERT_DELTA( sfxn( *the_pose ), 2.3075, 1e-3 );
	}

	void test_OneChainPairEnergyTest() {
		PoseOP onechain_pose = utility::pointer::make_shared< Pose >();
		core::import_pose::centroid_pose_from_pdb( *onechain_pose, "core/pose/onechain.pdb" );

		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( core::scoring::interchain_pair, 1.0 );
		TS_ASSERT_EQUALS( sfxn( *onechain_pose ), 0 );
	}
};

