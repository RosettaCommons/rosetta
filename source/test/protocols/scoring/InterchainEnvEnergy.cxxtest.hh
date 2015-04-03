// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/protocols/scoring/InterchainEnvEnergy.cxxtest.hh
/// @brief  test suite for protocols/scoring/InterchainEnvEnergy.cc
/// @author Monica Berrondo


// Test headers
#include <cxxtest/TestSuite.h>

#include <test/protocols/init_util.hh>


#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/types.hh>

// Unit headers

// Package headers

#include <test/UTracer.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <protocols/scoring/Interface.fwd.hh>
#include <utility/vector1.hh>


using basic::T;
using basic::Error;
using basic::Warning;

// using declarations
using namespace core;
using namespace core::pose;
using namespace core::conformation;
using namespace protocols::scoring;

/// @name InterchainEnvEnergy
/// @brief: test the Interchain Environment Energy between two proteins

class InterchainEnvEnergy : public CxxTest::TestSuite
{
public:
	PoseOP the_pose;
	InterchainEnvEnergy() {}

	void setUp() {
		protocols_init();

		the_pose = PoseOP( new Pose );
		core::import_pose::centroid_pose_from_pdb( *the_pose, "protocols/scoring/dock_in.pdb" );
	}

	void tearDown() {
		the_pose.reset();
	}

	void test_InterchainEnvEnergyTest() {
		using namespace core::pose;
		using namespace core::scoring;
		ScoreFunction sfxn;
		sfxn.set_weight( interchain_env, 1.0 );
		TS_ASSERT_DELTA( sfxn( *the_pose ), 6.2852, 1e-3 );
	}
};

