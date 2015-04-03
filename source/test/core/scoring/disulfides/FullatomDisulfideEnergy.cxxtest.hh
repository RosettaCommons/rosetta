// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/core/scoring/methods/PairEnergy.cxxtest.hh
/// @brief  test suite for core::scoring::PairEnergy.cc
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers

#include <platform/types.hh>

// Package Headers
#include <test/util/pdb1rpb.hh>
#include <test/util/deriv_funcs.hh>
#include <test/core/init_util.hh>


//Auto Headers
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <vector>


// --------------- Test Class --------------- //

// using declarations
using namespace core;
using namespace core::pose;
using namespace core::scoring;
using namespace core::scoring::methods;

class FullatomDisulfideEnergyTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void test_disulfE_deriv_check_w_total_flexibility()
	{
		core::pose::Pose pose = pdb1rpb_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( dslf_ss_dst, 0.5 );
		sfxn.set_weight( dslf_cs_ang, 0.5 );
		sfxn.set_weight( dslf_ss_dih, 0.5 );
		sfxn.set_weight( dslf_ca_dih, 0.5 );

		kinematics::MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv( pose, sfxn, movemap );
		adv.simple_deriv_check( true, 1e-5 ); // 1st argument true to make sure that start_score == start_func

	}

	void test_disulfE_deriv_check_w_partial_flexibility()
	{
		core::pose::Pose pose = pdb1rpb_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( dslf_ss_dst, 0.5 );
		sfxn.set_weight( dslf_cs_ang, 0.5 );
		sfxn.set_weight( dslf_ss_dih, 0.5 );
		sfxn.set_weight( dslf_ca_dih, 0.5 );

		kinematics::MoveMap movemap;
		movemap.set_bb( 15, true ); // right in the middle
		AtomDerivValidator adv( pose, sfxn, movemap );
		adv.simple_deriv_check( true, 1e-5 ); // 1st argument true to make sure that start_score == start_func

	}

};


