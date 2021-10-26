// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/kinematics/inverse/jump.cxxtest.hh
/// @brief  Utility functions for calculating jumps by knowing desired atom positions
/// @author Jack Maguire

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit Headers
#include <core/kinematics/inverse/jump.hh>
#include <core/kinematics/Jump.hh>
#include <core/id/AtomID.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <test/util/pose_funcs.hh>
#include <core/import_pose/import_pose.hh>

// Utility headers

// C/C++
#include <stdexcept>

//Auto Headers

#ifdef SERIALIZATION
// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


using namespace core::kinematics;
using namespace core::kinematics::inverse;

class InverseJumpTests : public CxxTest::TestSuite {
public:

	void setUp() {
		core_init();
	}

	void test_known_solution(){
		//start with a known solution, reset the jump, try and recapture the original jump

		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, "core/pack/4MRS__opm.pdb" , core::import_pose::PDB_file );

		////////////
		//check pose
		TS_ASSERT_EQUALS( pose.num_jump(), 1 );
		Jump const desired_jump = pose.jump( 1 );

		core::Size const first_resid_of_chain2 = pose.chain_sequence( 1 ).size() + 1;
		TS_ASSERT( pose.chain_sequence( 2 ).size() > 5 ); //5 is arbitrary, just want to make sure it's big

		///////////////////////
		//capture desired state
		AlignmentAtomArray desired_state;
		desired_state.atoms[ 0 ].set(
			core::id::AtomID( 2, first_resid_of_chain2 ),
			pose.conformation()
		);
		desired_state.atoms[ 1 ].set(
			core::id::AtomID( 3, first_resid_of_chain2 ),
			pose.conformation()
		);
		desired_state.atoms[ 2 ].set(
			core::id::AtomID( 2, first_resid_of_chain2+2 ),
			pose.conformation()
		);

		////////////
		//reset jump
		pose.set_jump( 1, Jump() );


		///////////////////////
		//recapture desired state
		Jump const new_jump = calculate_new_jump( pose.conformation(), 1, desired_state );

		//Vector is 0-indexed but Matrix is 1-indexed???
		for ( core::Size ii = 0; ii < 3; ++ii ) {
			TS_ASSERT_DELTA( desired_jump.get_translation()[ii], new_jump.get_translation()[ii], 0.01 );
			for ( core::Size jj = 1; jj <= 3; ++jj ) {
				TS_ASSERT_DELTA( desired_jump.get_rotation().row(jj)[ii], new_jump.get_rotation().row(jj)[ii], 0.01 );
			}
		}
	}

	void test_throw_upon_bad_input(){
		//repeating the previous test, but the first alignment atom is not downstream of the jump

		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, "core/pack/4MRS__opm.pdb" , core::import_pose::PDB_file );

		////////////
		//check pose
		TS_ASSERT_EQUALS( pose.num_jump(), 1 );
		Jump const desired_jump = pose.jump( 1 );

		core::Size const first_resid_of_chain2 = pose.chain_sequence( 1 ).size() + 1;

		///////////////////////
		//capture desired state
		AlignmentAtomArray desired_state;
		desired_state.atoms[ 0 ].set(
			core::id::AtomID( 2, 1 ), // THIS IS THE DIFFERENT LINE
			pose.conformation()
		);
		desired_state.atoms[ 1 ].set(
			core::id::AtomID( 3, first_resid_of_chain2 ),
			pose.conformation()
		);
		desired_state.atoms[ 2 ].set(
			core::id::AtomID( 2, first_resid_of_chain2+2 ),
			pose.conformation()
		);

		////////////
		//reset jump
		pose.set_jump( 1, Jump() );

		///////////////////////
		//recapture desired state
		bool thrown = false;
		try {
			calculate_new_jump(pose.conformation(),1,desired_state);
		} catch( ... ){
			thrown = true;
		}
		TS_ASSERT( thrown );
	}

};

