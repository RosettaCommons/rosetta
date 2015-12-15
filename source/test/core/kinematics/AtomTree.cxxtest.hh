// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/kinematics/AtomTree.cxxtest.hh
/// @brief  test suite for core::kinematics::AtomTree.cc
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit Headers
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/tree/JumpAtom.hh>
#include <core/kinematics/tree/BondedAtom.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <test/util/pose_funcs.hh>

// Utility headers
#include <basic/Tracer.hh>

// C/C++
#include <iostream>

//Auto Headers
#include <utility/vector1.hh>

#ifdef SERIALIZATION
// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


using namespace core::kinematics;
using namespace core::kinematics::tree;

class AtomTreeTests : public CxxTest::TestSuite {
public:

	void setUp() {
		core_init();
	}

	void see_if_tree_topologies_match( AtomCOP at1, AtomCOP at2 ) {
		using utility::pointer::dynamic_pointer_cast;
		TS_ASSERT_EQUALS( at1->n_children(), at2->n_children() );
		TS_ASSERT_EQUALS( (bool) dynamic_pointer_cast<   JumpAtom const > (at1), (bool) dynamic_pointer_cast<   JumpAtom const > (at2) );
		TS_ASSERT_EQUALS( (bool) dynamic_pointer_cast< BondedAtom const > (at1), (bool) dynamic_pointer_cast< BondedAtom const > (at2) );
		if ( at1->n_children() == at2->n_children() ) {
			for ( core::Size ii = 0; ii < at1->n_children(); ++ii ) {
				see_if_tree_topologies_match( at1->child( ii ), at2->child( ii ) );
			}
		}
	}

	void test_atom_tree_assignment_operator() {
		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		AtomTree const & at1 = trpcage.atom_tree();
		AtomTree at2;
		at2 = at1;
		see_if_tree_topologies_match( at1.root(), at2.root() );
	}

	void test_serialize_atom_tree() {
		TS_ASSERT( true );
#ifdef SERIALIZATION
		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		AtomTree at1 = trpcage.atom_tree();

		std::ostringstream oss;
		{
			cereal::BinaryOutputArchive arch( oss );
			arch( at1 );
		}

		AtomTree at2;
		std::istringstream iss( oss.str() );
		{
			cereal::BinaryInputArchive arch( iss );
			arch( at2 );
		}

		see_if_tree_topologies_match( at1.root(), at2.root() );
#endif

	}

};

