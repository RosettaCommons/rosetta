// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/loops/LoopsUtil.cxxtest.hh
/// @brief test suite for protocols/loops/util
/// @author Christopher Miles (cmiles@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
#include <test/util/pose_funcs.hh>

// Project headers
#include <core/types.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/id/NamedAtomID.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/util.hh>

#include <numeric/xyzVector.hh>

namespace {

using core::Size;
using core::pose::Pose;
using core::pose::PoseOP;
using protocols::loops::Loop;
using protocols::loops::Loops;

class LoopsUtilTest : public CxxTest::TestSuite {

private:
	core::pose::Pose query_pose_;

private:
	bool equal_torsions(const Pose& p1, const Pose& p2) {
		if ( p1.size() != p2.size() ) {
			return false;
		}

		for ( Size i = 1; i <= p1.size(); ++i ) {
			if ( p1.phi(i) != p2.phi(i) ) {
				return false;
			} else if ( p1.psi(i) != p2.psi(i) ) {
				return false;
			} else if ( p1.omega(i) != p2.omega(i) ) {
				return false;
			}
		}

		return true;
	}

public:
	void setUp() {
		protocols_init();
		query_pose_ = create_test_in_pdb_pose(); // 116 residues
	}

	void testSafeExtendLoopsAndIdealize() {
		PoseOP pose = core::import_pose::pose_from_file("protocols/loops/2GB3.pdb", core::import_pose::PDB_file);
		Pose other = *pose;

		Loops loops;
		protocols::loops::safe_set_extended_torsions_and_idealize_loops(loops, &other);
		TS_ASSERT(equal_torsions(*pose, other));

		loops.push_back(Loop(1,5));
		protocols::loops::safe_set_extended_torsions_and_idealize_loops(loops, &other);
		TS_ASSERT(!equal_torsions(*pose, other));
	}

	void test_loop_picking_aligned() {
		utility::vector1< core::Size > unaligned_residues;

		unaligned_residues.push_back(   1 );
		unaligned_residues.push_back(   2 );
		unaligned_residues.push_back(   3 );
		unaligned_residues.push_back(  10 );
		unaligned_residues.push_back(  11 );
		unaligned_residues.push_back( 116 );

		protocols::loops::LoopsOP myloops = protocols::loops::pick_loops_unaligned(
			query_pose_.size(),
			unaligned_residues,
			5 // min_loop_size
		);

		// expected loops are:
		// 1-5
		// 8-13
		// 112-116
		TS_ASSERT( myloops->size() == 3 );

		protocols::loops::Loops::const_iterator it = myloops->begin();
		TS_ASSERT( it->start() ==   1 );
		TS_ASSERT( it->stop()  ==   5 );
		++it;
		TS_ASSERT( it->start() ==   9 );
		TS_ASSERT( it->stop()  ==  13 );
		++it;
		TS_ASSERT( it->start() == 112 );
		TS_ASSERT( it->stop()  == 116 );
		++it;
	} // void test_loop_picking_aligned()

	void test_loop_picking_chainbreak() {
		utility::vector1< core::Size > unaligned_residues;

		using namespace protocols::loops;

		// translate a couple of residues to weird places to make picking
		// loops by chainbreak find something.
		core::id::NamedAtomID id( "CA", 81 );
		query_pose_.set_xyz( id, query_pose_.xyz(id) + 15.0 );

		LoopsOP myloops = pick_loops_chainbreak(
			query_pose_,
			5 // min_loop_size
		);

		// expected loops is for five residues around residue 81.
		//TR << myloops << std::endl;
		//for ( protocols::loops::Loops::const_iterator it = myloops.begin(),
		//    it_end = myloops.end();
		//   it != it_end; ++it
		//) {
		// TR << "start,stop = " << it->start() << "," << it->stop()
		//  << std::endl;
		//}

		TS_ASSERT( myloops->size() == 1 );

		protocols::loops::Loops::const_iterator it = myloops->begin();
		TS_ASSERT( it->start() == 78 );
		TS_ASSERT( it->stop()  == 83 );
		++it;

	} // test_save_and_restore
};

}  // anonymous namespace
