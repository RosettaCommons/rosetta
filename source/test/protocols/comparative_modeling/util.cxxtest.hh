// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/comparative_modeling/util.cxxtest.hh
/// @brief  test suite for comparative modeling utilities in
/// protocols/comparative_modeling/util.cc.
/// @author James Thompson

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

#include <basic/Tracer.hh>
#include <core/id/NamedAtomID.hh>

#include <utility/vector1.hh>

#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/comparative_modeling/util.hh>

#include <numeric/xyzVector.hh>

//Auto Headers


static basic::Tracer TR("test.protocols.comparative_modeling.util");

using namespace core;

class CM_Util_Tests : public CxxTest::TestSuite {
public:
	CM_Util_Tests() {}

	void setUp() {
		core_init();
		query_pose_ = create_test_in_pdb_pose(); // 116 residues
		//core::import_pose::pose_from_file( query_pose_, "core/io/test_in.pdb" , core::import_pose::PDB_file);
	}

	void tearDown() {}

	void test_loop_picking_aligned() {
		utility::vector1< core::Size > unaligned_residues;

		using namespace protocols::comparative_modeling;
		unaligned_residues.push_back(   1 );
		unaligned_residues.push_back(   2 );
		unaligned_residues.push_back(   3 );
		unaligned_residues.push_back(  10 );
		unaligned_residues.push_back(  11 );
		unaligned_residues.push_back( 116 );

		protocols::loops::LoopsOP myloops = pick_loops_unaligned(
			query_pose_.total_residue(),
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
		using namespace protocols::comparative_modeling;

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

private:
	core::pose::Pose query_pose_;
}; // class CM_Util_Tests
