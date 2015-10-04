// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/id/SequenceMapping.hh
/// @brief
/// @author James Thompson

// Test headers
#include <cxxtest/TestSuite.h>

#include <core/id/SequenceMapping.hh>

#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

#include <core/sequence/util.hh>
#include <core/conformation/signals/LengthEvent.hh>

#include <basic/Tracer.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/types.hh>
#include <core/conformation/signals/LengthEvent.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/sequence/DerivedSequenceMapping.fwd.hh>
#include <core/sequence/ScoringScheme.fwd.hh>
#include <core/sequence/SequenceAlignment.fwd.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.fwd.hh>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iosfwd>
#include <limits>
#include <string>
#include <vector>


static basic::Tracer TR("test.core.id.SequenceMapping");


class SequenceMappingTests : public CxxTest::TestSuite
{
private:

	utility::vector1<core::Size> startmap_;
	core::pose::PoseOP pose_;

public:
	SequenceMappingTests() {}

	// Shared initialization goes here.
	void setUp() {
		core_init();
		pose_ = core::pose::PoseOP( new core::pose::Pose);
		core::import_pose::pose_from_pdb( *pose_, "core/conformation/test_in.pdb" );

		startmap_.clear();
		// Non-1 start
		startmap_.push_back(2);
		startmap_.push_back(3);
		startmap_.push_back(4);
		startmap_.push_back(5);
		// Gap
		startmap_.push_back(9);
		startmap_.push_back(10);
		startmap_.push_back(11);
	}

	/// @brief test SequenceMapping observer interface
	void test_constructor() {
		core::id::SequenceMapping seqmap(startmap_);
		TR << "Mapping:\n" << seqmap.to_string() << std::endl;
		//TS_ASSERT_EQUALS( seqmap.size1(), 11 ); // size1 not set with this constructor
		TS_ASSERT_EQUALS( seqmap.size2(), 7 );
		TS_ASSERT_EQUALS( seqmap[1], 2 );
		TS_ASSERT_EQUALS( seqmap[4], 5 );
		TS_ASSERT_EQUALS( seqmap[5], 9 );
		TS_ASSERT_EQUALS( seqmap[7], 11);
		TS_ASSERT_EQUALS( seqmap[8], 0 ); // Off the end is zero
	}

	void test_reverse() {
		core::id::SequenceMapping seqmap(startmap_);
		seqmap.reverse();
		TR << "Mapping:\n" << seqmap.to_string() << std::endl;
		//TS_ASSERT_EQUALS( seqmap.size1(), 11 ); // size1 not set with this constructor
		TS_ASSERT_EQUALS( seqmap.size2(), 7 );
		TS_ASSERT_EQUALS( seqmap[2], 1 );
		TS_ASSERT_EQUALS( seqmap[5], 4 );
		TS_ASSERT_EQUALS( seqmap[9], 5 );
		TS_ASSERT_EQUALS( seqmap[11], 7);
		TS_ASSERT_EQUALS( seqmap[12], 0 ); // Off the end is zero
		TS_ASSERT_EQUALS( seqmap[1], 0 );
		TS_ASSERT_EQUALS( seqmap[6], 0 );
	}

	// Test the length event constructor
	// Map is from the before numbers to the after numbers
	void test_lengthevent_append() {
		// Add 4 residues after residue 31
		core::conformation::signals::LengthEvent le(&pose_->conformation(),  core::conformation::signals::LengthEvent::RESIDUE_APPEND, 31, 4, &pose_->residue(1) );
		core::id::SequenceMapping seqmap(le);
		TR << "Mapping LengthEvent APPEND:\n" << seqmap.to_string() << std::endl;
		TS_ASSERT_EQUALS( seqmap[30], 30 );
		TS_ASSERT_EQUALS( seqmap[31], 31 );
		// Residues 32, 33, 34, 35 inserted
		TS_ASSERT_EQUALS( seqmap[32], 36 );
		TS_ASSERT_EQUALS( seqmap[33], 37 );
		TS_ASSERT_EQUALS( seqmap[34], 38 );
		TS_ASSERT_EQUALS( seqmap[35], 39 );
		TS_ASSERT_EQUALS( seqmap[36], 40 );
		TS_ASSERT_EQUALS( seqmap[37], 41 );
	}

	// Test the length event constructor
	// Map is from the before numbers to the after numbers
	void test_lengthevent_prepend() {
		// Add 4 residues before residue 31
		core::conformation::signals::LengthEvent le(&pose_->conformation(),  core::conformation::signals::LengthEvent::RESIDUE_PREPEND, 31, 4, &pose_->residue(1) );
		core::id::SequenceMapping seqmap(le);
		TR << "Mapping LengthEvent PREPEND:\n" << seqmap.to_string() << std::endl;
		TS_ASSERT_EQUALS( seqmap[30], 30 );
		// Residues 31, 32, 33, 34 inserted
		TS_ASSERT_EQUALS( seqmap[31], 35 );
		TS_ASSERT_EQUALS( seqmap[32], 36 );
		TS_ASSERT_EQUALS( seqmap[33], 37 );
		TS_ASSERT_EQUALS( seqmap[34], 38 );
		TS_ASSERT_EQUALS( seqmap[35], 39 );
		TS_ASSERT_EQUALS( seqmap[36], 40 );
		TS_ASSERT_EQUALS( seqmap[37], 41 );
	}

	// Test the length event constructor
	// Map is from the before numbers to the after numbers
	void test_lengthevent_delete() {
		// Delete 4 residues starting at residue 31
		core::conformation::signals::LengthEvent le(&pose_->conformation(),  core::conformation::signals::LengthEvent::RESIDUE_DELETE, 31, -4, &pose_->residue(1) );
		core::id::SequenceMapping seqmap(le);
		TR << "Mapping LengthEvent DELETE:\n" << seqmap.to_string() << std::endl;
		TS_ASSERT_EQUALS( seqmap[30], 30 );
		TS_ASSERT_EQUALS( seqmap[31],  0 );
		TS_ASSERT_EQUALS( seqmap[32],  0 );
		TS_ASSERT_EQUALS( seqmap[33],  0 );
		TS_ASSERT_EQUALS( seqmap[34],  0 );
		TS_ASSERT_EQUALS( seqmap[35], 31 );
		TS_ASSERT_EQUALS( seqmap[36], 32 );
		TS_ASSERT_EQUALS( seqmap[37], 33 );
	}

	// Test the length event constructor
	// Map is from the before numbers to the after numbers
	void test_lengthevent_append_end() {
		// Add 4 residues at the very end of the pose
		core::Size size_before( pose_->conformation().size() - 4 ); // Conformation is the post-change structure
		core::conformation::signals::LengthEvent le(&pose_->conformation(),  core::conformation::signals::LengthEvent::RESIDUE_APPEND, size_before, 4, &pose_->residue(1) );
		core::id::SequenceMapping seqmap(le);
		TR << "Mapping LengthEvent APPEND end:\n" << seqmap.to_string() << std::endl;
		TS_ASSERT_EQUALS( seqmap[size_before -1], size_before-1 );
		TS_ASSERT_EQUALS( seqmap[size_before], size_before );
		TS_ASSERT_EQUALS( seqmap[size_before+1], 0 ); // didn't exist before, no correspondance after
	}

	// Test the length event constructor
	// Map is from the before numbers to the after numbers
	void test_lengthevent_prepend_end() {
		// Add 4 residues at the very end of the pose
		core::Size size_before( pose_->conformation().size() - 4 ); // Conformation is the post-change structure
		core::conformation::signals::LengthEvent le(&pose_->conformation(),  core::conformation::signals::LengthEvent::RESIDUE_PREPEND, size_before, 4, &pose_->residue(1) );
		core::id::SequenceMapping seqmap(le);
		TR << "Mapping LengthEvent PREPEND end:\n" << seqmap.to_string() << std::endl;
		TS_ASSERT_EQUALS( seqmap[size_before -1], size_before-1 );
		TS_ASSERT_EQUALS( seqmap[size_before], pose_->conformation().size() ); // Last residue is still last residue
		TS_ASSERT_EQUALS( seqmap[size_before+1], 0 ); // didn't exist before, no correspondance after
	}

	// Test the length event constructor
	// Map is from the before numbers to the after numbers
	void test_lengthevent_delete_end() {
		// Deleted 4 residues from the end
		core::Size start_pos( pose_->conformation().size() + 1 ); // Starting point is the now non-existant residue one past the end.
		core::conformation::signals::LengthEvent le(&pose_->conformation(),  core::conformation::signals::LengthEvent::RESIDUE_DELETE, start_pos, -4, &pose_->residue(1) );
		core::id::SequenceMapping seqmap(le);
		TR << "Mapping LengthEvent DELETE end:\n" << seqmap.to_string() << std::endl;
		TS_ASSERT_EQUALS( seqmap[ start_pos-1 ], pose_->conformation().size() ); // Last residue is in same location
		TS_ASSERT_EQUALS( seqmap[ start_pos ], 0 ); // Deleted residue is missing
		TS_ASSERT_EQUALS( seqmap[ start_pos +1 ], 0 ); // Deleted residue is missing
		TS_ASSERT_EQUALS( seqmap[ start_pos +2 ], 0 ); // Deleted residue is missing
		TS_ASSERT_EQUALS( seqmap[ start_pos +3 ], 0 ); // Deleted residue is missing
		TS_ASSERT_EQUALS( seqmap[ start_pos +4 ], 0 ); // Did not exist before, should not exist after
	}

	void test_transitive_mapping() {
		using namespace core::id;
		using namespace core::sequence;

		SequenceMapping map1, map2;
		map1.insert_aligned_residue_safe( 1, 1 );
		map1.insert_aligned_residue_safe( 2, 0 );
		map1.insert_aligned_residue_safe( 3, 2 );
		map1.insert_aligned_residue_safe( 4, 3 );
		map1.insert_aligned_residue_safe( 5, 5 );
		map1.insert_aligned_residue_safe( 6, 6 );

		map2.insert_aligned_residue_safe( 2, 3 );
		map2.insert_aligned_residue_safe( 3, 0 );
		map2.insert_aligned_residue_safe( 4, 4 );
		map2.insert_aligned_residue_safe( 5, 6 );
		map2.insert_aligned_residue_safe( 6, 7 );

		SequenceMapping trans_map = transitive_map( map1, map2 );
		TS_ASSERT( trans_map[1] == 0 );
		TS_ASSERT( trans_map[2] == 0 );
		TS_ASSERT( trans_map[3] == 3 );
		TS_ASSERT( trans_map[4] == 0 );
		TS_ASSERT( trans_map[5] == 6 );
		TS_ASSERT( trans_map[6] == 7 );
	} // test_transitive_map

};
