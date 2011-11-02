// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
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

#include <test/core/init_util.hh>

#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
// AUTO-REMOVED #include <basic/Tracer.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/id/SequenceMapping.hh>

#include <core/conformation/Residue.hh>

#include <protocols/comparative_modeling/ThreadingMover.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <utility/vector1.hh>


class ThreadingMover_Tests : public CxxTest::TestSuite {

//static basic::Tracer TR("test.protocols.comparative_modeling.ThreadingMover");

public:
ThreadingMover_Tests() {}

void setUp() {
	core_init();
}

void tearDown() {}

void test_basic_threading() {
	using namespace core::sequence;
	using namespace protocols::comparative_modeling;
	using core::Size;
	using core::Real;
	using core::pose::Pose;
	using core::import_pose::pose_from_pdb;
	using core::pose::make_pose_from_sequence;

	SequenceOP query( new Sequence( "LNE-DILILGCSAMGDEVLE-ESEFEPFIEEI-STKISGKKVALFG", "4fxn_", '-' ) );
	SequenceOP templ( new Sequence( "FEGFDLVLLGCSTWGDDSIELQDDFIPLFDSLEETGAQGRKVACFG", "1f4pA", '-' ) );
	SequenceAlignment align;
	align.add_sequence(templ);
	align.add_sequence(query);

	Pose template_pose;
	core::import_pose::pose_from_pdb( template_pose, "protocols/comparative_modeling/1f4pA.pdb" );

	// test basic threading
	Pose query_pose;
	core::pose::make_pose_from_sequence( query_pose, query->ungapped_sequence(), "fa_standard" );
	ThreadingMover mover( align, template_pose );
	mover.apply(query_pose);

	Real const TOLERANCE( 1e-5 );
	std::string const atom_name("CA");
	core::id::SequenceMapping map( align.sequence_mapping(1,2) );
	for ( Size ii = 1; ii <= query_pose.total_residue(); ++ii ) {
		if ( map[ii] != 0 ) {
			//Size const tmpl_ii( map[ii] );
			TS_ASSERT_DELTA(
				query_pose.residue(ii).xyz(atom_name),
				template_pose.residue(ii).xyz(atom_name),
				TOLERANCE
			);
		}
	}
} // test_basic_threading

}; // class ThreadingMover_Tests
