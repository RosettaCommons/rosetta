// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/comparative_modeling/util.cxxtest.hh
/// @brief  test suite for comparative modeling utilities in
/// protocols/comparative_modeling/util.cc.
/// @author James Thompson

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>

#include <core/pose/Pose.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/id/SequenceMapping.hh>

#include <core/conformation/Residue.hh>

#include <protocols/comparative_modeling/ThreadingMover.hh>

//Auto Headers
#include <core/types.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/sequence/Sequence.fwd.hh>
#include <string>


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
		using core::import_pose::pose_from_file;
		using core::pose::make_pose_from_sequence;

		SequenceOP query( new Sequence( "LNE-DILILGCSAMGDEVLE-ESEFEPFIEEI-STKISGKKVALFG", "4fxn_", 0 ) );
		SequenceOP templ( new Sequence( "FEGFDLVLLGCSTWGDDSIELQDDFIPLFDSLEETGAQGRKVACFG", "1f4pA", 0 ) );
		SequenceAlignment align;
		align.add_sequence(query);
		align.add_sequence(templ);

		Pose template_pose;
		core::import_pose::pose_from_file( template_pose, "protocols/comparative_modeling/1f4pA.pdb" , core::import_pose::PDB_file);

		// test basic threading
		Pose query_pose;
		core::pose::make_pose_from_sequence( query_pose, query->ungapped_sequence(), "fa_standard" );
		ThreadingMover mover( align, template_pose );
		mover.apply(query_pose);

		Real const TOLERANCE( 1e-5 );
		std::string const atom_name("CA");
		core::id::SequenceMapping map( align.sequence_mapping(1,2) );
		for ( Size ii = 1; ii <= query_pose.size(); ++ii ) {
			if ( map[ii] != 0 ) {
				//Size const tmpl_ii( map[ii] );
				core::Vector tgt_x = query_pose.residue(ii).xyz(atom_name);
				core::Vector src_x = query_pose.residue(ii).xyz(atom_name);
				TS_ASSERT_DELTA(tgt_x[0],src_x[0],TOLERANCE);
				TS_ASSERT_DELTA(tgt_x[1],src_x[1],TOLERANCE);
				TS_ASSERT_DELTA(tgt_x[2],src_x[2],TOLERANCE);
			}
		}
	} // test_basic_threading

}; // class ThreadingMover_Tests
