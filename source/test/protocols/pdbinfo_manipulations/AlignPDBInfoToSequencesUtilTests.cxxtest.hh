// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/pdbinfo_manipulations/AlignPDBInfoToSequencesTestsUtil.cxxtest.hh
/// @brief  tests for AlignPDBInfoToSequencesUtil
/// @author Dan Farrell (danpf@uw.edu)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Core Headers
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/sequence/SimpleScoringScheme.hh>
#include <core/sequence/SWAligner.hh>

// Project Headers
#include <protocols/pdbinfo_manipulations/AlignPDBInfoToSequencesUtil.hh>

class AlignPDBInfoToSequencesUtilTests : public CxxTest::TestSuite {

public:

	void setUp(){
		core_init();
	}

	void tearDown(){
	}

	void test_pad_sequence() {
		std::string const full_pose_seq("FGHILMN");
		core::sequence::SequenceOP const pose_seq(
			utility::pointer::make_shared< core::sequence::Sequence >( full_pose_seq, "pose", 1 )
		);
		std::string const full_target_seq("ABCDEFGHIJKLMNOP");
		core::sequence::SequenceOP const target_seq(
			utility::pointer::make_shared< core::sequence::Sequence >( full_target_seq, "target", 1 )
		);

		core::sequence::SWAligner sw_align;
		core::sequence::ScoringSchemeOP ss(  utility::pointer::make_shared< core::sequence::SimpleScoringScheme >(120, -999, -100, 0) );

		core::sequence::SequenceAlignment const chainSeq_to_target(
			sw_align.align( pose_seq, target_seq, ss )
		);
		core::sequence::SequenceOP current_aln(chainSeq_to_target.sequence(1));
		// Start must be 1 for pad sequences to work.
		current_aln->start(1);
		core::sequence::SequenceOP target_aln(chainSeq_to_target.sequence(2));
		protocols::pdbinfo_manipulations::pad_sequences(full_pose_seq, current_aln, full_target_seq, target_aln);
		TS_ASSERT_EQUALS(current_aln->sequence(), "-----FGHI--LMN--");
		TS_ASSERT_EQUALS(target_seq->sequence(), full_target_seq);
	}
};
