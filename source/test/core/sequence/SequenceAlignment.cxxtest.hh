// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/sequence/SequenceAlignment.hh
/// @brief
/// @author James Thompson

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
#include <test/UTracer.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/id/SequenceMapping.hh>

#include <core/sequence/SWAligner.hh>
#include <core/sequence/ScoringScheme.fwd.hh>
#include <core/sequence/SimpleScoringScheme.hh>


//Auto Headers
#include <platform/types.hh>
#include <core/types.hh>
#include <core/conformation/signals/LengthEvent.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/sequence/Aligner.hh>
#include <core/sequence/DP_Matrix.hh>
#include <core/sequence/DerivedSequenceMapping.fwd.hh>
#include <core/sequence/ScoringScheme.hh>
#include <core/sequence/Sequence.fwd.hh>
#include <core/sequence/SequenceAlignment.fwd.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/file/FileName.fwd.hh>
#include <utility/io/izstream.fwd.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <numeric/xyzVector.fwd.hh>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iosfwd>
#include <limits>
#include <map>
#include <string>
#include <vector>
#include <boost/functional/hash.hpp>


static basic::Tracer TR("test.core.sequence.SequenceAlignment");

class SequenceAlignmentTests : public CxxTest::TestSuite {
public:
	SequenceAlignmentTests() {}

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	void test_alignment_io() {
		using namespace core::sequence;

		SequenceAlignment align;
		align.read_from_file( "core/sequence/test.aln" );

		SequenceOP seq1 = align.sequence(1);
		SequenceOP seq2 = align.sequence(2);
		TS_ASSERT( seq1->id() == "1yr0A.pdb" );
		TS_ASSERT( seq1->start() == 1 );

		TS_ASSERT( seq2->id() == "T0374" );
		TS_ASSERT( seq2->start() == 2 );
	} // test_alignment_io

	void test_grishin_align_io() {
		using namespace core::sequence;
		utility::vector1< SequenceAlignment > alignments
			= core::sequence::read_grishin_aln_file( "core/sequence/test.grishin_aln" );

		TS_ASSERT( alignments[1].score() == 175.0 );
		TS_ASSERT( alignments[2].score() == 175.0 );

		TS_ASSERT( alignments[2].sequence(1)->id() == "T0288" );
		TS_ASSERT( alignments[1].sequence(2)->id() == "2FNEA_1" );
		TS_ASSERT( alignments[2].sequence(1)->id() == "T0288" );
		TS_ASSERT( alignments[2].sequence(2)->id() == "2FNEA_2" );

		TS_ASSERT( alignments[1].sequence(1)->sequence() == "KVTLQKDAQNLIGISIGGG-----AQPCLYIVQVFDNTPAALDGTVAAGDEITGVNGRSIKGKTKVEVAKMIQEVKGEVTIHYNK" );
		TS_ASSERT( alignments[1].sequence(2)->sequence() == "SITLERGPDG-LGFSIVGGYGSPHGDLPIYVKTVFAKGAASEDGRLKRGDQIIAVNGQSLEGVTHEEAVAILKRTKGTVTLMVLS" );
		TS_ASSERT( alignments[2].sequence(1)->sequence() == "SMVP--GKVTLQKDAQNLIGISIGGG-----AQPCLYIVQVFDNTPAALDGTVAAGDEITGVNGRSIKGKTKVEVAKMIQEVKGEVTIHYNKLQYYKV" );
		TS_ASSERT( alignments[2].sequence(2)->sequence() == "--MPQCKSITLERGPDG-LGFSIVGGYGSPHGDLPIYVKTVFAKGAASEDGRLKRGDQIIAVNGQSLEGVTHEEAVAILKRTKGTVTLMVLSSDETSV" );
	} // test_grishin_aln_io

	void test_general_aln_io() {
		using namespace core::sequence;

		core::Real const TOLERATED_ERROR_SCORE( 1e-4 );
		utility::vector1< SequenceAlignment > alignments
			= core::sequence::read_general_aln_file( "core/sequence/test.general_aln" );

		TS_ASSERT( alignments.size() == 2 );

		// first alignment
		SequenceAlignment first( alignments[1] );
		TS_ASSERT_DELTA( first.score(), -108.9146, TOLERATED_ERROR_SCORE );

		TS_ASSERT( first.sequence(1)->id() == "1ub0A" );
		TS_ASSERT( first.sequence(2)->id() == "1jxhA" );
		TS_ASSERT( first.sequence(1)->start() == 4 );
		TS_ASSERT( first.sequence(2)->start() == 6 );
		TS_ASSERT( first.sequence(1)->sequence() == "ALTIAGSDSGGGAGVQADLKVFFRFGVYGTSALTLVTAQNTLGVQRVHLLPPEVVYAQIESVAQDFPLHAAKTGALGDAAIVEAVAEAVRRFGVRPLVVDPVM---AKEAAAALKERLFPLADLVTPNRLEAEALLGRP-IRTLKEAEEAAKALLALGPKAVLLKGGHLEAVDLLATRGGVLRFSAPRVHTRNTHGTGCTLSAAIAALLAKGRPLAEAVAEAKAYLTRALKTAPSL--GHGHGPLDHW" );
		TS_ASSERT( first.sequence(2)->sequence() == "ALTIAGTDPSGGAGIQADLKTFSALGAYGCSVITALVAENTCGVQSVYRIEPDFVAAQLDSVFSDVRIDTTKIGMLAETDIVEAVAERLQRHHVRNVVLDTVMLLLSPSAIETLRVRLLPQVSLITPNLPEAAALLDAPHARTEQEMLAQGRALLAMGCEAVLMKG------DWLFTREGEQRF---RVNTKNTHGTGCTLSAALAALRPRHRSWGETVNEAKAWLSAALAQADTLEVGKGIGPVHHF" );

		SequenceAlignment second( alignments[2] );
		TS_ASSERT_DELTA( second.score(), 2.0669, TOLERATED_ERROR_SCORE );

		TS_ASSERT( second.sequence(1)->id() == "1ub0A" );
		TS_ASSERT( second.sequence(2)->id() == "1w78A" );
		TS_ASSERT( second.sequence(1)->start() == 168  );
		TS_ASSERT( second.sequence(2)->start() == 27   );
		TS_ASSERT( second.sequence(1)->sequence() == "LEAVDLLATRGGVLRFSAPRVHT-RNTHGTGCTLSAAIAALLAKG" );
		TS_ASSERT( second.sequence(2)->sequence() == "LERVSLVAARLGVLK-PAPFVFTVAGTNGKGTTCRTLESILMAAG" );
	}

	void test_alignment_mapping() {
		using core::id::SequenceMapping;
		using namespace core::sequence;

		SequenceOP seq1( new Sequence( "ABCD-F", "one", 1 ) );
		SequenceOP seq2( new Sequence( "A-CDEF", "two", 1 ) );

		SequenceAlignment align;
		align.add_sequence( seq1 );
		align.add_sequence( seq2 );

		core::id::SequenceMapping mapping = align.sequence_mapping( 1, 2 );

		TS_ASSERT( mapping[ 1 ] == 1 );
		TS_ASSERT( mapping[ 2 ] == 0 );
		TS_ASSERT( mapping[ 3 ] == 2 );
		TS_ASSERT( mapping[ 4 ] == 3 );
		TS_ASSERT( mapping[ 5 ] == 5 );

		core::id::SequenceMapping reverse = align.sequence_mapping( 2, 1 );
		TS_ASSERT( reverse[ 1 ] == 1 );
		TS_ASSERT( reverse[ 2 ] == 3 );
		TS_ASSERT( reverse[ 3 ] == 4 );
		TS_ASSERT( reverse[ 4 ] == 0 );
		TS_ASSERT( reverse[ 5 ] == 5 );
	} // test_alignment_mapping

	void test_simple_aligner() {
		using namespace core::sequence;
		SequenceOP seq1( new Sequence(
			"PKALIVYGSTTGNTEYTAETIARELADAGYEVDSRDAASVEAGGLFEGFDLVLLGCSTWGDDSIELQDDFIPLFDSLEETGAQGRKVACFGCGDSSWEYFCGAVDAIEEKLKNLGAEIVQDGLRIDGDPRAARDDIVGWAHDVRGAI",
			"1f4pA_full", 1
			) );
		SequenceOP seq2( new Sequence(
			"FEGFDLVLLGCSTWGDDSIELQDDFIPLFDSLEETGAQGRKVACFG", "1f4pA_frag", 46
			) );

		SWAligner local_aligner;
		ScoringSchemeOP ss( new SimpleScoringScheme( 6, 1, -4, -1 ) );
		SequenceAlignment local_align = local_aligner.align( seq1, seq2, ss );

		TS_ASSERT( local_align.score() == 276 ); // 46 identities, each worth +6
		TS_ASSERT( local_align.sequence(1)->sequence() ==
			"FEGFDLVLLGCSTWGDDSIELQDDFIPLFDSLEETGAQGRKVACFG"
		);
		TS_ASSERT( local_align.sequence(2)->sequence() ==
			"FEGFDLVLLGCSTWGDDSIELQDDFIPLFDSLEETGAQGRKVACFG"
		);
		TS_ASSERT( local_align.sequence(1)->start() == 46 );
		TS_ASSERT( local_align.sequence(2)->start() == 46 );
		TS_ASSERT( local_align.sequence(1)->id() == "1f4pA_full" );
		TS_ASSERT( local_align.sequence(2)->id() == "1f4pA_frag" );
		TS_ASSERT( local_align.identities() == 46 );
	} // test_simple_aligner

	void test_alignment_functions() {
		using namespace core::sequence;
		SequenceOP seq1( new Sequence( "ABCDEFGHIJ", "first",  1 ) );
		SequenceOP seq2( new Sequence( "----EFG-IJ", "second", 1 ) );

		core::Real const arbitrary_score( 57.3 );

		SequenceAlignment aln;
		aln.add_sequence(seq2);
		aln.add_sequence(seq1);
		aln.score( arbitrary_score );
		TS_ASSERT( aln.score() == arbitrary_score );
		TS_ASSERT( aln.identities() == 5 );
		TS_ASSERT( aln.gapped_positions() == 5 );
		TS_ASSERT( aln.max_gap_percentage() == 0.5 );

		// test copying
		{
			SequenceAlignment aln_copy( aln );

			TS_ASSERT( aln.score() == aln_copy.score() );
			TS_ASSERT( aln.identities() == aln_copy.identities() );
			TS_ASSERT( aln.sequence(1)->sequence() == aln_copy.sequence(1)->sequence() );
			TS_ASSERT( aln.sequence(2)->sequence() == aln_copy.sequence(2)->sequence() );
		}
	} // test_alignment_functions

	void test_mapping_to_alignment() {
		using namespace core::sequence;
		SequenceOP seq1( new Sequence( "ABCDEGH", "first",  1 ) );
		SequenceOP seq2( new Sequence( "CDEFGH",  "second", 3 ) ); // This sequence's *third* residue is C, it's first two residues are undefined
		core::id::SequenceMapping map;
		map.insert_aligned_residue_safe( 3, 3 ); // C - C
		map.insert_aligned_residue_safe( 4, 4 ); // D - D
		map.insert_aligned_residue_safe( 5, 5 ); // E - E
		map.insert_aligned_residue_safe( 0, 6 ); // 0 - F
		map.insert_aligned_residue_safe( 6, 7 ); // G - G
		map.insert_aligned_residue_safe( 7, 8 ); // H - H
		TS_ASSERT( map[1] == 0 );
		TS_ASSERT( map[2] == 0 );
		TS_ASSERT( map[3] == 3 );
		TS_ASSERT( map[4] == 4 );
		TS_ASSERT( map[5] == 5 );
		TS_ASSERT( map[6] == 7 );
		TS_ASSERT( map[7] == 8 );
		SequenceAlignment align = mapping_to_alignment( map, seq1, seq2 );
		TS_ASSERT( align.sequence(1)->sequence() == "ABCDE-GH" );
		TS_ASSERT( align.sequence(2)->sequence() == "--CDEFGH" );
	}

	void test_per_residue_scores() {
		using namespace core::sequence;
		SequenceOP seq1( new Sequence( "ABCD--F", "one", 1 ) );
		SequenceOP seq2( new Sequence( "A-CDEHF", "two", 1 ) );

		SequenceAlignment align;
		align.add_sequence( seq1 );
		align.add_sequence( seq2 );
		ScoringSchemeOP ss( new SimpleScoringScheme( 6, 1, -4, -1 ) );

		using core::Real;
		using utility::vector1;
		vector1< Real > scores = align.calculate_per_position_scores( ss );
		TS_ASSERT( align.size()   == 2 );
		TS_ASSERT( align.length() == 7 );
		TS_ASSERT( scores.size() == align.length() );
		TS_ASSERT( scores[1] ==  6 );
		TS_ASSERT( scores[2] == -4 );
		TS_ASSERT( scores[3] ==  6 );
		TS_ASSERT( scores[4] ==  6 );
		TS_ASSERT( scores[5] == -4 );
		TS_ASSERT( scores[6] == -1 );
		TS_ASSERT( scores[7] ==  6 );
	} // test_per_residue_scores

	void test_equality_operator() {
		using namespace core::sequence;
		Sequence seq1( "ABCD--F", "one",  1 );
		Sequence seq2( "A-CDEHF", "two",  1 );
		Sequence seq3( "ABCD--F", "one",  3 );
		Sequence seq4( "ABCD--F", "four", 1 );
		Sequence seq5( "ABCD--F", "one",  1 );

		TS_ASSERT( !(seq1 == seq2) );
		TS_ASSERT( !(seq1 == seq3) );
		TS_ASSERT( !(seq1 == seq4) );
		TS_ASSERT(   seq1 == seq5 );
	}

}; // SequenceAlignmentTests
