// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/id/SequenceMapping.hh
/// @brief
/// @author James Thompson

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
#include <test/UTracer.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/Sequence.fwd.hh>
#include <core/id/SequenceMapping.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/sequence/AnnotatedSequence.hh>


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
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/file/FileName.fwd.hh>
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


static basic::Tracer TR("test.core.sequence.SequenceUtil");

class SequenceUtilTests : public CxxTest::TestSuite {

public:

	SequenceUtilTests() {}

	// Shared initialization goes here.
	void setUp() {
		core_init_with_additional_options("-mute core.io.pdb");
	}

	void test_annotated_sequence() {
		using core::Size;
		using namespace core::sequence;

		AnnotatedSequence seq1( "L[LEU:NtermProteinFull]KHSISKSGFKQ[GLN:CtermProteinFull]M[MET:NtermProteinFull]ESKRNKNIRVTTPKRHIDIH[HIS:CtermProteinFull]");
		std::string control_sequence( "LKHSISKSGFKQMESKRNKNIRVTTPKRHIDIH" );

		TR.Info << seq1 << std::endl;
		TR.Info << seq1.one_letter_sequence() << std::endl;

		TS_ASSERT_EQUALS( seq1.one_letter_sequence(), control_sequence );
		TS_ASSERT( seq1.is_patched( 1 ) );
		TS_ASSERT( !seq1.is_patched( 5 ) );
		TS_ASSERT( seq1.is_patched( 12 ) );
		TS_ASSERT( seq1.is_patched( 13 ) );
		TS_ASSERT( !seq1.is_patched( 14 ) );
		TS_ASSERT_EQUALS( seq1.patch_str( 12 ), "GLN:CtermProteinFull" );
		TS_ASSERT_EQUALS( seq1.patch_str( 1 ), "LEU:NtermProteinFull" );
		TS_ASSERT_EQUALS( seq1.one_letter( 12 ), 'Q' );
	}


	void test_naive_alignment() {
		using core::Size;
		using namespace core::sequence;

		SequenceOP seq1( new Sequence(
			"NIRNIIINIMAHELSVINNHIKYINELFYKLDTNHNGSLSHREIYTVLASVGIKKWDINRILQALDINDRGNITYTEFMAGCYRWKNIESTFLKAAFNKIDKDEDGYISKSDIVSLVHDKVLDNNDIDNFFLSVHSIKKGIPREHIINKISFQEFKDYML",
			"seq1", 16
			) );
		SequenceOP seq2( new Sequence(
			"VINNHIKYINELFYKLDTNHNGSLSHREIYTVLASVGIKKWDINRILQALDINDRGNITYTEFMAGCYRWKNIESTFLKAAFNKIDKDEDGYISKSDIVSLVHDKVLDNNDIDNFFLSVHSIKKGIPREHIINKISFQEFKDYMLS",
			"seq2", 31
			) );
		SequenceAlignment aln = align_naive( seq1, seq2 );
		TS_ASSERT( aln.gapped_positions() == 0 );
	}

	void test_alignment_regen(
		core::sequence::SequenceAlignment & aln
	) {
		using core::Size;
		using core::id::SequenceMapping;
		using namespace core::sequence;
		SequenceOP seq1_copy = aln.sequence(1)->clone();
		SequenceOP seq2_copy = aln.sequence(2)->clone();
		seq1_copy->sequence( seq1_copy->ungapped_sequence() );
		seq2_copy->sequence( seq2_copy->ungapped_sequence() );

		SequenceMapping map = aln.sequence_mapping( 1, 2 );
		SequenceAlignment new_aln = mapping_to_alignment(
			map, seq1_copy, seq2_copy
		);

		TS_ASSERT( aln.size() == new_aln.size() );

		for ( Size ii = 1; ii <= aln.size(); ++ii ) {
			//std::cout << "comparing" << std::endl;
			//std::cout << aln.sequence(ii)->to_string() << std::endl;
			//std::cout << new_aln.sequence(ii)->to_string() << std::endl;
			TS_ASSERT_EQUALS(
				aln.sequence(ii)->sequence(), new_aln.sequence(ii)->sequence()
			);
		}
	}

	void test_mapping_into_sequence() {
		using namespace core::sequence;
		SequenceOP seq1( new Sequence(
			"MRLGDAAELCYNLTSSYLQIAAESDSIIAQTQRAINT--TKSILINETFPKWSPLNGEISFSYNGGKDCQVLLLLYLSCLWEYYIVKLSQSQFDGKFHRFPLTKLPTVFIDHDDTFKTLENFIEETSLRYSL----SLYESDRDK----------------------------CETMAEAFETFLQVFPETKAIVIGIRHTDPFGEHLKPIQKTDANWPDFYRLQPLLHWNLANIWSFLLYSNEPICELYRYGFTSLGNVEETLPNPHLRKDKNSTPLKLNFEWEIENRYKHNEVTKAEPIPIADEDLVKIENLHEDYYPGWYLVDDKLERAGRIKKK",
			"t395_", 1
			) );
		SequenceOP seq2( new Sequence(
			"MKTYHLNN-----------DIIVTQEQLDHWNEQLIKLETPQEIIAWSIVTFP----HLFQTTAFGLTGLVTIDMLSKLS-------------------EKYYMPELLFIDTLHHFPQTLTLKNEIEKKYYQPKNQTIHVYKPDGCESEADFASKYGDFLWEKDDDKYDYLAKVEPAHRAYKEL-----HISAVFTGRRKSQGSARSQLSIIEIDE-LNGILKINPLINWTFEQVKQYIDANNVPYNELLDLGYRSIGDYHSTQPVK-------------EGEDERAGRW-------TECGIH------------------------EASRFAQF---",
			"2oq2A", 1
			) );

		SequenceAlignment aln;
		aln.add_sequence( seq2 );
		aln.add_sequence( seq1 );

		test_alignment_regen( aln );

		SequenceAlignment aln2;
		aln2.add_sequence( SequenceOP( new Sequence(
			"NIRVIARVRPVTKEDGEGPEATNAV------TFDADDDSI-------------IHLLHKGKPVSFELDKVFSPQASQQDVFQEVQALVTSCIDGFNVCIFAYGQTGAGKTYTMEGTAENPGIN------QRALQLLFSEVQEKAS-------DWEYTITVSAAEIYNEVLRDLLGKEPQEKL----------------EIRLCPDGSGQLYVPGLTEFQVQSVDDINKVFEFGHTNRTTEFTNLNEHSSRSHALLIVTVRGVDCSTG-----LRTTGKLNLVDLAGSERVGSRLR-------------EAQHINKSLSALGDVIAALRSRQGH--------VPFRNSKLTYLLQDSLSGDSKTLMVVQVSPVEKNTSETLYSLKFAER",
			"t313_", 1
			) ) );

		aln2.add_sequence( SequenceOP( new Sequence(
			"NIRVYCRIRPALKNLEN--------SDTSLINVN------EFDDNSGVQSMEVTKIQNTAQVHEFKFDKIFDQQDTNVDVFKEVGQLVQSSLDGYNVCIFAYGQTGSGKTFTMLN--------PGDGIIPSTISHIFNW------INKLKTKGWDYKVNCEFIEIYNENIVDL---------LRSDNNNKEDTSIGLKHEIRHDQETKTTTITNVTSCKLESEEMVEIILKKANKLRSTASTASNEHSSRSHSIFIIHLSGS-----NAKTGAHSYGTLNLVDLAGSE-------RINVSQVVGDRLRETQNINKSLSCLGDVIHALG-----QPDSTKRHIPFRNSKLTYLLQYSLTGDSKTLMFVNISPSSSHINETLNSLRFASK",
			"1f9tA", 2
			) ) );

		//test_alignment_regen( aln2 );

		//std::cout << aln;
		//std::cout << new_aln;
	} // test_mapping_into_sequence

}; // SequenceUtilTests
