// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
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

#include <test/core/init_util.hh>
#include <basic/Tracer.hh>
#include <test/UTracer.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/Sequence.fwd.hh>
#include <core/id/SequenceMapping.hh>
#include <core/sequence/SequenceAlignment.hh>

#include <numeric/random/random.hh>

//Auto Headers
#include <utility/stream_util.hh>



static basic::Tracer TR("test.core.sequence.SequenceUtil");

class SequenceUtilTests : public CxxTest::TestSuite {

public:

SequenceUtilTests() {}

// Shared initialization goes here.
void setUp() {
	core_init();
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
	aln2.add_sequence( new Sequence(
		"NIRVIARVRPVTKEDGEGPEATNAV------TFDADDDSI-------------IHLLHKGKPVSFELDKVFSPQASQQDVFQEVQALVTSCIDGFNVCIFAYGQTGAGKTYTMEGTAENPGIN------QRALQLLFSEVQEKAS-------DWEYTITVSAAEIYNEVLRDLLGKEPQEKL----------------EIRLCPDGSGQLYVPGLTEFQVQSVDDINKVFEFGHTNRTTEFTNLNEHSSRSHALLIVTVRGVDCSTG-----LRTTGKLNLVDLAGSERVGSRLR-------------EAQHINKSLSALGDVIAALRSRQGH--------VPFRNSKLTYLLQDSLSGDSKTLMVVQVSPVEKNTSETLYSLKFAER",
		"t313_", 1
	) );

	aln2.add_sequence( new Sequence(
		"NIRVYCRIRPALKNLEN--------SDTSLINVN------EFDDNSGVQSMEVTKIQNTAQVHEFKFDKIFDQQDTNVDVFKEVGQLVQSSLDGYNVCIFAYGQTGSGKTFTMLN--------PGDGIIPSTISHIFNW------INKLKTKGWDYKVNCEFIEIYNENIVDL---------LRSDNNNKEDTSIGLKHEIRHDQETKTTTITNVTSCKLESEEMVEIILKKANKLRSTASTASNEHSSRSHSIFIIHLSGS-----NAKTGAHSYGTLNLVDLAGSE-------RINVSQVVGDRLRETQNINKSLSCLGDVIHALG-----QPDSTKRHIPFRNSKLTYLLQYSLTGDSKTLMFVNISPSSSHINETLNSLRFASK",
		"1f9tA", 2
	) );

	//test_alignment_regen( aln2 );

	//std::cout << aln;
	//std::cout << new_aln;
} // test_mapping_into_sequence

}; // SequenceUtilTests
