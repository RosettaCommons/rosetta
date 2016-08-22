// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/protocols/denovo_design/SegmentPairingTests.cxxtest.hh
/// @brief  test suite for protocols::denovo_design::SegmentPairing.cc
/// @author Tom Linsky (tlinsky@uw.edu)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Unit headers
#include <protocols/denovo_design/components/SegmentPairing.hh>

// Protocol headers
#include <protocols/denovo_design/architects/StrandArchitect.hh>
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/SegmentPairing.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/components/StructureDataFactory.hh>
#include <protocols/denovo_design/util.hh>

// Core headers

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

// Boost headers
#include <boost/assign.hpp>
#include <boost/lexical_cast.hpp>

// C++ headers

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.SegmentPairingTests.cxxtest" );

typedef utility::vector1< std::string > Motifs;

// --------------- Test Class --------------- //
class SegmentPairingTests : public CxxTest::TestSuite {
private:
	typedef protocols::denovo_design::components::SegmentPairing SegmentPairing;
	typedef protocols::denovo_design::components::StructureData StructureData;
	typedef protocols::denovo_design::components::StructureDataFactory StructureDataFactory;
	typedef protocols::denovo_design::components::StructureDataOP StructureDataOP;

public:

	// Shared initialization goes here.
	void setUp() {
		// load params for ligand
		protocols_init();

		// set preserve header always for "tomponent"
		basic::options::option[basic::options::OptionKeys::run::preserve_header].value(true);
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_helix_pairing_str()
	{
		using protocols::denovo_design::architects::UP;
		using protocols::denovo_design::architects::DOWN;
		using protocols::denovo_design::components::HelixPairing;
		using protocols::denovo_design::components::StrandPairing;
		using protocols::denovo_design::components::HelixSheetPairing;

		std::string const motifs = "16HA-4LX-15HA-4LX-17HA";
		StructureData const orig = StructureDataFactory::get_instance()->create_from_motifs( motifs, "" );

		std::string const h1 = "H01";
		std::string const h2 = "H02";
		std::string const h3 = "H03";
		std::string const l1 = "L01";

		TS_ASSERT_EQUALS( SegmentPairing::get_helix_pairings( orig ), "" );

		StructureData sd = orig;
		sd.add_pairing( StrandPairing( h1, l1, UP, UP, 0 ) );
		sd.add_pairing( HelixPairing( h3, h1, true ) );
		sd.add_pairing( HelixSheetPairing( h2, h1, h3 ) );
		TS_ASSERT_EQUALS( SegmentPairing::get_helix_pairings( sd ), "1-3.P" );

		sd = orig;
		sd.add_pairing( HelixPairing( h1, h2, false ) );
		TS_ASSERT_EQUALS( SegmentPairing::get_helix_pairings( sd ), "1-2.A" );

		sd = orig;
		sd.add_pairing( HelixPairing( h1, h2, true ) );
		sd.add_pairing( StrandPairing( h1, l1, UP, UP, 0 ) );
		sd.add_pairing( HelixPairing( h1, h3, false ) );
		TS_ASSERT_EQUALS( SegmentPairing::get_helix_pairings( sd ), "1-2.P;1-3.A" );
	}

	void test_hss_triplet_str()
	{
		using protocols::denovo_design::components::HelixSheetPairing;
		std::string const motif_str = "6EB-2LX-15HA-2LX-6EB";
		StructureData const orig = StructureDataFactory::get_instance()->create_from_motifs( motif_str, "" );

		std::string const s1 = "E01";
		std::string const h1 = "H01";
		std::string const s2 = "E02";

		TS_ASSERT_EQUALS( SegmentPairing::get_hss_triplets( orig ), "" );

		StructureData sd = orig;
		sd.add_pairing( HelixSheetPairing( h1, s1, s2 ) );
		TS_ASSERT_EQUALS( SegmentPairing::get_hss_triplets( sd ), "1,1-2" );

		sd = orig;
		sd.add_pairing( HelixSheetPairing( h1, s2, s1 ) );
		TS_ASSERT_EQUALS( SegmentPairing::get_hss_triplets( sd ), "1,2-1" );

		sd = orig;
		sd.add_pairing( HelixSheetPairing( h1, s1, s2 ) );
		sd.add_pairing( HelixSheetPairing( h1, s2, s1 ) );
		TS_ASSERT_EQUALS( SegmentPairing::get_hss_triplets( sd ), "1,1-2;1,2-1" );
	}

	void test_strand_pairing() {
		using protocols::denovo_design::components::StrandPairing;
		using protocols::denovo_design::architects::UP;
		using protocols::denovo_design::architects::DOWN;

		std::string const motif_str = "1LX-5EB-2LX-6EB-2LX-6EB-2LX-4EB-1LX";
		StructureData const orig = StructureDataFactory::get_instance()->create_from_motifs( motif_str, "" );

		std::string const s1_id = "E01";
		std::string const s2_id = "E02";
		std::string const s3_id = "E03";
		std::string const s4_id = "E04";

		TS_ASSERT_EQUALS( SegmentPairing::get_strand_pairings( orig ), "" );

		// UP/UP negative shift
		StructureData sd = orig;
		sd.add_pairing( StrandPairing( s1_id, s2_id, UP, UP, -1 ) );
		TS_ASSERT_EQUALS( SegmentPairing::get_strand_pairings( sd ), "1-2.P.-1" );

		// swapped order -- given pairing is identical to the one above
		sd = orig;
		sd.add_pairing( StrandPairing( s2_id, s1_id, UP, UP, -1 ) );
		TS_ASSERT_EQUALS( SegmentPairing::get_strand_pairings( sd ), "1-2.P.1" );

		// UP/UP positive shift
		sd = orig;
		sd.add_pairing( StrandPairing( s3_id, s4_id, UP, UP, 1 ) );
		TS_ASSERT_EQUALS( SegmentPairing::get_strand_pairings( sd ), "3-4.P.1" );

		// swapped order
		sd = orig;
		sd.add_pairing( StrandPairing( s4_id, s3_id, UP, UP, 1 ) );
		TS_ASSERT_EQUALS( SegmentPairing::get_strand_pairings( sd ), "3-4.P.-1" );

		// DOWN/DOWN negative shift
		sd = orig;
		sd.add_pairing( StrandPairing( s1_id, s2_id, DOWN, DOWN, -1 ) );
		TS_ASSERT_EQUALS( SegmentPairing::get_strand_pairings( sd ), "1-2.P.0" );

		// swapped order
		sd = orig;
		sd.add_pairing( StrandPairing( s2_id, s1_id, DOWN, DOWN, -1 ) );
		TS_ASSERT_EQUALS( SegmentPairing::get_strand_pairings( sd ), "1-2.P.-2" );

		// DOWN/DOWN positive shift
		sd = orig;
		sd.add_pairing( StrandPairing( s3_id, s4_id, DOWN, DOWN, 1 ) );
		TS_ASSERT_EQUALS( SegmentPairing::get_strand_pairings( sd ), "3-4.P.1" );

		// swapped order
		sd = orig;
		sd.add_pairing( StrandPairing( s4_id, s3_id, DOWN, DOWN, 1 ) );
		TS_ASSERT_EQUALS( SegmentPairing::get_strand_pairings( sd ), "3-4.P.3" );

		// UP/DOWN negative shift
		sd = orig;
		sd.add_pairing( StrandPairing( s1_id, s2_id, UP, DOWN, -1 ) );
		TS_ASSERT_EQUALS( SegmentPairing::get_strand_pairings( sd ), "1-2.A.-1" );

		// swapped order
		sd = orig;
		sd.add_pairing( StrandPairing( s2_id, s1_id, DOWN, UP, -1 ) );
		TS_ASSERT_EQUALS( SegmentPairing::get_strand_pairings( sd ), "1-2.A.1" );

		// UP/DOWN positive shift
		sd = orig;
		sd.add_pairing( StrandPairing( s3_id, s4_id, UP, DOWN, 1 ) );
		TS_ASSERT_EQUALS( SegmentPairing::get_strand_pairings( sd ), "3-4.A.1" );

		// swapped order
		sd = orig;
		sd.add_pairing( StrandPairing( s4_id, s3_id, DOWN, UP, 1 ) );
		TS_ASSERT_EQUALS( SegmentPairing::get_strand_pairings( sd ), "3-4.A.-1" );

		// DOWN/UP negative shift
		sd = orig;
		sd.add_pairing( StrandPairing( s1_id, s2_id, DOWN, UP, -1 ) );
		TS_ASSERT_EQUALS( SegmentPairing::get_strand_pairings( sd ), "1-2.A.0" );

		// swapped order
		sd = orig;
		sd.add_pairing( StrandPairing( s2_id, s1_id, UP, DOWN, -1 ) );
		TS_ASSERT_EQUALS( SegmentPairing::get_strand_pairings( sd ), "1-2.A.-2" );

		// UP/DOWN positive shift
		sd = orig;
		sd.add_pairing( StrandPairing( s3_id, s4_id, DOWN, UP, 1 ) );
		TS_ASSERT_EQUALS( SegmentPairing::get_strand_pairings( sd ), "3-4.A.1" );

		// swapped order
		sd = orig;
		sd.add_pairing( StrandPairing( s4_id, s3_id, DOWN, UP, 1 ) );
		TS_ASSERT_EQUALS( SegmentPairing::get_strand_pairings( sd ), "3-4.A.-1" );
	}

	void test_bulge_strandpairing() {
		using protocols::denovo_design::abego_vector;
		using protocols::denovo_design::count_bulges;
		using protocols::denovo_design::architects::UP;
		using protocols::denovo_design::architects::DOWN;
		using protocols::denovo_design::components::StrandPairing;

		std::string const motif_str = "1LX-5EB-2LG-4EB-1LX";
		StructureDataOP perm( new StructureData( StructureDataFactory::get_instance()->create_from_motifs( motif_str ) ) );
		TS_ASSERT( perm );
		TS_ASSERT_THROWS_NOTHING( perm->check_consistency() );

		std::string const s1 = "E01";
		std::string const s2 = "E02";

		// add pairing info
		perm->clear_pairings();
		perm->add_pairing( StrandPairing( s1, s2, UP, DOWN, 0 ) );

		// add bulge
		utility::vector1< std::string > e1_abego = abego_vector( perm->segment( s1 ).abego() );
		e1_abego[ 2 ] = "A";
		perm->set_abego( s1, e1_abego );
		TS_ASSERT_EQUALS( perm->segment( s1 ).abego()[ 1 ], 'A' );
		TR << *perm << std::endl;

		// count bulges
		TS_ASSERT_EQUALS( count_bulges( *perm, s1 ), 1 );
		TS_ASSERT_EQUALS( count_bulges( *perm, s2 ), 0 );

		// get strand_pairings
		TS_ASSERT_EQUALS( SegmentPairing::get_strand_pairings( *perm ), "1-2.A.0" );

		// switch to parallel
		perm->clear_pairings();
		perm->add_pairing( StrandPairing( s1, s2, UP, UP, 0 ) );
		TS_ASSERT_EQUALS( SegmentPairing::get_strand_pairings( *perm ), "1-2.P.0" );

		// flipped orientations
		perm->clear_pairings();
		perm->add_pairing( StrandPairing( s1, s2, DOWN, UP, 0 ) );
		TS_ASSERT_EQUALS( SegmentPairing::get_strand_pairings( *perm ), "1-2.A.0" );

		// flip again
		perm->clear_pairings();
		perm->add_pairing( StrandPairing( s1, s2, DOWN, DOWN, 0 ) );
		TS_ASSERT_EQUALS( SegmentPairing::get_strand_pairings( *perm ), "1-2.P.0" );

		// reset and try flipping the order
		perm->clear_pairings();
		perm->add_pairing( StrandPairing( s2, s1, DOWN, UP, 0 ) );
		TS_ASSERT_EQUALS( SegmentPairing::get_strand_pairings( *perm ), "1-2.A.0" );

		// switch to parallel
		perm->clear_pairings();
		perm->add_pairing( StrandPairing( s2, s1, UP, UP, 0 ) );
		TS_ASSERT_EQUALS( SegmentPairing::get_strand_pairings( *perm ), "1-2.P.0" );

		perm->clear_pairings();
		perm->add_pairing( StrandPairing( s2, s1, UP, DOWN, 0 ) );
		TS_ASSERT_EQUALS( SegmentPairing::get_strand_pairings( *perm ), "1-2.A.0" );

		perm->clear_pairings();
		perm->add_pairing( StrandPairing( s2, s1, DOWN, DOWN, 0 ) );
		TS_ASSERT_EQUALS( SegmentPairing::get_strand_pairings( *perm ), "1-2.P.0" );
	}
};
