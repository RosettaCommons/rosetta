// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/denovo_design/SheetConstraintGenerator.cxxtest.hh
/// @brief  test suite for protocols::denovo_design::SheetConstraintGenerator.cc
/// @author Tom Linsky (tlinsky@uw.edu)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
#include <test/util/pose_funcs.hh>

// Unit headers
#include <protocols/fldsgn/SheetConstraintGenerator.hh>

// Devel headers

// Protocol headers
#include <protocols/denovo_design/architects/StrandArchitect.hh>
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/SegmentPairing.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/components/StructureDataFactory.hh>
#include <protocols/denovo_design/movers/BuildDeNovoBackboneMover.hh>
#include <protocols/denovo_design/util.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>
#include <protocols/fldsgn/topology/StrandPairing.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/stream_util.hh>

// Boost headers
#include <boost/assign.hpp>
#include <boost/lexical_cast.hpp>

// C++ headers

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.SheetConstraintGeneratorTests.cxxtest" );

namespace mytest {

typedef utility::vector1< std::string > Motifs;

protocols::denovo_design::components::StructureDataOP
structuredata_from_motifs( Motifs const & motifs )
{
	using namespace protocols::denovo_design;
	using namespace protocols::denovo_design::components;

	StructureDataOP perm( new StructureData( "test" ) );
	core::Size segment_num = 1;
	Motifs::const_iterator prev = motifs.end();
	Motifs::const_iterator next = ++motifs.begin();
	for ( Motifs::const_iterator m=motifs.begin(); m!=motifs.end(); ++m, ++segment_num, ++next ) {
		std::string const segment_name = boost::lexical_cast< std::string >( segment_num );

		std::string lower = "";
		if ( prev != motifs.end() ) {
			lower = boost::lexical_cast< std::string >( segment_num - 1 );
		}
		std::string upper = "";
		if ( next != motifs.end() ) {
			upper = boost::lexical_cast< std::string >( segment_num + 1 );
		}

		Segment newsegment( segment_name );
		newsegment.parse_motif( *m );
		newsegment.set_movable_group( segment_num );
		newsegment.delete_lower_padding();
		newsegment.delete_upper_padding();
		newsegment.set_lower_segment( lower );
		newsegment.set_upper_segment( upper );
		perm->add_segment( newsegment );
		prev = m;
	}
	return perm;
}

class SheetConstraintGeneratorTest : public protocols::fldsgn::SheetConstraintGenerator {
public:
	std::string
	get_secstruct_( core::pose::Pose const & pose ) const { return get_secstruct(pose); }

	std::string
	get_strandpairings_( core::pose::Pose const & pose ) const { return get_strandpairings(pose); }
};

} // namespace mytest

// --------------- Test Class --------------- //
class SheetConstraintGeneratorTests : public CxxTest::TestSuite {
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


	void reset( protocols::denovo_design::components::StructureData & perm )
	{
		using protocols::denovo_design::architects::UP;
		using protocols::denovo_design::architects::DOWN;
		using protocols::denovo_design::components::StrandPairing;
		using protocols::denovo_design::components::StrandPairingOP;
		perm.clear_pairings();
		perm.add_pairing( StrandPairing( "2", "4", UP, DOWN, 0 ) );
		perm.set_data_str( "2", "paired_strands", ",4" );
		perm.set_data_str( "4", "paired_strands", "2," );
		perm.set_data_int( "2", "orientation", 1 );
		perm.set_data_int( "4", "orientation", 0 );
	}

	void test_user_secstruct() {
		core::pose::Pose pose = create_trpcage_ideal_pose();
		for ( core::Size resid=1; resid<=pose.size(); ++resid ) {
			pose.set_secstruct( resid, 'E' );
		}

		mytest::SheetConstraintGeneratorTest sheet_rcg;
		sheet_rcg.set_use_dssp( false );

		TS_ASSERT_EQUALS( sheet_rcg.get_secstruct_(pose), "EEEEEEEEEEEEEEEEEEEE" );

		sheet_rcg.set_use_dssp( true );
		TS_ASSERT_EQUALS( sheet_rcg.get_secstruct_(pose), "LHHHHHHLLLHHHHLLLLLL" );

		std::string const test_secstruct = "LEEEEEEEELLEEEEEEEEL";
		sheet_rcg.set_secstruct( test_secstruct );
		TS_ASSERT_EQUALS( sheet_rcg.get_secstruct_(pose), test_secstruct );

		sheet_rcg.set_use_dssp( false );
		TS_ASSERT_EQUALS( sheet_rcg.get_secstruct_(pose), test_secstruct );
	}

	void test_spair_strings() {
		using protocols::denovo_design::components::StrandPairing;
		using protocols::denovo_design::components::StructureData;
		using protocols::denovo_design::components::StructureDataOP;
		using protocols::denovo_design::components::StructureDataFactory;
		using protocols::denovo_design::architects::UP;
		using protocols::denovo_design::architects::DOWN;

		core::pose::Pose pose = create_trpcage_ideal_pose();

		mytest::Motifs const motifs = boost::assign::list_of
			("1LX")("8EB")("2LG")("8EB")("1LX");

		StructureDataOP perm = mytest::structuredata_from_motifs( motifs );
		TS_ASSERT( perm );
		StructureData & sd = *perm;
		StructureDataFactory::get_instance()->save_into_pose( pose, sd );

		mytest::SheetConstraintGeneratorTest sheet_rcg;
		sheet_rcg.set_secstruct( "LEEEEEEEELLEEEEEEEEL" );
		std::string const spair_str = "1-2.A.0";
		sheet_rcg.set_strand_pairs( spair_str );
		TS_ASSERT_EQUALS( sheet_rcg.get_strandpairings_(pose), spair_str );

		sheet_rcg.set_strand_pairs( "" );
		TS_ASSERT_EQUALS( sheet_rcg.get_strandpairings_(pose), "" );

		// add pairings
		sd.add_pairing( StrandPairing( "2", "4", UP, DOWN, 0 ) );
		StructureDataFactory::get_instance()->save_into_pose( pose, sd );
		TS_ASSERT_EQUALS( sheet_rcg.get_strandpairings_(pose), spair_str );
	}

	void test_bulge_strandpairing() {
		using protocols::fldsgn::topology::SS_Info2;
		using protocols::fldsgn::topology::StrandPairingSet;
		using protocols::fldsgn::topology::StrandPairings;
		using protocols::denovo_design::movers::SetPoseSecstructFromStructureDataMover;
		using namespace protocols::denovo_design;
		using namespace protocols::denovo_design::components;

		mytest::Motifs const motifs = boost::assign::list_of
			("1LX")("5EB")("2LG")("4EB")("1LX");

		TR << motifs << std::endl;
		StructureDataOP perm = mytest::structuredata_from_motifs( motifs );
		TS_ASSERT( perm );
		TS_ASSERT_THROWS_NOTHING( perm->check_consistency() );
		TR << "Created " << *perm << std::endl;

		// add pairing info
		reset( *perm );

		// add bulge
		utility::vector1< std::string > e1_abego = abego_vector( perm->segment( "2" ).abego() );
		core::Size const abegoidx = 3;
		e1_abego[ abegoidx ] = "A";
		perm->set_abego( "2", e1_abego );
		TS_ASSERT_EQUALS( perm->segment( "2" ).abego()[ abegoidx - 1 ], 'A' );

		// count bulges
		TS_ASSERT_EQUALS( count_bulges( *perm, "2" ), 1 );
		TS_ASSERT_EQUALS( count_bulges( *perm, "4" ), 0 );

		core::pose::PoseOP dummy = protocols::denovo_design::construct_dummy_pose( "VAL", perm->pose_length() );

		// cache into pose
		StructureDataFactory::get_instance()->save_into_pose( *dummy, *perm );

		StructureData const & dummy_sd = StructureDataFactory::get_instance()->get_from_pose( *dummy );
		TS_ASSERT_THROWS_NOTHING( dummy_sd.check_consistency() );
		SetPoseSecstructFromStructureDataMover().apply( *dummy );

		// create generator
		mytest::SheetConstraintGeneratorTest sheet_rcg;
		sheet_rcg.set_strand_pairs( "1-2.A.0" );
		std::string const secstruct = sheet_rcg.get_secstruct_( *dummy );
		TS_ASSERT_EQUALS( secstruct, "LEEEEELLEEEEL" );

		std::string spair_str = sheet_rcg.get_strandpairings_( *dummy );
		TS_ASSERT_EQUALS( spair_str, "1-2.A.0" );

		sheet_rcg.set_strand_pairs( "" );
		spair_str = sheet_rcg.get_strandpairings_( *dummy );
		TS_ASSERT_EQUALS( spair_str, "1-2.A.0" );
		TR << secstruct << " " << spair_str << " abego " << dummy_sd.abego() << std::endl;

		protocols::fldsgn::topology::SS_Info2_OP ssinfo( new SS_Info2( *dummy, secstruct ) );
		StrandPairingSet spairset( spair_str, ssinfo, abego_vector( dummy_sd.abego() ) );
		StrandPairings spairs = spairset.strand_pairings();
		TS_ASSERT_EQUALS( spairs.size(), 1 );
		for ( StrandPairings::const_iterator sp=spairs.begin(); sp!=spairs.end(); ++sp ) {
			TS_ASSERT_EQUALS( (*sp)->begin1(), 2 );
			TS_ASSERT_EQUALS( (*sp)->end1(), 6 );
			TS_ASSERT_EQUALS( (*sp)->begin2(), 12 );
			TS_ASSERT_EQUALS( (*sp)->end2(), 9 );
		}

		protocols::fldsgn::ResiduePairs const pairs = sheet_rcg.compute_residue_pairs( spairs );
		TS_ASSERT_EQUALS( pairs.size(), 4 );
		TR << pairs << std::endl;

		std::set< protocols::fldsgn::ResiduePair > const respairset( pairs.begin(), pairs.end() );
		protocols::fldsgn::ResiduePairs const correct_pairs = boost::assign::list_of
			( protocols::fldsgn::ResiduePair( 2, 12, 'A' ) )
			( protocols::fldsgn::ResiduePair( 3, 11, 'A' ) )
			( protocols::fldsgn::ResiduePair( 5, 10, 'A' ) )
			( protocols::fldsgn::ResiduePair( 6, 9, 'A' ) );
		for ( protocols::fldsgn::ResiduePairs::const_iterator p=correct_pairs.begin(); p!=correct_pairs.end(); ++p ) {
			TR << "Comparing " << *p << " with set" << std::endl;
			TS_ASSERT_DIFFERS( respairset.find( *p ), respairset.end() );
		}

		// generate constraints and apply to pose
		core::Size const csts_per_pair = 6;
		core::scoring::constraints::ConstraintCOPs const all_csts = sheet_rcg.apply( *dummy );
		TS_ASSERT_EQUALS( all_csts.size(), csts_per_pair * correct_pairs.size() );
		TR << "Dummy has " << all_csts.size() << " constraints." << std::endl;

		typedef std::pair< core::Size, core::Size > ResPair;
		std::map< ResPair, core::Size > res_counts;
		for ( core::scoring::constraints::ConstraintCOPs::const_iterator c=all_csts.begin(); c!=all_csts.end(); ++c ) {
			TS_ASSERT( *c );
			utility::vector1< core::Size > const residues = (*c)->residues();
			core::Size const nres = residues.size();
			TS_ASSERT_EQUALS( nres, 2 );
			ResPair pair( residues[1], residues[2] );
			if ( residues[1] > residues[2] ) {
				pair = std::make_pair( residues[2], residues[1] );
			}
			std::map< ResPair, core::Size >::iterator rinfo = res_counts.find( pair );
			if ( rinfo == res_counts.end() ) {
				res_counts[ pair ] = 1;
			} else {
				rinfo->second += 1;
			}
		}

		// should be 6 csts for each pair
		for ( std::map< ResPair, core::Size >::const_iterator pair=res_counts.begin(); pair!=res_counts.end(); ++pair ) {
			TS_ASSERT_EQUALS( pair->second, csts_per_pair );
		}

		// parallel
		StrandPairingSet spairset2( "1-2.P.0", ssinfo, abego_vector( dummy_sd.abego() ) );
		StrandPairings spairs2 = spairset2.strand_pairings();
		TS_ASSERT_EQUALS( spairs2.size(), 1 );
		for ( StrandPairings::const_iterator sp=spairs2.begin(); sp!=spairs2.end(); ++sp ) {
			TS_ASSERT_EQUALS( (*sp)->begin1(), 2 );
			TS_ASSERT_EQUALS( (*sp)->end1(), 6 );
			TS_ASSERT_EQUALS( (*sp)->begin2(), 9 );
			TS_ASSERT_EQUALS( (*sp)->end2(), 12 );
		}

		protocols::fldsgn::ResiduePairs const pairs2 = sheet_rcg.compute_residue_pairs( spairs2 );
		TS_ASSERT_EQUALS( pairs2.size(), 4 );
		TR << pairs2 << std::endl;

		std::set< protocols::fldsgn::ResiduePair > const respairset2( pairs2.begin(), pairs2.end() );
		protocols::fldsgn::ResiduePairs const correct_pairs2 = boost::assign::list_of
			( protocols::fldsgn::ResiduePair( 2, 9, 'P' ) )
			( protocols::fldsgn::ResiduePair( 3, 10, 'P' ) )
			( protocols::fldsgn::ResiduePair( 5, 11, 'P' ) )
			( protocols::fldsgn::ResiduePair( 6, 12, 'P' ) );
		for ( protocols::fldsgn::ResiduePairs::const_iterator p=correct_pairs2.begin(); p!=correct_pairs2.end(); ++p ) {
			TS_ASSERT_DIFFERS( respairset2.find( *p ), respairset2.end() );
		}

	}
};
