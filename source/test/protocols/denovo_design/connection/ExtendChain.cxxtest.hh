// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*- // vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/protocols/denovo_design/ExtendChain.cxxtest.hh
/// @brief  test suite for protocols::denovo_design::connection::ExtendChain
/// @author Tom Linsky (tlinsky@uw.edu)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Unit headers
#include <protocols/denovo_design/connection/ExtendChain.hh>

// Protocol headers
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/forge/methods/chainbreak_eval.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/moves/Mover.hh>

// Core headers
#include <core/io/pdb/file_data.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>

// Boost headers

// unit test utility functions
#include <protocols/denovo_design/test_utils.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.connection.ExtendChain.cxxtest" );

// --------------- Test Class --------------- //
class ExtendChainTests : public CxxTest::TestSuite {
public:

	// Shared initialization goes here.
	void setUp() {
		protocols_init();

		// set preserve header always
		basic::options::option[basic::options::OptionKeys::run::preserve_header].value(true);
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	core::scoring::ScoreFunctionOP create_scorefxn() const
	{
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction() );
		scorefxn->set_weight( core::scoring::vdw, 1.0 );
		scorefxn->set_weight( core::scoring::rg, 1.0 );
		scorefxn->set_weight( core::scoring::rama, 0.1 );
		scorefxn->set_weight( core::scoring::hs_pair, 1.0 );
		scorefxn->set_weight( core::scoring::ss_pair, 1.0 );
		scorefxn->set_weight( core::scoring::rsigma, 1.0 );
		scorefxn->set_weight( core::scoring::coordinate_constraint, 1.0 );
		return scorefxn;
	}

	// test connectchains - interaction with components
	void test_input()
	{
		using namespace protocols::denovo_design;
		using namespace protocols::denovo_design::components;
		using namespace protocols::denovo_design::connection;

		basic::datacache::DataMap datamap;
		protocols::moves::Movers_map movers;
		protocols::filters::Filters_map filters;
		core::pose::Pose input_pose;
		ExtendChain extend;

		// should throw RosettaScriptsOption exception due to both segments being missing
		std::stringstream xml;
		xml << "<ExtendChain name=extend motif=\"2LG-6EB-1LX\" />";
		utility::tag::TagOP tag = utility::tag::Tag::create( xml );
		try {
			extend.parse_my_tag( tag, datamap, filters, movers, input_pose );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_RosettaScriptsOption const & e ) {
			TR << "Successfully caught exception ";
			e.show( TR );
			TR.flush();
			TS_ASSERT( true );
		}

		std::stringstream xml2;
		xml2 << "<ExtendChain name=extend motif=\"2LG-6EB-1LX\" segment1=\"pose.1\" segment2=\"pose.2\" />";
		tag = utility::tag::Tag::create( xml2 );
		try {
			extend.parse_my_tag( tag, datamap, filters, movers, input_pose );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_RosettaScriptsOption const & e ) {
			TR << "Successfully caught exception ";
			e.show( TR );
			TR.flush();
			TS_ASSERT( true );
		}

		std::stringstream xml3;
		xml3 << "<ExtendChain name=extend motif=\"2LG-6EB-1LX\" segment2=\"pose.2\" />";
		tag = utility::tag::Tag::create( xml3 );
		try {
			extend.parse_my_tag( tag, datamap, filters, movers, input_pose );
			TS_ASSERT( true );
		} catch ( utility::excn::EXCN_RosettaScriptsOption const & e ) {
			TR << "Caught exception ";
			e.show( TR );
			TR.flush();
			TS_ASSERT( false );
		}

		std::stringstream xml4;
		xml4 << "<ExtendChain name=extend motif=\"2LG-6EB-1LX\" segment1=\"pose.2\" />";
		tag = utility::tag::Tag::create( xml4 );
		try {
			extend.parse_my_tag( tag, datamap, filters, movers, input_pose );
			TS_ASSERT( true );
		} catch ( utility::excn::EXCN_RosettaScriptsOption const & e ) {
			TR << "Caught exception ";
			e.show( TR );
			TR.flush();
			TS_ASSERT( false );
		}

// must have motif specified!
		std::stringstream xml5;
		xml5 << "<ExtendChain name=extend segment1=\"pose.2\" />";
		tag = utility::tag::Tag::create( xml5 );
		try {
			extend.parse_my_tag( tag, datamap, filters, movers, input_pose );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_RosettaScriptsOption const & e ) {
			TR << "Caught exception ";
			e.show( TR );
			TR.flush();
			TS_ASSERT( true );
		}

// must have motif specified, not length!
		std::stringstream xml6;
		xml6 << "<ExtendChain name=extend length=10 segment1=\"pose.2\" />";
		tag = utility::tag::Tag::create( xml6 );
		try {
			extend.parse_my_tag( tag, datamap, filters, movers, input_pose );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_RosettaScriptsOption const & e ) {
			TR << "Caught exception ";
			e.show( TR );
			TR.flush();
			TS_ASSERT( true );
		}
	}

	void test_append() {
		using namespace protocols::denovo_design;
		using namespace protocols::denovo_design::components;
		using namespace protocols::denovo_design::connection;

		basic::datacache::DataMap datamap;
		protocols::moves::Movers_map movers;
		protocols::filters::Filters_map filters;
		core::pose::Pose input_pose;
		ExtendChain extend;

		// should throw RosettaScriptsOption exception due to both segments being missing
		std::stringstream xml;
		xml << "<ExtendChain name=extend motif=\"2LG-6HA\" chain1=\"2\" />";
		utility::tag::TagOP tag = utility::tag::Tag::create( xml );
		extend.parse_my_tag( tag, datamap, filters, movers, input_pose );
		extend.set_do_remodel( false );

		core::io::pdb::build_pose_from_pdb_as_is( input_pose, "protocols/denovo_design/connection/test_pdbcomp_BridgeChains.pdb" );

		StructureDataOP sd = StructureData::create_from_pose( input_pose, "extend" );
		TS_ASSERT( sd );
		TS_ASSERT( sd->check_consistency() );
		StructureDataOP orig = sd->clone();

		extend.setup_permutation( *sd );
		TS_ASSERT( sd->has_segment( "extend" ) );
		TS_ASSERT_EQUALS( extend.lower_segment_id( *sd ), "sheet1.s2" );
		TS_ASSERT_EQUALS( extend.upper_segment_id( *sd ), "" );
		TS_ASSERT_EQUALS( extend.loop_lower( *sd ), "extend");
		TS_ASSERT_EQUALS( extend.loop_upper( *sd ), "" );
		TS_ASSERT_EQUALS( sd->segment( "extend" ).lower_segment(), "sheet1.s2" );
		TS_ASSERT_EQUALS( sd->segment( "extend" ).upper_segment(), "" );
		TS_ASSERT_EQUALS( sd->segment( "sheet1.s2" ).lower_segment(), "" );
		TS_ASSERT_EQUALS( sd->segment( "sheet1.s2" ).upper_segment(), "extend" );
		TS_ASSERT_EQUALS( sd->segment( "sheet1.s1" ).lower_segment(), "" );
		check_unwanted_movement( *orig, *sd );
		orig = sd->clone();

		extend.apply_permutation( *sd );
		TS_ASSERT( sd->check_consistency() );
		// show have same # of chains that we started with
		TS_ASSERT_EQUALS( orig->pose()->conformation().num_chains(), sd->pose()->conformation().num_chains() );
		// should have 8 more residues
		TS_ASSERT_EQUALS( input_pose.total_residue() + 8, sd->pose_length() );
		// should not have any linear chainbreak
		TS_ASSERT_DELTA( protocols::forge::methods::linear_chainbreak( *sd->pose(), sd->segment( "sheet1.s2" ).cterm_resi() ), 0.0, 1e-1 );
		check_unwanted_movement( *orig, *sd );
	}

	void test_prepend() {
		using namespace protocols::denovo_design;
		using namespace protocols::denovo_design::components;
		using namespace protocols::denovo_design::connection;

		basic::datacache::DataMap datamap;
		protocols::moves::Movers_map movers;
		protocols::filters::Filters_map filters;
		core::pose::Pose input_pose;
		ExtendChain extend;

		// should throw RosettaScriptsOption exception due to both segments being missing
		std::stringstream xml;
		xml << "<ExtendChain name=extend motif=\"2LG-6HA\" chain2=\"2\" />";
		utility::tag::TagOP tag = utility::tag::Tag::create( xml );
		extend.parse_my_tag( tag, datamap, filters, movers, input_pose );
		extend.set_do_remodel( false );

		core::io::pdb::build_pose_from_pdb_as_is( input_pose, "protocols/denovo_design/connection/test_pdbcomp_BridgeChains.pdb" );

		StructureDataOP sd = StructureData::create_from_pose( input_pose, "extend" );
		TS_ASSERT( sd );
		TS_ASSERT( sd->check_consistency() );
		StructureDataOP orig = sd->clone();

		extend.setup_permutation( *sd );
		TS_ASSERT( sd->has_segment( "extend" ) );
		TS_ASSERT_EQUALS( extend.upper_segment_id( *sd ), "sheet1.s2" );
		TS_ASSERT_EQUALS( extend.lower_segment_id( *sd ), "" );
		TS_ASSERT_EQUALS( extend.loop_upper( *sd ), "extend");
		TS_ASSERT_EQUALS( extend.loop_lower( *sd ), "" );
		TS_ASSERT_EQUALS( sd->segment( "extend" ).upper_segment(), "sheet1.s2" );
		TS_ASSERT_EQUALS( sd->segment( "extend" ).lower_segment(), "" );
		TS_ASSERT_EQUALS( sd->segment( "sheet1.s2" ).upper_segment(), "" );
		TS_ASSERT_EQUALS( sd->segment( "sheet1.s2" ).lower_segment(), "extend" );
		TS_ASSERT_EQUALS( sd->segment( "sheet1.s1" ).lower_segment(), "" );
		check_unwanted_movement( *orig, *sd );
		orig = sd->clone();

		extend.apply_permutation( *sd );
		TS_ASSERT( sd->check_consistency() );
		// show have same # of chains that we started with
		TS_ASSERT_EQUALS( orig->pose()->conformation().num_chains(), sd->pose()->conformation().num_chains() );
		// should have 8 more residues
		TS_ASSERT_EQUALS( input_pose.total_residue() + 8, sd->pose_length() );
		// should not have any linear chainbreak
		TS_ASSERT_DELTA( protocols::forge::methods::linear_chainbreak( *sd->pose(), sd->segment( "extend" ).cterm_resi() ), 0.0, 1e-1 );
		check_unwanted_movement( *orig, *sd );
	}

};

