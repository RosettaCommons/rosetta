// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/denovo_design/ExtendChainMover.cxxtest.hh
/// @brief  test suite for protocols::denovo_design::connection::ExtendChainMover
/// @author Tom Linsky (tlinsky@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Unit headers
#include <protocols/denovo_design/movers/ExtendChainMover.hh>

// Protocol headers
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/components/StructureDataFactory.hh>
#include <protocols/denovo_design/util.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/moves/Mover.hh>

// Core headers
#include <core/io/pdb/build_pose_as_is.hh>
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

static basic::Tracer TR( "protocols.denovo_design.connection.ExtendChainMover.cxxtest" );

// --------------- Test Class --------------- //
class ExtendChainMoverTests : public CxxTest::TestSuite {
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
		using namespace protocols::denovo_design::movers;

		basic::datacache::DataMap datamap;
		protocols::moves::Movers_map movers;
		protocols::filters::Filters_map filters;
		core::pose::Pose input_pose;
		ExtendChainMover extend;

		// should throw RosettaScriptsOption exception due to both segments being missing
		std::stringstream xml;
		xml << "<ExtendChain name=extend motif=\"2LG-6EB-1LX\" />";
		utility::tag::TagOP tag = utility::tag::Tag::create( xml );
		TS_ASSERT_THROWS( extend.parse_my_tag( tag, datamap, filters, movers, input_pose ), utility::excn::EXCN_RosettaScriptsOption );

		std::stringstream xml2;
		xml2 << "<ExtendChain name=extend motif=\"2LG-6EB-1LX\" segment=\"pose.1\" chain=\"2\" />";
		tag = utility::tag::Tag::create( xml2 );
		extend = ExtendChainMover();
		TS_ASSERT_THROWS( extend.parse_my_tag( tag, datamap, filters, movers, input_pose ), utility::excn::EXCN_RosettaScriptsOption );

		std::stringstream xml3;
		xml3 << "<ExtendChain name=extend motif=\"2LG-6EB-1LX\" segment=\"pose.2\" />";
		tag = utility::tag::Tag::create( xml3 );
		extend = ExtendChainMover();
		TS_ASSERT_THROWS_NOTHING( extend.parse_my_tag( tag, datamap, filters, movers, input_pose ) );

		std::stringstream xml4;
		xml4 << "<ExtendChain name=extend motif=\"2LG-6EB-1LX\" chain=\"2\" />";
		tag = utility::tag::Tag::create( xml4 );
		extend = ExtendChainMover();
		TS_ASSERT_THROWS_NOTHING( extend.parse_my_tag( tag, datamap, filters, movers, input_pose ) );

		// must have motif specified!
		std::stringstream xml5;
		xml5 << "<ExtendChain name=extend segment=\"pose.2\" />";
		tag = utility::tag::Tag::create( xml5 );
		extend = ExtendChainMover();
		TS_ASSERT_THROWS( extend.parse_my_tag( tag, datamap, filters, movers, input_pose ), utility::excn::EXCN_RosettaScriptsOption );

		// must have motif specified, not length!
		std::stringstream xml6;
		xml6 << "<ExtendChain name=extend length=10 segment1=\"pose.2\" />";
		tag = utility::tag::Tag::create( xml6 );
		extend = ExtendChainMover();
		TS_ASSERT_THROWS( extend.parse_my_tag( tag, datamap, filters, movers, input_pose ), utility::excn::EXCN_RosettaScriptsOption );
	}

	void test_append() {
		using namespace protocols::denovo_design;
		using namespace protocols::denovo_design::components;
		using namespace protocols::denovo_design::movers;
		StructureDataFactory const & factory = *StructureDataFactory::get_instance();

		core::pose::Pose pose;
		basic::datacache::DataMap datamap;
		protocols::moves::Movers_map movers;
		protocols::filters::Filters_map filters;
		ExtendChainMover extend;

		// should throw RosettaScriptsOption exception due to both segments being missing
		std::string const ext_id = "extend";
		std::stringstream xml;
		xml << "<ExtendChain name=" << ext_id << " motif=\"2LG-6HA\" chain=\"2\" />";
		utility::tag::TagOP tag = utility::tag::Tag::create( xml );
		extend.parse_my_tag( tag, datamap, filters, movers, pose );
		extend.set_dry_run( true );

		core::io::pdb::build_pose_from_pdb_as_is( pose, "protocols/denovo_design/connection/test_pdbcomp_BridgeChains.pdb" );
		core::pose::Pose const input_pose = pose;
		StructureData const orig = factory.get_from_pose( pose );
		TS_ASSERT_THROWS_NOTHING( orig.check_pose_consistency( input_pose ) );

		extend.apply( pose );
		TS_ASSERT( factory.has_cached_data( pose ) );
		StructureData const sd = factory.get_from_pose( pose );

		TS_ASSERT( sd.has_segment( ext_id ) );
		TS_ASSERT_EQUALS( sd.segment( ext_id ).lower_segment(), "sheet1.s2" );
		TS_ASSERT_EQUALS( sd.segment( ext_id ).upper_segment(), "" );
		TS_ASSERT_EQUALS( sd.segment( "sheet1.s2" ).lower_segment(), "" );
		TS_ASSERT_EQUALS( sd.segment( "sheet1.s2" ).upper_segment(), ext_id );
		TS_ASSERT_THROWS_NOTHING( sd.check_pose_consistency( pose ) );

		check_unwanted_movement( orig, input_pose, sd, pose );

		// show have same # of chains that we started with
		TS_ASSERT_EQUALS( input_pose.conformation().num_chains(), pose.conformation().num_chains() );
		// should have 8 more residues
		TS_ASSERT_EQUALS( input_pose.size() + 8, pose.size() );
		// should not have any linear chainbreak
		TS_ASSERT_DELTA( linear_chainbreak( pose, sd.segment( "sheet1.s2" ).upper() ), 0.0, 1e-1 );
	}

	void test_prepend() {
		using namespace protocols::denovo_design;
		using namespace protocols::denovo_design::components;
		using namespace protocols::denovo_design::movers;

		basic::datacache::DataMap datamap;
		protocols::moves::Movers_map movers;
		protocols::filters::Filters_map filters;
		core::pose::Pose pose;
		ExtendChainMover extend;

		std::string const e_id = "extend";
		std::stringstream xml;
		xml << "<ExtendChain name=" << e_id << " motif=\"2LG-6HA\" prepend=\"1\" chain=\"2\" />";
		utility::tag::TagOP tag = utility::tag::Tag::create( xml );
		extend.parse_my_tag( tag, datamap, filters, movers, pose );
		extend.set_dry_run( true );

		core::io::pdb::build_pose_from_pdb_as_is( pose, "protocols/denovo_design/connection/test_pdbcomp_BridgeChains.pdb" );
		core::pose::Pose const input_pose = pose;

		StructureDataFactory const & factory = *StructureDataFactory::get_instance();
		StructureData const orig = factory.create_from_pose( input_pose );
		TS_ASSERT_THROWS_NOTHING( orig.check_pose_consistency( input_pose ) );

		extend.apply( pose );
		TS_ASSERT( factory.has_cached_data( pose ) );
		StructureData const sd = factory.get_from_pose( pose );
		TS_ASSERT_THROWS_NOTHING( sd.check_pose_consistency( pose ) );

		TS_ASSERT( sd.has_segment( e_id ) );
		TS_ASSERT_EQUALS( sd.segment( e_id ).upper_segment(), "sheet1.s2" );
		TS_ASSERT_EQUALS( sd.segment( e_id ).lower_segment(), "" );
		TS_ASSERT_EQUALS( sd.segment( "sheet1.s2" ).upper_segment(), "" );
		TS_ASSERT_EQUALS( sd.segment( "sheet1.s2" ).lower_segment(), e_id );
		TS_ASSERT_EQUALS( sd.segment( "sheet1.s1" ).lower_segment(), "" );

		check_unwanted_movement( orig, input_pose, sd, pose );
		// show have same # of chains that we started with
		TS_ASSERT_EQUALS( input_pose.conformation().num_chains(), pose.conformation().num_chains() );
		// should have 8 more residues
		TS_ASSERT_EQUALS( input_pose.size() + 8, pose.size() );
		// should not have any linear chainbreak
		TS_ASSERT_DELTA( linear_chainbreak( pose, sd.segment( e_id ).upper() ), 0.0, 1e-1 );
	}

};

