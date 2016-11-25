// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rosetta_scripts/MultiplePoseMover.cxxtest.hh
/// @brief test suite for protocols::rosetta_scripts::MultiplePoseMover
/// @author Luki Goldschmidt <lugo@uw.edu>

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Package headers
#include <core/pose/Pose.fwd.hh>

#include <protocols/filters/Filter.hh>
#include <protocols/filters/BasicFilters.hh>
#include <protocols/filters/FilterFactory.hh>
#include <protocols/moves/MoverFactory.hh>
#include <protocols/rosetta_scripts/MultiplePoseMover.hh>
#include <protocols/rosetta_scripts/ParsedProtocol.hh>
#include <protocols/rosetta_scripts/RosettaScriptsParser.hh>
#include <protocols/simple_moves/MutateResidue.hh>

// Utility headers
#include <util/pose_funcs.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/tag/Tag.hh>
#include <utility/exit.hh>
#include <basic/Tracer.hh>

// Numberic headers

// C++ headers
#include <string>

static THREAD_LOCAL basic::Tracer TR("protocols.rosetta_scripts.MultiplePoseMover.cxxtest");

////////////////////////////////////////////////////////////////////////
// Dummy Mover that generates derivitives of poses from create_trpcage_ideal_pose()

#include <test/protocols/moves/DummyMover.hh>

class DummyMultipleOutputMover : public DummyMover {

public:
	DummyMultipleOutputMover() : npose_(0), pos_(1), max_poses_(100) {};

	protocols::moves::MoverOP clone() const {
		return protocols::moves::MoverOP( new DummyMultipleOutputMover(*this) );
	}

	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) {
		if ( tag->hasOption("pos") ) {
			pos_ = tag->getOption<int>("pos");
		}
		if ( tag->hasOption("max_poses") ) {
			max_poses_ = tag->getOption<int>("max_poses");
		}
	}

	void apply( core::pose::Pose & pose ) {
		base_pose_ = pose;
		generate_pose( pose );
	}

	core::pose::PoseOP get_additional_output() {
		if ( npose_ >= max_poses_ ) {
			return NULL;
		}
		if ( base_pose_.size() < 1 ) {
			base_pose_ = create_trpcage_ideal_pose();
		}
		core::pose::PoseOP p( new core::pose::Pose(base_pose_) );
		generate_pose( *p );
		return p;
	}

	void generate_pose( core::pose::Pose & pose ) {
		// Replace residue at position pos_ to one from set (in sequence)
		const char *sequence = get_sequence();
		int i = npose_ % strlen(sequence);

		protocols::simple_moves::MutateResidue mutate_mover(pos_, sequence[i]);
		mutate_mover.apply(pose);
		++npose_;
	}

	const char *get_sequence() {
		return "ACDEFGHIKLMNPQRSTVWY";
	}

private:
	int npose_;
	int pos_;
	int max_poses_;
	core::pose::Pose base_pose_;
};

class DummyMultipleOutputMoverCreator : public protocols::moves::MoverCreator {
public:
	protocols::moves::MoverOP create_mover() const override { return protocols::moves::MoverOP( new DummyMultipleOutputMover ); }
	std::string keyname() const override { return mover_name(); }
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const override {}
	static std::string mover_name() { return "DummyMultipleOutputMover"; }
};

class DummyFilter : public protocols::filters::TrueFilter { };

class DummyFilterCreator : public protocols::filters::FilterCreator {
public:
	protocols::filters::FilterOP create_filter() const override { return protocols::filters::FilterOP( new DummyFilter ); }
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const override {}
	std::string keyname() const override { return "DummyFilter"; }
};

class DummyHalfFilter : public protocols::filters::Filter {
	int mutable i_;
public:
	DummyHalfFilter() : Filter(), i_(0) { }
	bool apply( core::pose::Pose const & ) const {
		bool r = (bool)((++i_) % 2);
		TR << "DummyHalfFilter at " << i_ << ", returing: " << r << std::endl;
		return r;
	}
	protocols::filters::FilterOP clone() const { return protocols::filters::FilterOP( new DummyHalfFilter ); }
	protocols::filters::FilterOP fresh_instance() const { return protocols::filters::FilterOP( new DummyHalfFilter ); }
	protocols::filters::FilterOP create_filter() const { return protocols::filters::FilterOP( new DummyHalfFilter ); }
};

class DummyHalfFilterCreator : public protocols::filters::FilterCreator {
public:
	protocols::filters::FilterOP create_filter() const override { return protocols::filters::FilterOP( new DummyHalfFilter ); }
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const override {}
	std::string keyname() const override { return "DummyHalfFilter"; }
};

//using namespace basic::resource_manager;
using namespace protocols::moves;
using namespace protocols::rosetta_scripts;

////////////////////////////////////////////////////////////////////////
// Tests

class MultiplePoseMoverTests : public CxxTest::TestSuite {

private:
	// basic::resource_manager::FallbackConfigurationRegistrator< DummyResourceFallbackConfiguration1Creator > reg_fallback_dummy1;

public:

	void setUp() {
		protocols_init();

		static bool first_run = true;
		if ( first_run ) {
			using protocols::filters::FilterCreatorOP;
			protocols::moves::MoverFactory::get_instance()->factory_register( MoverCreatorOP( new DummyMultipleOutputMoverCreator ) );
			protocols::filters::FilterFactory::get_instance()->factory_register( FilterCreatorOP( new DummyFilterCreator ) );
			protocols::filters::FilterFactory::get_instance()->factory_register( FilterCreatorOP( new DummyHalfFilterCreator ) );
			first_run = false;
		}
	}

	void test_DummyMultipleOutputMover_get_additional_output() {

		/// Test get_additional_output() from DummyMultipleOutputMover,
		/// check provided pose sequence (first residue only)

		TR << "Testing DummyMultipleOutputMover::get_additional_output()" << std::endl;

		DummyMultipleOutputMover mover;
		const char *sequence = mover.get_sequence();

		for ( int i = 0; i < 20; ++i ) {
			core::pose::PoseOP pose = mover.get_additional_output();
			std::string pose_sequence( pose->sequence() );
			TR << "Pose " << (i+1) << " sequence: " << pose_sequence << std::endl;
			TS_ASSERT( pose_sequence[0] == sequence[i]);
		}

	}

	void test_tag_parent() {

		/// Test parent tag access

		TR << "Testing tag parent access" << std::endl;

		std::string xmlfile =
			"<ROSETTASCRIPTS>\n"
			"  <MOVERS>\n"
			"    <DummyMultipleOutputMover name=mover1/>\n"
			"    <MultiplePoseMover name=mover2 max_input_poses=10>\n"
			"      <ROSETTASCRIPTS>\n"
			"        <MOVERS>\n"
			"        </MOVERS>\n"
			"        <PROTOCOLS>\n"
			"          <Add mover_name=null/>\n"
			"        </PROTOCOLS>\n"
			"      </ROSETTASCRIPTS>\n"
			"    </MultiplePoseMover>\n"
			"  </MOVERS>\n"
			"  <PROTOCOLS>\n"
			"    <Add mover_name=mover1/>\n"
			"    <Add mover_name=mover2/>\n"
			"  </PROTOCOLS>\n"
			"</ROSETTASCRIPTS>\n";

		std::istringstream script_stream( xmlfile );
		utility::tag::TagCOP script_tags = utility::tag::Tag::create( script_stream );

		try {
			// Get second ROSETTASCRIPTS tag
			utility::tag::TagCOP tag1 = script_tags;
			utility::tag::TagCOP tag2 = tag1->getTag("MOVERS");
			utility::tag::TagCOP tag3 = tag2->getTag("MultiplePoseMover");
			utility::tag::TagCOP tag4 = tag3->getTag("ROSETTASCRIPTS");

			TS_ASSERT( tag4->getName() == "ROSETTASCRIPTS" && tag4 != tag1 );

			// Check parents
			TS_ASSERT( tag4->getParent().lock() == tag3 );
			TS_ASSERT( tag3->getParent().lock() == tag2 );
			TS_ASSERT( tag1->getParent().expired() );

		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Raised exception: " << e.msg() << std::endl;
			TS_ASSERT( false );
		}
	}

	void test_import() {

		/// Simple import of movers and filters from parent

		std::string xmlfile =
			"<ROSETTASCRIPTS>\n"
			"  <TASKOPERATIONS>\n"
			"    <InitializeFromCommandline name=init/>\n"
			"  </TASKOPERATIONS>\n"
			"  <FILTERS>\n"
			"    <DummyFilter name=filter1/>\n"
			"    <DummyFilter name=filter2/>\n"
			"  </FILTERS>\n"
			"  <MOVERS>\n"
			"    <DummyMultipleOutputMover name=mover1/>\n"
			"    <DummyMultipleOutputMover name=mover2/>\n"
			"    <MultiplePoseMover name=mpm>\n"
			"      <ROSETTASCRIPTS>\n"
			"        <IMPORT movers=mover2 filters=filter2,filter1 taskoperations=init/>\n"
			"        <MOVERS/>\n"
			"        <PROTOCOLS>\n"
			"          <Add mover_name=null/>\n"
			"        </PROTOCOLS>\n"
			"      </ROSETTASCRIPTS>\n"
			"    </MultiplePoseMover>\n"
			"  </MOVERS>\n"
			"  <PROTOCOLS>\n"
			"    <Add mover_name=mpm/>\n"
			"  </PROTOCOLS>\n"
			"</ROSETTASCRIPTS>\n";

		TR << "Testing IMPORT script tag using:\n" << xmlfile << std::endl;

		std::istringstream script_stream( xmlfile );
		utility::tag::TagCOP script_tags = utility::tag::Tag::create( script_stream );

		try {
			RosettaScriptsParser parser;
			MoverOP mover = parser.parse_protocol_tag( script_tags );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Raised exception: " << e.msg() << std::endl;
			TS_ASSERT( false );
		}

	}

	void test_import_multilevel() {

		/// Multi-level import, third level imports from level 2 and 1

		std::string xmlfile =
			"<ROSETTASCRIPTS>\n"
			"  <TASKOPERATIONS>\n"
			"    <InitializeFromCommandline name=init/>\n"
			"  </TASKOPERATIONS>\n"
			"  <FILTERS>\n"
			"    <DummyFilter name=filter1/>\n"
			"    <DummyFilter name=filter2/>\n"
			"  </FILTERS>\n"
			"  <MOVERS>\n"
			"    <DummyMultipleOutputMover name=mover1/>\n"
			"    <DummyMultipleOutputMover name=mover2/>\n"
			"    <MultiplePoseMover name=mpm1>\n"
			"      <ROSETTASCRIPTS>\n"
			"        <IMPORT movers=mover2 filters=filter1 taskoperations=init/>\n"
			"        <MOVERS>\n"
			"          <DummyMultipleOutputMover name=mover3/>\n"
			"          <MultiplePoseMover name=mpm2>\n"
			"            <ROSETTASCRIPTS>\n"
			"              <IMPORT movers=mover2 filters=filter2 taskoperations=init/>\n"
			"              <PROTOCOLS>\n"
			"                <Add mover_name=mover2 filter=filter2/>\n"
			"              </PROTOCOLS>\n"
			"            </ROSETTASCRIPTS>\n"
			"          </MultiplePoseMover>\n"
			"        </MOVERS>\n"
			"        <PROTOCOLS>\n"
			"          <Add mover_name=mpm2/>\n"
			"        </PROTOCOLS>\n"
			"      </ROSETTASCRIPTS>\n"
			"    </MultiplePoseMover>\n"
			"  </MOVERS>\n"
			"  <PROTOCOLS>\n"
			"    <Add mover_name=mover1/>\n"
			"    <Add mover_name=mpm1/>\n"
			"  </PROTOCOLS>\n"
			"</ROSETTASCRIPTS>\n";

		TR << "Testing multi-level IMPORT script tag using:\n" << xmlfile << std::endl;

		std::istringstream script_stream( xmlfile );
		utility::tag::TagCOP script_tags = utility::tag::Tag::create( script_stream );

		try {
			RosettaScriptsParser parser;
			MoverOP mover = parser.parse_protocol_tag( script_tags );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Raised exception: " << e.msg() << std::endl;
			TS_ASSERT( false );
		}

	}

	void test_input_limit() {

		/// Test input limit

		std::string xmlfile =
			"<ROSETTASCRIPTS>\n"
			"  <MOVERS>\n"
			"    <DummyMultipleOutputMover name=mover1/>\n"
			"    <MultiplePoseMover name=mover2 max_input_poses=5/>\n"
			"  </MOVERS>\n"
			"  <PROTOCOLS>\n"
			"    <Add mover_name=mover1/>\n"
			"    <Add mover_name=mover2/>\n"
			"  </PROTOCOLS>\n"
			"</ROSETTASCRIPTS>\n";

		TR << "Testing max_input_poses limit using:\n" << xmlfile << std::endl;

		std::istringstream script_stream( xmlfile );
		utility::tag::TagCOP script_tags = utility::tag::Tag::create( script_stream );

		try {
			RosettaScriptsParser parser;
			core::pose::Pose pose( create_trpcage_ideal_pose() );

			MoverOP mover = parser.parse_protocol_tag( pose, script_tags );
			mover->apply(pose);

			// Expecting 5 poses (1 + 4 additional ones)
			int i = 1;
			while ( mover->get_additional_output() ) ++i;
			TS_ASSERT_EQUALS( i, 5 );

		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Raised exception: " << e.msg() << std::endl;
			TS_ASSERT( false );
		}
	}

	void test_output_limit() {

		/// Test output limit

		std::string xmlfile =
			"<ROSETTASCRIPTS>\n"
			"  <MOVERS>\n"
			"    <DummyMultipleOutputMover name=mover1 max_poses=7/>\n"
			"    <MultiplePoseMover name=mover2 max_output_poses=5/>\n"
			"  </MOVERS>\n"
			"  <PROTOCOLS>\n"
			"    <Add mover_name=mover1/>\n"
			"    <Add mover_name=mover2/>\n"
			"  </PROTOCOLS>\n"
			"</ROSETTASCRIPTS>\n";

		TR << "Testing max_output_poses limit using:\n" << xmlfile << std::endl;

		std::istringstream script_stream( xmlfile );
		utility::tag::TagCOP script_tags = utility::tag::Tag::create( script_stream );

		try {
			RosettaScriptsParser parser;
			core::pose::Pose pose( create_trpcage_ideal_pose() );

			MoverOP mover = parser.parse_protocol_tag( pose, script_tags );
			mover->apply(pose);

			// Expecting 5 poses (1 + 4 additional ones)
			int i = 1;
			while ( mover->get_additional_output() ) ++i;
			TS_ASSERT_EQUALS( i, 5 );

		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Raised exception: " << e.msg() << std::endl;
			TS_ASSERT( false );
		}
	}

	void test_selector() {

		/// Test selector

		std::string xmlfile =
			"<ROSETTASCRIPTS>\n"
			"  <MOVERS>\n"
			"    <DummyMultipleOutputMover name=mover1/>\n"
			"    <MultiplePoseMover name=mover2 max_input_poses=10>\n"
			"      <SELECT>\n"
			"        <TopNByProperty n=5>\n"
			"          <EnergyReporter/>\n"
			"        </TopNByProperty>\n"
			"      </SELECT>\n"
			"      <ROSETTASCRIPTS>\n"
			"        <PROTOCOLS>\n"
			"          <Add mover_name=null/>\n"
			"        </PROTOCOLS>\n"
			"      </ROSETTASCRIPTS>\n"
			"    </MultiplePoseMover>\n"
			"  </MOVERS>\n"
			"  <PROTOCOLS>\n"
			"    <Add mover_name=mover1/>\n"
			"    <Add mover_name=mover2/>\n"
			"  </PROTOCOLS>\n"
			"</ROSETTASCRIPTS>\n";

		TR << "Testing SELECT script tag using:\n" << xmlfile << std::endl;

		std::istringstream script_stream( xmlfile );
		utility::tag::TagCOP script_tags = utility::tag::Tag::create( script_stream );

		try {
			RosettaScriptsParser parser;
			core::pose::Pose pose( create_trpcage_ideal_pose() );

			MoverOP mover = parser.parse_protocol_tag( pose, script_tags );
			mover->apply(pose);

			// Expecting 5 poses (1 + 4 additional ones) due to TopNByProperty n=10 above
			int i = 1;
			while ( mover->get_additional_output() ) ++i;
			TS_ASSERT( i == 5 );

		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Raised exception: " << e.msg() << std::endl;
			TS_ASSERT( false );
		}

	}

	void test_logical_selector_and() {

		/// Test AND selector

		std::string xmlfile =
			"<ROSETTASCRIPTS>\n"
			"  <MOVERS>\n"
			"    <DummyMultipleOutputMover name=mover1/>\n"
			"    <MultiplePoseMover name=mover2 max_input_poses=10>\n"
			"      <SELECT>\n"
			"        <AndSelector>\n"
			"          <TopNByProperty n=5>\n"
			"            <EnergyReporter/>\n"
			"          </TopNByProperty>\n"
			"          <TopNByProperty n=3>\n"
			"            <EnergyReporter/>\n"
			"          </TopNByProperty>\n"
			"        </AndSelector>\n"
			"      </SELECT>\n"
			"      <ROSETTASCRIPTS>\n"
			"        <PROTOCOLS>\n"
			"          <Add mover_name=null/>\n"
			"        </PROTOCOLS>\n"
			"      </ROSETTASCRIPTS>\n"
			"    </MultiplePoseMover>\n"
			"  </MOVERS>\n"
			"  <PROTOCOLS>\n"
			"    <Add mover_name=mover1/>\n"
			"    <Add mover_name=mover2/>\n"
			"  </PROTOCOLS>\n"
			"</ROSETTASCRIPTS>\n";

		TR << "Testing AndSelector using:\n" << xmlfile << std::endl;

		std::istringstream script_stream( xmlfile );
		utility::tag::TagCOP script_tags = utility::tag::Tag::create( script_stream );

		try {
			RosettaScriptsParser parser;
			core::pose::Pose pose( create_trpcage_ideal_pose() );

			MoverOP mover = parser.parse_protocol_tag( pose, script_tags );
			mover->apply(pose);

			// Expecting 3 poses (1 + 2 additional ones) due to n=3 above
			int i = 1;
			while ( mover->get_additional_output() ) ++i;
			TS_ASSERT( i == 3 );

		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Raised exception: " << e.msg() << std::endl;
			TS_ASSERT( false );
		}

	}

	void test_logical_selector_or_and() {

		/// Test OR/AND selector combo

		std::string xmlfile =
			"<ROSETTASCRIPTS>\n"
			"  <MOVERS>\n"
			"    <DummyMultipleOutputMover name=mover1/>\n"
			"    <MultiplePoseMover name=mover2 max_input_poses=10>\n"
			"      <SELECT>\n"
			"        <OrSelector>\n"
			"          <AndSelector>\n"
			"            <TopNByProperty n=7>\n"
			"              <EnergyReporter term=fa_atr/>\n"
			"            </TopNByProperty>\n"
			"            <TopNByProperty n=3>\n"
			"              <EnergyReporter term=fa_atr/>\n"
			"            </TopNByProperty>\n"
			"          </AndSelector>\n"
			"          <TopNByProperty n=5>\n"
			"            <EnergyReporter term=fa_atr/>\n"
			"          </TopNByProperty>\n"
			"        </OrSelector>\n"
			"      </SELECT>\n"
			"      <ROSETTASCRIPTS>\n"
			"        <PROTOCOLS>\n"
			"          <Add mover_name=null/>\n"
			"        </PROTOCOLS>\n"
			"      </ROSETTASCRIPTS>\n"
			"    </MultiplePoseMover>\n"
			"  </MOVERS>\n"
			"  <PROTOCOLS>\n"
			"    <Add mover_name=mover1/>\n"
			"    <Add mover_name=mover2/>\n"
			"  </PROTOCOLS>\n"
			"</ROSETTASCRIPTS>\n";

		TR << "Testing OrSelector / AndSelector using:\n" << xmlfile << std::endl;

		std::istringstream script_stream( xmlfile );
		utility::tag::TagCOP script_tags = utility::tag::Tag::create( script_stream );

		try {
			RosettaScriptsParser parser;
			core::pose::Pose pose( create_trpcage_ideal_pose() );

			MoverOP mover = parser.parse_protocol_tag( pose, script_tags );
			mover->apply(pose);

			// Expecting 5 poses (1 + 2 additional ones) due to (n=7 && n=3) || n=5 above
			int i = 1;
			while ( mover->get_additional_output() ) ++i;
			TS_ASSERT( i == 5 );

		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Raised exception: " << e.msg() << std::endl;
			TS_ASSERT( false );
		}

	}

	void test_selector_filter_reporter() {

		/// Test selector

		std::string xmlfile =
			"<ROSETTASCRIPTS>\n"
			"  <FILTERS>\n"
			"    <DummyFilter name=filter1/>\n"
			"  </FILTERS>\n"
			"  <MOVERS>\n"
			"    <DummyMultipleOutputMover name=mover1/>\n"
			"    <MultiplePoseMover name=mover2 max_input_poses=10>\n"
			"      <SELECT>\n"
			"        <TopNByProperty n=5>\n"
			"          <FilterReporter filter=filter1/>\n"
			"        </TopNByProperty>\n"
			"      </SELECT>\n"
			"      <ROSETTASCRIPTS>\n"
			"        <PROTOCOLS>\n"
			"          <Add mover_name=null/>\n"
			"        </PROTOCOLS>\n"
			"      </ROSETTASCRIPTS>\n"
			"    </MultiplePoseMover>\n"
			"  </MOVERS>\n"
			"  <PROTOCOLS>\n"
			"    <Add mover_name=mover1/>\n"
			"    <Add mover_name=mover2/>\n"
			"  </PROTOCOLS>\n"
			"</ROSETTASCRIPTS>\n";

		TR << "Testing FilterReporter using:\n" << xmlfile << std::endl;

		std::istringstream script_stream( xmlfile );
		utility::tag::TagCOP script_tags = utility::tag::Tag::create( script_stream );

		try {
			RosettaScriptsParser parser;
			core::pose::Pose pose( create_trpcage_ideal_pose() );

			MoverOP mover = parser.parse_protocol_tag( pose, script_tags );
			mover->apply(pose);

			// Expecting 5 poses (1 + 4 additional ones) due to TopNByProperty n=10 above
			int i = 1;
			while ( mover->get_additional_output() ) ++i;
			TS_ASSERT( i == 5 );

		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Raised exception: " << e.msg() << std::endl;
			TS_ASSERT( false );
		}

	}

	void test_selector_filter() {

		/// Test filter selector and max output
		/// This also tests uncached / on demand / lazy pull mode

		std::string xmlfile =
			"<ROSETTASCRIPTS>\n"
			"  <MOVERS>\n"
			"    <DummyMultipleOutputMover name=mover1 max_poses=10/>\n"
			"    <MultiplePoseMover name=mover2>\n"
			"      <SELECT>\n"
			"        <Filter>\n"
			"          <DummyHalfFilter name=filter1/>\n"
			"        </Filter>\n"
			"      </SELECT>\n"
			"      <ROSETTASCRIPTS>\n"
			"        <PROTOCOLS>\n"
			"          <Add mover_name=null/>\n"
			"        </PROTOCOLS>\n"
			"      </ROSETTASCRIPTS>\n"
			"    </MultiplePoseMover>\n"
			"  </MOVERS>\n"
			"  <PROTOCOLS>\n"
			"    <Add mover_name=mover1/>\n"
			"    <Add mover_name=mover2/>\n"
			"  </PROTOCOLS>\n"
			"</ROSETTASCRIPTS>\n";

		TR << "Testing Filter Selector using:\n" << xmlfile << std::endl;

		std::istringstream script_stream( xmlfile );
		utility::tag::TagCOP script_tags = utility::tag::Tag::create( script_stream );

		try {
			RosettaScriptsParser parser;
			core::pose::Pose pose( create_trpcage_ideal_pose() );

			MoverOP mover = parser.parse_protocol_tag( pose, script_tags );
			mover->apply(pose);

			// Expecting 5 poses (1 + 4 additional ones) due to max_poses=10
			// and DummyHalfFilter, which drops every other pose
			int i = 1;
			while ( mover->get_additional_output() ) ++i;
			TR << "i = " << i << std::endl;
			TS_ASSERT( i == 5 );

		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Raised exception: " << e.msg() << std::endl;
			TS_ASSERT( false );
		}

	}

	void test_invalid_import_mover() {

		/// Invalid import, mover2 is not defined

		std::string xmlfile =
			"<ROSETTASCRIPTS>\n"
			"  <FILTERS>\n"
			"    <DummyFilter name=filter1/>\n"
			"    <DummyFilter name=filter2/>\n"
			"  </FILTERS>\n"
			"  <MOVERS>\n"
			"    <DummyMultipleOutputMover name=mover1/>\n"
			"    <MultiplePoseMover name=mpm>\n"
			"      <ROSETTASCRIPTS>\n"
			"        <IMPORT movers=mover2/>\n"
			"        <MOVERS/>\n"
			"        <PROTOCOLS>\n"
			"          <Add mover_name=null/>\n"
			"        </PROTOCOLS>\n"
			"      </ROSETTASCRIPTS>\n"
			"    </MultiplePoseMover>\n"
			"  </MOVERS>\n"
			"  <PROTOCOLS>\n"
			"    <Add mover_name=mpm/>\n"
			"  </PROTOCOLS>\n"
			"</ROSETTASCRIPTS>\n";

		TR << "Testing invalid IMPORT script tag using:\n" << xmlfile << std::endl;

		std::istringstream script_stream( xmlfile );
		utility::tag::TagCOP script_tags = utility::tag::Tag::create( script_stream );

		try {
			RosettaScriptsParser parser;
			MoverOP mover = parser.parse_protocol_tag( script_tags );
		} catch ( utility::excn::EXCN_Msg_Exception const & e ) {
			std::string expected_error =
				"Exception in MultiplePoseMover with name \"mpm\": Failed to import mover2 from MOVERS";
			if ( e.msg() != expected_error ) {
				std::cout << e.msg() << std::endl;
				std::cout << expected_error << std::endl;
			}
			TS_ASSERT( e.msg() == expected_error );
		}
	}

	void test_invalid_import_mover_recursion() {

		/// Import mover into itself

		std::string xmlfile =
			"<ROSETTASCRIPTS>\n"
			"  <FILTERS/>\n"
			"  <MOVERS>\n"
			"    <MultiplePoseMover name=mpm>\n"
			"      <ROSETTASCRIPTS>\n"
			"        <IMPORT movers=mpm/>\n"
			"        <MOVERS/>\n"
			"        <PROTOCOLS>\n"
			"          <Add mover_name=null/>\n"
			"        </PROTOCOLS>\n"
			"      </ROSETTASCRIPTS>\n"
			"    </MultiplePoseMover>\n"
			"  </MOVERS>\n"
			"  <PROTOCOLS>\n"
			"    <Add mover_name=mpm/>\n"
			"  </PROTOCOLS>\n"
			"</ROSETTASCRIPTS>\n";

		TR << "Testing invalid IMPORT script tag with recursion using:\n" << xmlfile << std::endl;

		std::istringstream script_stream( xmlfile );
		utility::tag::TagCOP script_tags = utility::tag::Tag::create( script_stream );

		try {
			RosettaScriptsParser parser;
			MoverOP mover = parser.parse_protocol_tag( script_tags );
		} catch ( utility::excn::EXCN_Msg_Exception const & e ) {
			std::string expected_error =
				"Exception in MultiplePoseMover with name \"mpm\": Cannot import mover mpm into itself; recursion detected";
			if ( e.msg() != expected_error ) {
				std::cout << e.msg() << std::endl;
				std::cout << expected_error << std::endl;
			}
			TS_ASSERT( e.msg() == expected_error );
		}
	}

	void test_invalid_import_filter() {

		/// Invalid import, filter3 cannot be imported because it's in a sibling and not a parent

		std::string xmlfile =
			"<ROSETTASCRIPTS>\n"
			"  <FILTERS>\n"
			"    <DummyFilter name=filter1/>\n"
			"    <DummyFilter name=filter2/>\n"
			"  </FILTERS>\n"
			"  <MOVERS>\n"
			"    <DummyMultipleOutputMover name=mover1/>\n"
			"    <MultiplePoseMover name=mpm>\n"
			"      <ROSETTASCRIPTS>\n"
			"        <IMPORT filters=filter3/>\n"
			"        <MOVERS/>\n"
			"        <PROTOCOLS>\n"
			"          <Add mover_name=null/>\n"
			"        </PROTOCOLS>\n"
			"      </ROSETTASCRIPTS>\n"
			"    </MultiplePoseMover>\n"
			"    <MultiplePoseMover name=mpm2>\n"
			"      <ROSETTASCRIPTS>\n"
			"        <MOVERS/>\n"
			"        <FILTERS>\n"
			"          <DummyFilter name=filter3/>\n"
			"        </FILTERS>\n"
			"        <PROTOCOLS>\n"
			"          <Add mover_name=null/>\n"
			"        </PROTOCOLS>\n"
			"      </ROSETTASCRIPTS>\n"
			"    </MultiplePoseMover>\n"
			"  </MOVERS>\n"
			"  <PROTOCOLS>\n"
			"    <Add mover_name=mpm/>\n"
			"  </PROTOCOLS>\n"
			"</ROSETTASCRIPTS>\n";

		TR << "Testing invalid IMPORT script tag using:\n" << xmlfile << std::endl;

		std::istringstream script_stream( xmlfile );
		utility::tag::TagCOP script_tags = utility::tag::Tag::create( script_stream );

		try {
			RosettaScriptsParser parser;
			MoverOP mover = parser.parse_protocol_tag( script_tags );
		} catch ( utility::excn::EXCN_Msg_Exception const & e ) {
			std::string expected_error =
				"Exception in MultiplePoseMover with name \"mpm\": Failed to import filter3 from FILTERS";
			if ( e.msg() != expected_error ) {
				std::cout << e.msg() << std::endl;
				std::cout << expected_error << std::endl;
			}
			TS_ASSERT( e.msg() == expected_error );
		}
	}

	void test_mover_chaining() {

		/// Test mover chaining from DummyMultipleOutputMover -> MultiplePoseMover -> MultiplePoseMover

		std::string xmlfile =
			"<ROSETTASCRIPTS>\n"
			"  <MOVERS>\n"
			"    <DummyMultipleOutputMover name=mover1/>\n"
			"    <MultiplePoseMover name=mover2 max_input_poses=10>\n"
			"      <SELECT>\n"
			"        <TopNByProperty n=5>\n"
			"          <EnergyReporter/>\n"
			"        </TopNByProperty>\n"
			"      </SELECT>\n"
			"      <ROSETTASCRIPTS>\n"
			"        <PROTOCOLS>\n"
			"          <Add mover_name=null/>\n"
			"        </PROTOCOLS>\n"
			"      </ROSETTASCRIPTS>\n"
			"    </MultiplePoseMover>\n"
			"    <MultiplePoseMover name=mover3>\n"
			"      <ROSETTASCRIPTS>\n"
			"        <PROTOCOLS>\n"
			"          <Add mover_name=null/>\n"
			"        </PROTOCOLS>\n"
			"      </ROSETTASCRIPTS>\n"
			"    </MultiplePoseMover>\n"
			"  </MOVERS>\n"
			"  <PROTOCOLS>\n"
			"    <Add mover_name=mover1/>\n"
			"    <Add mover_name=mover2/>\n"
			"    <Add mover_name=mover3/>\n"
			"  </PROTOCOLS>\n"
			"</ROSETTASCRIPTS>\n";

		// This is an easy example with one distribute/collect branch,
		// and additional n:n pose passing (mover2 -> mover3)

		// mover1 generates poses ad infinitum
		// mover2 accepts up to 10 poses from mover1 (max_input_poses=10)
		//        selects 5 poses (TopNBy n=5)
		//        outputs 5 poses
		// mover3 accepts all poses from mover2 (5 total)
		//        applies null mover to 5 poses
		//        outputs 5 poses

		TR << "Testing mover chaining support using:\n" << xmlfile << std::endl;

		std::istringstream script_stream( xmlfile );
		utility::tag::TagCOP script_tags = utility::tag::Tag::create( script_stream );

		try {
			RosettaScriptsParser parser;
			core::pose::Pose pose( create_trpcage_ideal_pose() );

			MoverOP mover = parser.parse_protocol_tag( pose, script_tags );
			mover->apply(pose);

			// Expecting 5 poses (see comment above)
			int i = 1;
			core::pose::PoseOP p;
			std::set<std::string> sequences;

			TR << "Pose " << i << ": " << pose.sequence() << std::endl;
			sequences.insert(pose.sequence());

			while ( ( p = mover->get_additional_output() ) ) {
				std::string sequence(p->sequence());
				bool const exists( sequences.find( sequence ) != sequences.end() );
				if ( exists ) {
					TR << "Error: sequence present twice! => " << sequence << std::endl;
					continue;
				}
				++i;
				TR << "Pose " << i << ": " << sequence << std::endl;
				sequences.insert(sequence);
			}

			TS_ASSERT( i == 5 );

		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Raised exception: " << e.msg() << std::endl;
			TS_ASSERT( false );
		}

	}

	void test_mover_chaining_two_levels() {

		/// Test mover chaining from DummyMultipleOutputMover -> MultiplePoseMover ( DummyMultipleOutputMover ) -> MultiplePoseMover

		std::string xmlfile =
			"<ROSETTASCRIPTS>\n"
			"  <MOVERS>\n"
			"    <DummyMultipleOutputMover name=mover1a/>\n"
			"    <MultiplePoseMover name=mover1b max_input_poses=10>\n"
			"      <ROSETTASCRIPTS>\n"
			"        <MOVERS>\n"
			"          <DummyMultipleOutputMover name=mover2a pos=2/>\n"
			"          <MultiplePoseMover name=mover2b max_input_poses=5>\n"
			"            <ROSETTASCRIPTS>\n"
			"              <PROTOCOLS>\n"
			"                <Add mover_name=null/>\n"
			"              </PROTOCOLS>\n"
			"            </ROSETTASCRIPTS>\n"
			"          </MultiplePoseMover>\n"
			"        </MOVERS>\n"
			"        <PROTOCOLS>\n"
			"          <Add mover_name=mover2a/>\n"
			"          <Add mover_name=mover2b/>\n"
			"        </PROTOCOLS>\n"
			"      </ROSETTASCRIPTS>\n"
			"    </MultiplePoseMover>\n"
			"    <MultiplePoseMover name=mover1c>\n"
			"      <ROSETTASCRIPTS>\n"
			"        <PROTOCOLS>\n"
			"          <Add mover_name=null/>\n"
			"        </PROTOCOLS>\n"
			"      </ROSETTASCRIPTS>\n"
			"    </MultiplePoseMover>\n"
			"  </MOVERS>\n"
			"  <PROTOCOLS>\n"
			"    <Add mover_name=mover1a/>\n"
			"    <Add mover_name=mover1b/>\n"
			"  </PROTOCOLS>\n"
			"</ROSETTASCRIPTS>\n";

		// This is a fairly complex example with TWO distribute branches:
		// 1 -> 10, followed by 1 -> 5 ===> 1 -> 50 overall

		// mover1  generates poses ad infinitum
		// mover2  accepts up to 10 poses from mover1 (max_input_poses=10)
		//         feeds each pose to mover2b
		// mover2b calls mover2a on each of the 10 poses
		//         accepts 5 poses for each pose (50 total)
		//         outputs 50 poses, which mover2 passes on from mover2b
		// mover3  accepts all poses from mover2 (50 total)
		//         applies null mover to 50 poses
		//         outputs 50 poses

		TR << "Testing mover chaining support with 2 levels of MultiplePoseMover using:\n" << xmlfile << std::endl;

		std::istringstream script_stream( xmlfile );
		utility::tag::TagCOP script_tags = utility::tag::Tag::create( script_stream );

		try {
			RosettaScriptsParser parser;
			core::pose::Pose pose( create_trpcage_ideal_pose() );

			MoverOP mover = parser.parse_protocol_tag( pose, script_tags );
			mover->apply(pose);

			// Expecting 50 poses (see comment above)
			int i = 1;
			core::pose::PoseOP p;
			std::set<std::string> sequences;

			TR << "Pose " << i << ": " << pose.sequence() << std::endl;
			sequences.insert(pose.sequence());

			while ( ( p = mover->get_additional_output() ) ) {
				std::string sequence(p->sequence());
				bool const exists( sequences.find( sequence ) != sequences.end() );
				if ( exists ) {
					TR << "Error: sequence present twice! => " << sequence << std::endl;
					continue;
				}
				++i;
				TR << "Pose " << i << ": " << sequence << std::endl;
				sequences.insert(sequence);
			}

			TS_ASSERT( i == 50 );

		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Raised exception:: " << e.msg() << std::endl;
			TS_ASSERT( false );
		}

	}

	void test_mover_chaining_distribute_twice() {

		/// Test mover chaining from DummyMultipleOutputMover -> MultiplePoseMover ( DummyMultipleOutputMover ) -> MultiplePoseMover

		std::string xmlfile =
			"<ROSETTASCRIPTS>\n"
			"  <MOVERS>\n"
			"    <DummyMultipleOutputMover name=mover1 max_poses=10/>\n"
			"    <DummyMultipleOutputMover name=mover1b max_poses=10 pos=2/>\n"
			"    <MultiplePoseMover name=mover2 max_input_poses=3>\n"
			"      <ROSETTASCRIPTS>\n"
			"        <IMPORT movers=mover1b/>\n"
			"        <PROTOCOLS>\n"
			"          <Add mover_name=mover1b/>\n"
			"        </PROTOCOLS>\n"
			"      </ROSETTASCRIPTS>\n"
			"    </MultiplePoseMover>\n"
			"    <MultiplePoseMover name=mover3>\n"
			"      <ROSETTASCRIPTS>\n"
			"        <PROTOCOLS>\n"
			"          <Add mover_name=null/>\n"
			"        </PROTOCOLS>\n"
			"      </ROSETTASCRIPTS>\n"
			"    </MultiplePoseMover>\n"
			"  </MOVERS>\n"
			"  <PROTOCOLS>\n"
			"    <Add mover_name=mover1/>\n"
			"    <Add mover_name=mover2/>\n"
			"    <Add mover_name=mover3/>\n"
			"  </PROTOCOLS>\n"
			"</ROSETTASCRIPTS>\n";

		// This is a fairly complex example with TWO sequential distribute branches:
		// 1 -> 3, followed by 1 -> 10 ===> 1 -> 30 overall

		// mover1  generates poses up to 10 poses
		// mover2  accepts up to 3 poses from mover1 (max_input_poses=3)
		//         feeds each pose to another instance of mover1
		//         resulting in 10 output poses for each input pose
		// mover3  accepts all poses from mover2 (30 total)
		//         applies null mover to 30 poses
		//         outputs 30 poses

		TR << "Testing mover chaining support with 2 x distribute MultiplePoseMover using:\n" << xmlfile << std::endl;

		std::istringstream script_stream( xmlfile );
		utility::tag::TagCOP script_tags = utility::tag::Tag::create( script_stream );

		try {
			RosettaScriptsParser parser;
			core::pose::Pose pose( create_trpcage_ideal_pose() );

			MoverOP mover = parser.parse_protocol_tag( pose, script_tags );
			mover->apply(pose);

			// Expecting 30 poses (see comment above)
			int i = 1;
			core::pose::PoseOP p;
			std::set<std::string> sequences;

			TR << "Pose " << i << ": " << pose.sequence() << std::endl;
			sequences.insert(pose.sequence());

			while ( ( p = mover->get_additional_output() ) ) {
				std::string sequence(p->sequence());
				bool const exists( sequences.find( sequence ) != sequences.end() );
				if ( exists ) {
					TR << "Error: sequence present twice! => " << sequence << std::endl;
					continue;
				}
				++i;
				TR << "Pose " << i << ": " << sequence << std::endl;
				sequences.insert(sequence);
			}

			TS_ASSERT( i == 30 );

		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Raised exception:: " << e.msg() << std::endl;
			TS_ASSERT( false );
		}

	}

	void test_additionaloutputwrapper() {

		/// Test the MultipleOutputWrapper

		std::string xmlfile =
			"<ROSETTASCRIPTS>\n"
			"  <MOVERS>\n"
			"    <MultipleOutputWrapper name=mover1 max_output_poses=10>\n"
			"      <DummyMultipleOutputMover name=mover1b/>\n"
			"    </MultipleOutputWrapper>\n"
			"    <MultiplePoseMover name=mover2>\n"
			"      <ROSETTASCRIPTS>\n"
			"        <PROTOCOLS>\n"
			"          <Add mover_name=null/>\n"
			"        </PROTOCOLS>\n"
			"      </ROSETTASCRIPTS>\n"
			"    </MultiplePoseMover>\n"
			"  </MOVERS>\n"
			"  <PROTOCOLS>\n"
			"    <Add mover_name=mover1/>\n"
			"    <Add mover_name=mover2/>\n"
			"  </PROTOCOLS>\n"
			"</ROSETTASCRIPTS>\n";

		TR << "Testing additional output wrapper:\n" << xmlfile << std::endl;

		std::istringstream script_stream( xmlfile );
		utility::tag::TagCOP script_tags = utility::tag::Tag::create( script_stream );

		try {
			RosettaScriptsParser parser;
			core::pose::Pose pose( create_trpcage_ideal_pose() );

			MoverOP mover = parser.parse_protocol_tag( pose, script_tags );
			mover->apply(pose);

			// Expecting 10 poses, and all should be identical here since a fresh instance of
			// DummyMultipleOutputMover is used in each iteration of MultipleOutputWrapper
			int total_poses = 0, unique_poses = 1;
			core::pose::PoseOP p;
			std::set<std::string> sequences;

			++total_poses;
			TR << "Pose " << total_poses << ": " << pose.sequence() << std::endl;
			sequences.insert(pose.sequence());

			while ( ( p = mover->get_additional_output() ) ) {
				++total_poses;
				std::string sequence(p->sequence());
				TR << "Pose " << total_poses << ": " << pose.sequence() << std::endl;

				bool const exists( sequences.find( sequence ) != sequences.end() );
				if ( exists ) {
					continue;
				}
				++unique_poses;
				sequences.insert(sequence);
			}

			TS_ASSERT( unique_poses == 1 );
			TS_ASSERT( total_poses == 10 );

		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Raised exception:: " << e.msg() << std::endl;
			TS_ASSERT( false );
		}

	}

};
