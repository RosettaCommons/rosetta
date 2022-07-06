// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rosetta_scripts/ParsedProtocol.cxxtest.hh
/// @brief test suite for protocols::rosetta_scripts::ParsedProtocol
/// @author Jack Maguire, jackmaguire1444@gmail.com

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Package headers

#include <protocols/filters/Filter.hh>
#include <protocols/rosetta_scripts/ParsedProtocol.hh>
#include <protocols/moves/NullMover.hh>

#include <core/pose/annotated_sequence.hh>

// Basic headers

// Utility headers
#include <utility/exit.hh>
#include <basic/Tracer.hh>

// Numberic headers

// C++ headers

#include <core/pose/Pose.hh> // AUTO IWYU For Pose

static basic::Tracer TR("protocols.rosetta_scripts.ParsedProtocol.cxxtest");

using namespace protocols;
using namespace protocols::moves;
using namespace protocols::filters;

class ReportCounter : public protocols::filters::Filter {
	//Not 100% best practices here...
public:
	ReportCounter( core::Size * report_count ):
		report_count_( report_count )
	{}

	void report( std::ostream & os, core::pose::Pose const & ) const override {
		++(*report_count_);
		os << "report(): Filter reports " << (*report_count_) << " calls to this point." << std::endl;
	}

	core::Real report_sm( core::pose::Pose const & ) const override {
		++(*report_count_);
		TR << "report_sm(): Filter reports " << (*report_count_) << " calls to this point." << std::endl;
		return 0.0;
	}

	bool apply( core::pose::Pose const & ) const override {
		++(*report_count_);
		TR << "apply(): Filter reports " << (*report_count_) << " calls to this point." << std::endl;
		return true;
	}

	core::Real score( core::pose::Pose & ) override {
		++(*report_count_);
		TR << "score(): Filter reports " << (*report_count_) << " calls to this point." << std::endl;
		return 0;
	}


	FilterOP clone() const override {
		return FilterOP( new ReportCounter( *this ) );
	};

	FilterOP fresh_instance() const override {
		//return FilterOP( new ReportCounter() );
		runtime_assert( false );
		return nullptr;
	};

	//mutable core::Size report_count_ = 0;
	mutable core::Size * report_count_ = nullptr;
};

class ParsedProtocolTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init_with_additional_options( "-mute protocols.rosetta_scripts.ParsedProtocol.REPORT" );
	}

	void test_filter_reporting() {

		using namespace protocols::rosetta_scripts;
		using PP = ParsedProtocol; //for scoping

		core::Size step1_count = 0;
		core::Size step2_count = 0;
		core::Size step3_count = 0;

		PP::ParsedProtocolStep const step1(
			NullMoverOP( new NullMover() ), //mover
			"step1", //name
			FilterOP( new ReportCounter( &step1_count ) ), //filter
			PP::FilterReportTime::AT_END
		);

		PP::ParsedProtocolStep const step2(
			NullMoverOP( new NullMover() ), //mover
			"step2", //name
			FilterOP( new ReportCounter( &step2_count ) ), //filter
			PP::FilterReportTime::AFTER_APPLY
		);

		PP::ParsedProtocolStep const step3(
			NullMoverOP( new NullMover() ), //mover
			"step3", //name
			FilterOP( new ReportCounter( &step3_count ) ), //filter
			PP::FilterReportTime::NONE
		);

		ParsedProtocol pp;
		pp.add_step( step1 );
		pp.add_step( step2 );
		pp.add_step( step3 );

		core::pose::Pose test_pose;
		core::pose::make_pose_from_sequence(test_pose, "MENTENAI", "fa_standard", false);
		pp.apply( test_pose );

		TS_ASSERT_EQUALS( step1_count, 2 ); //Apply + end report_sm + end report, though end report is suppressed since report tracer is muted.  Also wasteful.
		TS_ASSERT_EQUALS( step2_count, 2 ); //Apply + report_sm + report, though report is suppressed since report tracer is muted. Wasteful if you ask me
		TS_ASSERT_EQUALS( step3_count, 1 );
	}

	void test_filter_reporting_commandline() {

		using namespace protocols::rosetta_scripts;
		using PP = ParsedProtocol; //for scoping

		core::Size step1_count = 0;
		core::Size step2_count = 0;
		core::Size step3_count = 0;

		PP::ParsedProtocolStep const step1(
			NullMoverOP( new NullMover() ), //mover
			"step1", //name
			FilterOP( new ReportCounter( &step1_count ) ), //filter
			PP::FilterReportTime::AT_END,
			true
		);

		PP::ParsedProtocolStep const step2(
			NullMoverOP( new NullMover() ), //mover
			"step2", //name
			FilterOP( new ReportCounter( &step2_count ) ), //filter
			PP::FilterReportTime::AFTER_APPLY,
			true
		);

		PP::ParsedProtocolStep const step3(
			NullMoverOP( new NullMover() ), //mover
			"step3", //name
			FilterOP( new ReportCounter( &step3_count ) ), //filter
			PP::FilterReportTime::NONE,
			true
		);

		ParsedProtocol pp;
		pp.add_step( step1 );
		pp.add_step( step2 );
		pp.add_step( step3 );

		core::pose::Pose test_pose;
		core::pose::make_pose_from_sequence(test_pose, "MENTENAI", "fa_standard", false);
		pp.apply( test_pose );

		// commandline option should result in all of these effectively being PP::FilterReportTime::NONE
		TS_ASSERT_EQUALS( step1_count, 1 );
		TS_ASSERT_EQUALS( step2_count, 1 );
		TS_ASSERT_EQUALS( step3_count, 1 );
	}

};
