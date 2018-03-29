// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/moves/FilterReportAsPoseExtraScoresMover.cxxtest.hh
/// @brief  test suite for protocols::moves::FilterReportAsPoseExtraScoresMover.cc
/// @author Steven Lewis smlewi@gmail.com


// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

// Unit headers
#include <protocols/moves/FilterReportAsPoseExtraScoresMover.hh>

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>

#include <protocols/pose_creation/MakePolyXMover.hh>
#include <protocols/score_filters/ScoreTypeFilter.hh>

// Package headers

#include <basic/Tracer.hh>

//Auto Headers

static basic::Tracer TR("protocols.moves.FilterReportAsPoseExtraScoresMover.cxxtest");

///////////////////////////////////////////////////////////////////////////
/// @name FilterReportAsPoseExtraScoresMoverTest
/// @brief: test FilterReportAsPoseExtraScoresMover (FRAPESM) - does a score get added to the Pose from a Filter?
/// @details Runs the FRAPESM with a scoring Filter, alters the Pose to change the score, reruns the FRAPESM; then checks that the Pose has the correct ExtraScores.
/// @author Steven Lewis July 2016
///////////////////////////////////////////////////////////////////////////
class FilterReportAsPoseExtraScoresMoverTests : public CxxTest::TestSuite {

public:

	//SF, scoretype, and filter can be shared
	core::scoring::ScoreFunctionOP sf;
	core::scoring::ScoreType ref;

	core::pose::Pose pose;

	protocols::moves::MoverOP polyx;
	protocols::filters::FilterOP scorefilter;

	std::string const before;
	std::string const after;

	FilterReportAsPoseExtraScoresMoverTests() :
		before("before"),
		after("after")
	{
		//TR << "CTOR" << std::endl;
	}

	void setUp() {
		//TR << "top setup" << std::endl;
		core_init();
		//these were in ctor, but failed because init hadn't happened; I think I can't init until here.

		//SF and Filter shareable
		sf = core::scoring::ScoreFunctionOP(new core::scoring::ScoreFunction());
		core::scoring::ScoreType ref(core::scoring::score_type_from_name("ref"));
		sf->set_weight(ref, 1);
		scorefilter = protocols::filters::FilterOP(
			new protocols::score_filters::ScoreTypeFilter(
			sf,
			ref,
			8675309 /*garbage number, with style*/));

		polyx = protocols::moves::MoverOP(new protocols::pose_creation::MakePolyXMover());

		//tiny, fast, yay!
		pose = create_twores_1ubq_pose(); //I think this fxn is global namespace
		//TR << "bottom setup" << std::endl;
	}

	void tearDown() {
		//TR << "teardown" << std::endl;
	}
	/// @details ok, I'm cheating, this is less unit-y than it could be; MoverB tests setters, MoverA tests complex ctor, (nothing tests parse_my_tag), and then the getters are tested here too.
	void test_FilterReportAsPoseExtraScoresMover() {
		//TR << "top main func" << std::endl;
		//note construct mover B(efore) with setters
		protocols::moves::FilterReportAsPoseExtraScoresMover moverB;
		moverB.set_filter(scorefilter);
		moverB.set_report_as(before);

		moverB.apply(pose); //should set an extra score
		core::Real s_before(sf->score(pose)); //score before we alter the pose

		polyx->apply(pose); //should reset the sequence to fiddle with scores

		//construct mover A(fter) with ctor
		protocols::moves::FilterReportAsPoseExtraScoresMover moverA(scorefilter, after);
		moverA.apply(pose); //should set an extra score
		core::Real s_after(sf->score(pose)); //score after we fiddle with it

		TS_ASSERT(core::pose::hasPoseExtraScore(pose, before));
		TS_ASSERT(core::pose::hasPoseExtraScore(pose, after));

		core::Real PES_stored(0);
		TS_ASSERT(core::pose::getPoseExtraScore(pose, before, PES_stored)); //returns TF on existence, PES_stored pass-by-ref
		TS_ASSERT_DELTA(s_before, PES_stored, 1e-6);

		TS_ASSERT(core::pose::getPoseExtraScore(pose, after, PES_stored)); //returns TF on existence, PES_stored pass-by-ref
		TS_ASSERT_DELTA(s_after, PES_stored, 1e-6);

		TR << "before: " << s_before << " after: " << s_after << std::endl;

		//test getters
		//ANDREW CHECK THESE LINES - is this the right way to check that a shallow-copy OP is the same?
		TS_ASSERT_EQUALS(scorefilter.get(), moverA.get_filter().get());
		TS_ASSERT_EQUALS(scorefilter.get(), moverB.get_filter().get());
		//END ANDREW CHECK THESE LINES
		TS_ASSERT_EQUALS(after, moverA.get_report_as());
		TS_ASSERT_EQUALS(before, moverB.get_report_as());

		TR << "test_OneFilterReportAsPoseExtraScoresMover completed!! " << std::endl;
	}//  void test_FilterReportAsPoseExtraScoresMover() {

};

