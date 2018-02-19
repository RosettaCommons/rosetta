// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/simple_filters/LongestContinuousApolarSegmentFilterTests.cxxtest.hh
/// @brief  Unit tests for the LongestContinuousApolarSegmentFilter.
/// @author Yang Hsia (yhsia@uw.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers


// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>

// Protocol Headers
#include <protocols/simple_filters/LongestContinuousApolarSegmentFilter.hh>

// Basic Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("LongestContinuousApolarSegmentFilterTests");


class LongestContinuousApolarSegmentFilterTests : public CxxTest::TestSuite {

public:

	core::pose::Pose pose1_; //1-chain pose
	core::pose::Pose pose2_; //2-chain pose

public:

	void setUp(){
		core_init(); //                               1234567890123456789012345678901234567890
		//core::pose::make_pose_from_sequence( pose1_, "QGRDSNKVFWDRHTLIDEVWLLIFQDRE", "fa_standard" ); //7 res polar, 3 apolar, 4 polar, 2 apolar, 2 polar, 6 apolar, 4 polar.
		core::pose::make_pose_from_sequence( pose1_, "AAAAAAASSSAAAASSAASSSSSSAAAA", "fa_standard" ); //7 res apolar, 3 polar, 4 apolar, 2 polar, 2 apolar, 6 polar, 4 apolar.
		//                                            1234567890123456789012345678901234567890
		//core::pose::make_pose_from_sequence( pose2_, "QDREVKERDSSKKKSRDVKERDSSKKKSRKERDSSKKKSR", "fa_standard" ); //4 polar, 1 apolar, 12 polar, 1 apolar, 22 polar.
		core::pose::make_pose_from_sequence( pose2_, "AAAASAAAAAAAAAAAASAAAAAAAAAAAAAAAAAAAAAA", "fa_standard" ); //4 apolar, 1 polar, 12 apolar, 1 polar, 22 apolar.
		pose2_.append_pose_by_jump( pose1_, 1 );
	}

	void tearDown(){}


	void test_singlechain_pose_withtermini(){
		TR << "Starting LongestContinuousApolarSegmentFilterTests::test_singlechain_pose_withtermini()." << std::endl;
		protocols::simple_filters::LongestContinuousApolarSegmentFilter myfilter;
		myfilter.set_exclude_chain_termini(false);

		bool const passfail1( myfilter.apply( pose1_) );
		core::Real const val1( myfilter.report_sm( pose1_ ) );
		myfilter.set_cutoff(7);
		bool const passfail2( myfilter.apply( pose1_ ) );
		core::Real const val2( myfilter.report_sm( pose1_ ) );

		std::stringstream desc;
		myfilter.report(desc, pose1_);

		TS_ASSERT( passfail1 == false );
		TS_ASSERT( passfail2 == true );
		TS_ASSERT_DELTA( val1, 7.0, 1e-6 );
		TS_ASSERT_DELTA( val2, 7.0, 1e-6 );
		TS_ASSERT_EQUALS( desc.str(), "The LongestContinuousApolarSegmentFilter reports that the longest stretch of apolar residues in the pose, including chain terminal segments, is 7 residues long.  It runs from ALA1 through ALA7.\n" );

		TR << "Finished LongestContinuousApolarSegmentFilterTests::test_singlechain_pose_withtermini()." << std::endl;
	}

	void test_singlechain_pose_withouttermini(){
		TR << "Starting LongestContinuousApolarSegmentFilterTests::test_singlechain_pose_withouttermini()." << std::endl;
		protocols::simple_filters::LongestContinuousApolarSegmentFilter myfilter;
		myfilter.set_exclude_chain_termini(true);

		bool const passfail1( myfilter.apply( pose1_) );
		core::Real const val1( myfilter.report_sm( pose1_ ) );
		myfilter.set_cutoff(3);
		bool const passfail2( myfilter.apply( pose1_ ) );
		core::Real const val2( myfilter.report_sm( pose1_ ) );

		std::stringstream desc;
		myfilter.report(desc, pose1_);

		TS_ASSERT( passfail1 == true );
		TS_ASSERT( passfail2 == false );
		TS_ASSERT_DELTA( val1, 4.0, 1e-6 );
		TS_ASSERT_DELTA( val2, 4.0, 1e-6 );
		TS_ASSERT_EQUALS( desc.str(), "The LongestContinuousApolarSegmentFilter reports that the longest stretch of apolar residues in the pose, excluding chain terminal segments, is 4 residues long.  It runs from ALA11 through ALA14.\n" );

		TR << "Finished LongestContinuousApolarSegmentFilterTests::test_singlechain_pose_withouttermini()." << std::endl;
	}

	void test_singlechain_pose_withselection(){
		TR << "Starting LongestContinuousApolarSegmentFilterTests::test_singlechain_pose_withselection()." << std::endl;
		protocols::simple_filters::LongestContinuousApolarSegmentFilter myfilter;
		myfilter.set_exclude_chain_termini(true);
		myfilter.set_filter_out_high(false);
		core::select::residue_selector::ResidueIndexSelectorOP selector( new core::select::residue_selector::ResidueIndexSelector );
		selector->set_index( "4,15-28" );
		myfilter.set_residue_selector(selector);

		bool const passfail1( myfilter.apply( pose1_) );
		core::Real const val1( myfilter.report_sm( pose1_ ) );
		myfilter.set_cutoff(2);
		bool const passfail2( myfilter.apply( pose1_ ) );
		core::Real const val2( myfilter.report_sm( pose1_ ) );

		std::stringstream desc;
		myfilter.report(desc, pose1_);

		TS_ASSERT( passfail1 == false );
		TS_ASSERT( passfail2 == true );
		TS_ASSERT_DELTA( val1, 2.0, 1e-6 );
		TS_ASSERT_DELTA( val2, 2.0, 1e-6 );
		TS_ASSERT_EQUALS( desc.str(), "The LongestContinuousApolarSegmentFilter reports that the longest stretch of apolar residues in the pose, excluding chain terminal segments, is 2 residues long.  It runs from ALA17 through ALA18.\n" );

		TR << "Finished LongestContinuousApolarSegmentFilterTests::test_singlechain_pose_withselection()." << std::endl;
	}

	void test_multichain_pose_withtermini(){
		TR << "Starting LongestContinuousApolarSegmentFilterTests::test_multichain_pose_withtermini()." << std::endl;
		protocols::simple_filters::LongestContinuousApolarSegmentFilter myfilter;
		myfilter.set_exclude_chain_termini(false);

		bool const passfail1( myfilter.apply( pose2_) );
		core::Real const val1( myfilter.report_sm( pose2_ ) );
		myfilter.set_cutoff(22);
		bool const passfail2( myfilter.apply( pose2_ ) );
		core::Real const val2( myfilter.report_sm( pose2_ ) );

		std::stringstream desc;
		myfilter.report(desc, pose2_);

		TS_ASSERT( passfail1 == false );
		TS_ASSERT( passfail2 == true );
		TS_ASSERT_DELTA( val1, 22.0, 1e-6 );
		TS_ASSERT_DELTA( val2, 22.0, 1e-6 );
		TS_ASSERT_EQUALS( desc.str(), "The LongestContinuousApolarSegmentFilter reports that the longest stretch of apolar residues in the pose, including chain terminal segments, is 22 residues long.  It runs from ALA19 through ALA40.\n" );

		TR << "Finished LongestContinuousApolarSegmentFilterTests::test_multichain_pose_withtermini()." << std::endl;
	}

	void test_multichain_pose_withouttermini(){
		TR << "Starting LongestContinuousApolarSegmentFilterTests::test_multichain_pose_withouttermini()." << std::endl;
		protocols::simple_filters::LongestContinuousApolarSegmentFilter myfilter;
		myfilter.set_exclude_chain_termini(true);

		bool const passfail1( myfilter.apply( pose2_) );
		core::Real const val1( myfilter.report_sm( pose2_ ) );
		myfilter.set_cutoff(12);
		bool const passfail2( myfilter.apply( pose2_ ) );
		core::Real const val2( myfilter.report_sm( pose2_ ) );

		std::stringstream desc;
		myfilter.report(desc, pose2_);

		TS_ASSERT( passfail1 == false );
		TS_ASSERT( passfail2 == true );
		TS_ASSERT_DELTA( val1, 12.0, 1e-6 );
		TS_ASSERT_DELTA( val2, 12.0, 1e-6 );
		TS_ASSERT_EQUALS( desc.str(), "The LongestContinuousApolarSegmentFilter reports that the longest stretch of apolar residues in the pose, excluding chain terminal segments, is 12 residues long.  It runs from ALA6 through ALA17.\n" );
		TR << "Finished LongestContinuousApolarSegmentFilterTests::test_multichain_pose_withouttermini()." << std::endl;
	}

	void test_multichain_pose_withselection(){
		TR << "Starting LongestContinuousApolarSegmentFilterTests::test_multichain_pose_withselection()." << std::endl;
		protocols::simple_filters::LongestContinuousApolarSegmentFilter myfilter;
		myfilter.set_exclude_chain_termini(true);
		myfilter.set_filter_out_high(false);
		core::select::residue_selector::ResidueIndexSelectorOP selector( new core::select::residue_selector::ResidueIndexSelector );
		selector->set_index( "41-68" ); //Just the second chain
		myfilter.set_residue_selector(selector);

		bool const passfail1( myfilter.apply( pose2_) );
		core::Real const val1( myfilter.report_sm( pose2_ ) );
		myfilter.set_cutoff(4);
		bool const passfail2( myfilter.apply( pose2_ ) );
		core::Real const val2( myfilter.report_sm( pose2_ ) );

		std::stringstream desc;
		myfilter.report(desc, pose2_);

		TS_ASSERT( passfail1 == false );
		TS_ASSERT( passfail2 == true );
		TS_ASSERT_DELTA( val1, 4.0, 1e-6 );
		TS_ASSERT_DELTA( val2, 4.0, 1e-6 );
		TS_ASSERT_EQUALS( desc.str(), "The LongestContinuousApolarSegmentFilter reports that the longest stretch of apolar residues in the pose, excluding chain terminal segments, is 4 residues long.  It runs from ALA51 through ALA54.\n" );

		TR << "Finished LongestContinuousApolarSegmentFilterTests::test_multichain_pose_withselection()." << std::endl;
	}



};



