// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/pack/task/residue_selector/OperateOnResidueSubset.cxxtest.hh
/// @brief  test suite for core::pack::task::residue_selector::OperateOnResidueSubset
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/core/pack/task/residue_selector/DummySelectors.hh>

// Package headers
#include <core/pack/task/operation/OperateOnResidueSubset.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/pack/task/operation/TaskOperationFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

// Project headers
#include <core/pose/Pose.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>

// C++ headers
#include <string>

using namespace core::pack::task;
using namespace core::pack::task::operation;
using namespace core::pack::task::residue_selector;

class OperateOnResidueSubsetTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	// @brief make sure that when we register a residue selector, we can later get it back
	void test_disable_odd_residues_init_by_ctor() {
		ResidueSelectorOP odd_rs( new OddResidueSelector );
		ResLvlTaskOperationOP prev_repacking( new PreventRepackingRLT );
		OperateOnResidueSubset op_on_subset( prev_repacking, odd_rs );

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		PackerTaskOP task = TaskFactory::create_packer_task( trpcage );

		op_on_subset.apply( trpcage, *task );
		for ( core::Size ii = 1; ii <= trpcage.total_residue(); ++ii ) {
			TS_ASSERT_EQUALS( task->residue_task(ii).being_packed(), ii % 2 == 0 );
		}
	}

	void test_disable_odd_residues_init_w_functions() {
		ResidueSelectorOP odd_rs( new OddResidueSelector );
		ResLvlTaskOperationOP prev_repacking( new PreventRepackingRLT );
		OperateOnResidueSubset op_on_subset;
		op_on_subset.op( prev_repacking );
		op_on_subset.selector( odd_rs );

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		PackerTaskOP task = TaskFactory::create_packer_task( trpcage );

		op_on_subset.apply( trpcage, *task );
		for ( core::Size ii = 1; ii <= trpcage.total_residue(); ++ii ) {
			TS_ASSERT_EQUALS( task->residue_task(ii).being_packed(), ii % 2 == 0 );
		}
	}

	/// @brief Test parse_tag initializing the selector from an option
	void test_oors_parse_tag_1() {
		std::string tag_string = "<OperateOnResidueSubset name=disable_odd selector=odd> <PreventRepackingRLT/> </OperateOnResidueSubset>";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;
		ResidueSelectorOP odd_rs( new OddResidueSelector );

		dm.add( "ResidueSelector", "odd", odd_rs );

		OperateOnResidueSubset op_on_subset;
		try {
			op_on_subset.parse_tag( tag, dm );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Exception!" << e.msg() << std::endl;
			TS_ASSERT( false ); // this parsing should succeed
		}

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		PackerTaskOP task = TaskFactory::create_packer_task( trpcage );

		op_on_subset.apply( trpcage, *task );
		for ( core::Size ii = 1; ii <= trpcage.total_residue(); ++ii ) {
			TS_ASSERT_EQUALS( task->residue_task(ii).being_packed(), ii % 2 == 0 );
		}
	}

	void test_oors_parse_tag_2() {
		std::string tag_string = "<OperateOnResidueSubset name=disable_1to10> <Index resnums=\"1-10\"/> <PreventRepackingRLT/> </OperateOnResidueSubset>";
		std::stringstream ss( tag_string );
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		OperateOnResidueSubset op_on_subset;
		try {
			op_on_subset.parse_tag( tag, dm );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Exception!" << e.msg() << std::endl;
			TS_ASSERT( false ); // this parsing should succeed
		}

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		PackerTaskOP task = TaskFactory::create_packer_task( trpcage );

		op_on_subset.apply( trpcage, *task );
		for ( core::Size ii = 1; ii <= trpcage.total_residue(); ++ii ) {
			TS_ASSERT_EQUALS( task->residue_task(ii).being_packed(), ii > 10 );
		}
	}



};
