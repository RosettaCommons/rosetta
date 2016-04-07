// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/select/residue_selector/TaskSelector.cxxtest.hh
/// @brief test suite for core::select::residue_selector::TaskSelector
/// @author Tom Linsky (tlinsky@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/core/select/residue_selector/DummySelectors.hh>

// Package headers
#include <protocols/residue_selectors/TaskSelector.hh>
#include <protocols/toolbox/task_operations/DesignAroundOperation.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/dssp/Dssp.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>

// C++ headers
#include <boost/assign.hpp>
#include <string>

using namespace core::select::residue_selector;

static THREAD_LOCAL basic::Tracer TR( "test.core.select.TaskSelectorTests" );

class TaskSelectorTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void test_TaskSelector_parse_my_tag() {
		using namespace protocols::residue_selectors;
		std::stringstream ss;
		ss << "<Task name=\"allloops\" />";
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );
		basic::datacache::DataMap dm;

		TaskSelectorOP rs( new TaskSelector );
		TS_ASSERT_THROWS( rs->parse_my_tag( tag, dm ), utility::excn::EXCN_RosettaScriptsOption );

		std::stringstream ssgood;
		ssgood << "<Task name=\"good_tos\" task_operations=\"test_op,test_op2\" packable=\"1\" designable=\"0\" />";
		tag->read( ssgood );
		TS_ASSERT_THROWS( rs->parse_my_tag( tag, dm ), utility::excn::EXCN_RosettaScriptsOption );

		protocols::toolbox::task_operations::DesignAroundOperationOP des_around( new protocols::toolbox::task_operations::DesignAroundOperation );
		des_around->include_residue( 2 );
		dm.add( "task_operations", "test_op", des_around );
		dm.add( "task_operations", "test_op2", des_around );
		TS_ASSERT_THROWS_NOTHING( rs->parse_my_tag( tag, dm ) );

		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		core::select::residue_selector::ResidueSubset subset = rs->apply( trpcage );
		trpcage.dump_pdb( "trpcage.pdb" );
		std::set< core::Size > const residues = boost::assign::list_of (1)(2)(3)(4)(5)(6)(19);
		for ( core::Size resid=1; resid<=trpcage.total_residue(); ++resid ) {
			bool const found = ( residues.find( resid ) != residues.end() );
			TS_ASSERT( found == subset[resid] );
			TR << resid << " " << subset[resid] << std::endl;
		}

		// select frozen
		rs->set_select_designable( false );
		rs->set_select_packable( false );
		rs->set_select_fixed( true );
		subset = rs->apply( trpcage );
		for ( core::Size resid=1; resid<=trpcage.total_residue(); ++resid ) {
			bool const not_found = ( residues.find( resid ) == residues.end() );
			TS_ASSERT( not_found == subset[resid] );
			TR << resid << " " << subset[resid] << std::endl;
		}

	}

};
