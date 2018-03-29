// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// @author Brian Koepnick (koepnick@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/UTracer.hh>

// Unit headers
#include <protocols/parser/BluePrint.hh>
#include <protocols/task_operations/DesignAroundOperation.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/PackerTask_.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/annotated_sequence.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <string>
#include <map>

static basic::Tracer TR("protocols.task_operations.DesignAroundOperationTests.cxxtest");

using namespace protocols::task_operations;
using namespace core;
using namespace core::pack::task;

class DesignAroundOperationTests : public CxxTest::TestSuite {

public:
	core::pose::PoseOP test_dimer_pose_ = nullptr;

	void setUp(){
		protocols_init();
		test_dimer_pose_ = create_2res_1ten_2res_trp_cage_poseop(); //dimer structure
	}

	void test_it_restricts_outside_pack_radius() {

		PackerTaskOP initial_task( new PackerTask_( *test_dimer_pose_ ) );
		// repacking residues should be pose.size()
		auto rr = initial_task->repacking_residues();
		Size n_repacking_init = std::count( rr.begin(), rr.end(), true );
		TS_ASSERT( n_repacking_init == test_dimer_pose_->size() );

		TaskFactory factory;
		DesignAroundOperationOP dao( new DesignAroundOperation );
		dao->include_residue( 1 );
		factory.push_back( dao );
		core::pack::task::PackerTaskOP task( factory.create_task_and_apply_taskoperations( *test_dimer_pose_ ) );

		rr = task->repacking_residues();
		Size n_repacking_final = std::count( rr.begin(), rr.end(), true );

		TS_ASSERT( n_repacking_final < n_repacking_init );
	}

	void test_it_restricts_outside_design_radius() {
		PackerTaskOP initial_task( new PackerTask_( *test_dimer_pose_ ) );
		// designing residues should be pose.size()
		auto dr = initial_task->designing_residues();
		Size n_designing_init = std::count( dr.begin(), dr.end(), true );
		TS_ASSERT( n_designing_init == test_dimer_pose_->size() );

		TaskFactory factory;
		DesignAroundOperationOP dao( new DesignAroundOperation );
		dao->include_residue( 1 );
		factory.push_back( dao );
		core::pack::task::PackerTaskOP task( factory.create_task_and_apply_taskoperations( *test_dimer_pose_ ) );

		dr = task->designing_residues();
		Size n_designing_final = std::count( dr.begin(), dr.end(), true );

		TS_ASSERT( n_designing_final < n_designing_init );
	}

	// We just need to set the design radius here, because we will
	// automatically reset the repack radius to match.
	void test_it_doesnt_restrict_for_big_enough_radius() {
		PackerTaskOP initial_task( new PackerTask_( *test_dimer_pose_ ) );
		// designing residues should be pose.size()
		auto dr = initial_task->designing_residues();
		auto rr = initial_task->repacking_residues();
		Size n_designing_init = std::count( dr.begin(), dr.end(), true );
		Size n_repacking_init = std::count( rr.begin(), rr.end(), true );
		TS_ASSERT( n_designing_init == test_dimer_pose_->size() );
		TS_ASSERT( n_repacking_init == test_dimer_pose_->size() );

		TaskFactory factory;
		DesignAroundOperationOP dao( new DesignAroundOperation );
		dao->include_residue( 1 );
		dao->design_shell( 1000000 );
		factory.push_back( dao );
		core::pack::task::PackerTaskOP task( factory.create_task_and_apply_taskoperations( *test_dimer_pose_ ) );

		dr = task->designing_residues();
		rr = task->repacking_residues();
		Size n_designing_final = std::count( dr.begin(), dr.end(), true );
		Size n_repacking_final = std::count( rr.begin(), rr.end(), true );

		TS_ASSERT( n_designing_final == test_dimer_pose_->size() );
		TS_ASSERT( n_repacking_final == test_dimer_pose_->size() );
	}

};
