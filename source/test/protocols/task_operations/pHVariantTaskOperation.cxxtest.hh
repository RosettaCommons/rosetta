// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/task_operation/pHVariantTaskOperation.cxxtest.hh
/// @brief  Unit test for pHVariantTaskOperation: Allows protonation variants during packing for ionizable side chains
/// @author Rebecca Alford (rfalford12@gmail.com)
/// @author Labonte <JWLabonte@jhu.edu>

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <protocols/task_operations/pHVariantTaskOperation.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/chemical/ResidueType.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Utility, etc Headers
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR("pHVariantTaskOperation");

using namespace protocols::task_operations;
using namespace core::pack::task;

class pHVariantTaskOperationTest : public CxxTest::TestSuite {

public:

	void setUp(){
		core_init_with_additional_options( "-pH_mode" );
	}

	void tearDown(){

	}

	void test_pHVariantTaskOp() {

		// Sequence: AADEHKYAA
		core::pose::Pose pose;
		core::import_pose::pose_from_file(
			pose, "protocols/task_operations/pose_with_ionizable_side_chains.pdb", core::import_pose::PDB_file );

		TaskFactory pHVariant_factory;
		pHVariant_factory.push_back( utility::pointer::make_shared< pHVariantTaskOperation >() );
		PackerTaskOP task( pHVariant_factory.create_task_and_apply_taskoperations(pose) );

		// The task should only permit design at protonatable positions, so residues 3 through 7.
		TS_ASSERT( ! task->residue_task( 1 ).being_designed() );
		TS_ASSERT( ! task->residue_task( 2 ).being_designed() );
		TS_ASSERT( task->residue_task( 3 ).being_designed() );
		TS_ASSERT( task->residue_task( 4 ).being_designed() );
		TS_ASSERT( task->residue_task( 5 ).being_designed() );
		TS_ASSERT( task->residue_task( 6 ).being_designed() );
		TS_ASSERT( task->residue_task( 7 ).being_designed() );
		TS_ASSERT( ! task->residue_task( 8 ).being_designed() );
		TS_ASSERT( ! task->residue_task( 9 ).being_designed() );

		// Residues 3, 4, and 5 have two additional protonation states compared to residues 6 and 7.
		TS_ASSERT_EQUALS( task->residue_task( 3 ).allowed_residue_types().size(), 3 );
		utility::vector1< std::string > allowed_types3;
		for ( auto type : task->residue_task( 3 ).allowed_residue_types() ) {
			allowed_types3.push_back( type->name() );
		}
		TS_ASSERT( allowed_types3.has_value( "ASP" ) );
		TS_ASSERT( allowed_types3.has_value( "ASP_P1" ) );
		TS_ASSERT( allowed_types3.has_value( "ASP_P2" ) );

		TS_ASSERT_EQUALS( task->residue_task( 4 ).allowed_residue_types().size(), 3 );
		utility::vector1< std::string > allowed_types4;
		for ( auto type : task->residue_task( 4 ).allowed_residue_types() ) {
			allowed_types4.push_back( type->name() );
		}
		TS_ASSERT( allowed_types4.has_value( "GLU" ) );
		TS_ASSERT( allowed_types4.has_value( "GLU_P1" ) );
		TS_ASSERT( allowed_types4.has_value( "GLU_P2" ) );

		TS_ASSERT_EQUALS( task->residue_task( 5 ).allowed_residue_types().size(), 3 );
		utility::vector1< std::string > allowed_types5;
		for ( auto type : task->residue_task( 5 ).allowed_residue_types() ) {
			allowed_types5.push_back( type->name() );
		}
		TS_ASSERT( allowed_types5.has_value( "HIS" ) );
		TS_ASSERT( allowed_types5.has_value( "HIS_D" ) );
		TS_ASSERT( allowed_types5.has_value( "HIS_P" ) );

		TS_ASSERT_EQUALS( task->residue_task( 6 ).allowed_residue_types().size(), 2 );
		utility::vector1< std::string > allowed_types6;
		for ( auto type : task->residue_task( 6 ).allowed_residue_types() ) {
			allowed_types6.push_back( type->name() );
		}
		TS_ASSERT( allowed_types6.has_value( "LYS" ) );
		TS_ASSERT( allowed_types6.has_value( "LYS_D" ) );

		TS_ASSERT_EQUALS( task->residue_task( 7 ).allowed_residue_types().size(), 2 );
		utility::vector1< std::string > allowed_types7;
		for ( auto type : task->residue_task( 7 ).allowed_residue_types() ) {
			allowed_types7.push_back( type->name() );
		}
		TS_ASSERT( allowed_types7.has_value( "TYR" ) );
		TS_ASSERT( allowed_types7.has_value( "TYR_D" ) );
	}
};
