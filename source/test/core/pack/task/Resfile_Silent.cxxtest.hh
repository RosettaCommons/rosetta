// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/Resfile_Silent.cxxtest.hh
/// @brief  test suite for resfile with silent file
/// @author Nobuasu Koga

// Test headers
#include <cxxtest/TestSuite.h>

#include "platform/types.hh"

#include <test/core/init_util.hh>

// Package Headers
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>

#include <string>

//Auto Headers
#include <core/chemical/ResidueType.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/silent/EnergyNames.fwd.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.io.Resfile_Silent.cxxtest");

using namespace core;

class Resfile_Silent : public CxxTest::TestSuite{

	pack::task::operation::ReadResfileOP rrop;
	chemical::ResidueTypeSetCOP residue_set;

public:

	typedef core::import_pose::pose_stream::SilentFilePoseInputStream SilentFilePoseInputStream;
	typedef core::import_pose::pose_stream::SilentFilePoseInputStreamOP SilentFilePoseInputStreamOP;
	typedef core::pack::task::ResidueLevelTask ResidueLevelTask;
	typedef core::pack::task::PackerTask PackerTask;
	typedef core::pack::task::PackerTaskOP PackerTaskOP;
	typedef core::pack::task::TaskFactory TaskFactory;
	typedef core::pack::task::operation::ReadResfile ReadResfile;
	typedef core::pose::Pose Pose;
	typedef core::pose::PoseOP PoseOP;

public:

	Resfile_Silent() {};

	// Shared initialization goes here.
	void setUp() {
		core_init_with_additional_options( "-no_optH" );
		residue_set = chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD );
		rrop = new ReadResfile( "core/pack/task/test_in.resfile" );
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	// ------------------------------------------ //
 	/// @brief test how RESFILE IO

	bool compare_packertasks( PackerTask const & tA, PackerTask const & tB )
	{
		TR << tA << std::endl;
		TS_ASSERT_EQUALS( tA.total_residue(), tA.total_residue() );
		for ( int i=1, it_end = tA.total_residue(); i <= it_end; ++i){
			TS_ASSERT_EQUALS( tA.pack_residue( i ), tB.pack_residue( i ) );
			TS_ASSERT_EQUALS( tA.design_residue( i ), tB.design_residue( i ) );
			ResidueLevelTask::ResidueTypeCOPListConstIter itA( tA.residue_task( i ).allowed_residue_types_begin() );
			ResidueLevelTask::ResidueTypeCOPListConstIter itB( tB.residue_task( i ).allowed_residue_types_begin() );
			while( itA != tA.residue_task( i ).allowed_residue_types_end() )
				{
					TS_ASSERT_EQUALS( (*itA)->name(), (*itB)->name() );
					++itA;
					++itB;
				}
		}
		TR << tB << std::endl;
		return true;
	}

	void test_silent_resfile_io()
	{
		// set packer task from pdb
		Pose pdb_pose;
		std::string pdb_file_name( "core/pack/task/test_in.pdb" );
		core::import_pose::pose_from_pdb( pdb_pose, pdb_file_name );
		PackerTaskOP task_pdb( TaskFactory::create_packer_task( pdb_pose ) );
		rrop->apply( pdb_pose, *task_pdb );

		// set packer task from silent file
		PoseOP silent_pose;
		const std::string silent_file_name( "core/pack/task/test_in.silent" );
		SilentFilePoseInputStreamOP silent_input = new SilentFilePoseInputStream( silent_file_name );
		silent_pose =	silent_input->get_all_poses( *residue_set )[ 1 ];
		PackerTaskOP task_silent( TaskFactory::create_packer_task( *silent_pose ) );
		rrop->apply( *silent_pose, *task_silent );

		if( compare_packertasks( *task_pdb, *task_silent ) ){
			TR << "END_OF_TEST" << std::endl;
		}
	}


};

