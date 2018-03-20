// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/moves/MoversTest.cxxtest.hh
/// @brief  tests for container Movers classes.
/// @author Sergey Lyskov

#ifndef INCLUDED_UMoverTest_HH
#define INCLUDED_UMoverTest_HH


// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/types.hh>

#include <protocols/moves/Mover.hh>

#include <basic/Tracer.hh>
#include <test/UTracer.hh>


#define TEST_MOVER(mover, fileIn, fileOut) one_mover_test(__FILE__, __LINE__, protocols::moves::MoverOP( new mover() ), fileIn, fileOut);
#define TEST_MOVER_OP(mover_op, fileIn, fileOut) one_mover_test(__FILE__, __LINE__, mover, fileIn, fileOut);

namespace test {


// using declarations
using namespace core;
using namespace core::pose;
using namespace protocols::moves;


//using basic::T;
//static basic::Tracer TR_I("cxxtest");

///////////////////////////////////////////////////////////////////////////
/// @name UMoverTest
/// @brief: class for Movers unified testing
/// @author Sergey Lyskov
///////////////////////////////////////////////////////////////////////////
class UMoverTest /*public CxxTest::TestSuite*/ {

public:
	chemical::ResidueTypeSetCAP residue_set;

	UMoverTest() {}

	void setUp() {
		core_init_with_additional_options( "-no_optH" );

		residue_set = chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD );
	}


	/// @brief: Helper function, that execute test on a given Mover object
	///
	void one_mover_test(const char * /*file*/, unsigned /*line*/, protocols::moves::MoverOP mover,
						const char *fileIn, const char *fileOut, const char *nativeFileIn=0 ) {
		//TR_I << "MoversTest: Testing " << mover->type().c_str() << "...\n";

		core::init::init_random_generators(1000, "mt19937");

		const std::string original_file_name( fileIn );

		//std::cout << " Testing: " << mover->type() << "..." << std::endl;

		PoseOP pose( new Pose );
		PoseOP n_pose( new Pose ); // native pose
		core::import_pose::pose_from_file( *pose, original_file_name , core::import_pose::PDB_file);
		if ( nativeFileIn ) {
			const std::string native_file_name( nativeFileIn );
			core::import_pose::pose_from_file( *n_pose, native_file_name , core::import_pose::PDB_file);
		} else {
			n_pose = pose;
		}

		mover->set_input_pose( pose );
		mover->set_native_pose( n_pose );
		mover->test_move( *pose );

		// We use a UTracer here because it allows for fuzzy comparison of real numbers in the pose (e.g. energies)
		UTracer UT( fileOut );
		UT.abs_tolerance( 1e-6 );
		pose->dump_pdb( UT );

	}

	void tearDown() {
	}

};

} // namespace test

#endif
