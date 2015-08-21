// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/RestrictNonSurfaceToRepackingOperation.cxxtest.hh
/// @brief  test for the task operation RestrictNonSurfaceToRepackingOperation
/// @author Ron Jacak ronj@unc.edu

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit header
#include <protocols/toolbox/task_operations/RestrictNonSurfaceToRepackingOperation.hh>

// project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>

// used for extra verification during initial test writing. leaving in just in case I have to test again
//#include <ObjexxFCL/format.hh>
//#include <core/scoring/ScoreFunction.hh>
//#include <core/scoring/ScoreFunctionFactory.hh>
//#include <core/scoring/Energies.hh>
//#include <core/scoring/TenANeighborGraph.hh>


// utility headers
#include <basic/Tracer.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <utility/vector1.hh>


// C++ headers


static basic::Tracer TR("test.protocols.toolbox.task_operations.RestrictNonSurfaceToRepackingOperationTests");


// --------------- Test Class --------------- //

class RestrictNonSurfaceToRepackingOperationTests : public CxxTest::TestSuite {

	core::pose::Pose pose;
	//core::scoring::ScoreFunctionOP scorefxn;

public:

	// --------------- Fixtures --------------- //

	//ctor sets up the pose once
	RestrictNonSurfaceToRepackingOperationTests(){
		core_init_with_additional_options( "-no_optH -mute core.io core.init core.scoring core.mm core.pack.task" );
		core::import_pose::pose_from_pdb( pose, "protocols/util/test_in.pdb" );
	}
	virtual ~RestrictNonSurfaceToRepackingOperationTests() {}
	static RestrictNonSurfaceToRepackingOperationTests* createSuite() {
		return new RestrictNonSurfaceToRepackingOperationTests();
	}
	static void destroySuite( RestrictNonSurfaceToRepackingOperationTests *suite ) { delete suite; }

	void setUp() {}
	void tearDown() {}

	// ------------- Helper Functions ------------- //


	// --------------- Test Cases --------------- //

	/// @brief this test applies the RestrictNonSurfaceToRepackingOperation task operation to a PackerTask and verifies
	/// that all of the residues that have more neighbors than the cutoff have been set to repacking only and that all
	/// the positions with neighbors fewer than the cutoff are designable.

	void test_RestrictNonSurfaceToRepackingOperation() {

		using namespace core;
		using namespace core::pack::task;
		using protocols::toolbox::task_operations::RestrictNonSurfaceToRepackingOperation;
		using core::pack::task::operation::TaskOperationCOP;
		//using namespace ObjexxFCL::format;

		TR << "Running test_RestrictNonSurfaceToRepackingOperation..." << std::endl;

		core::pack::task::TaskFactory tf;
		tf.push_back( TaskOperationCOP( new RestrictNonSurfaceToRepackingOperation() ) ); // default to cutoff of 16

		// the following call should work even if the pose hasn't been scored
		core::pack::task::PackerTaskOP task = tf.create_task_and_apply_taskoperations( pose );
		//std::cout << *task << std::endl; // have to dereference the OP.  printing this only for debugging

		// verify the nb count calculation using the tenA neighbor graph
		//scorefxn = core::scoring::get_score_function();
		//(*scorefxn)(pose);

		//for ( Size ii=1; ii <= pose.total_residue(); ++ii ) {
		// std::cout << "resid: " << I(3,ii)
		//   << ", nbs: " << pose.energies().tenA_neighbor_graph().get_node( ii )->num_neighbors_counting_self() << std::endl;
		//}

		// The first 5 residues of the util/test_in.pdb file have the following neighbor counts. By default, the task operation uses
		// 16 as the cutoff.  If a residue has more than 16 nbs, it should be limited to packing only. design should be enabled for
		// positions with 16 or fewer neighbors only.
		//
		// resid:   1, nbs: 17
		// resid:   2, nbs: 11
		// resid:   3, nbs: 8
		// resid:   4, nbs: 16
		// resid:   5, nbs: 18
		//
		// So a "correct" packertask would look as follows:
		// resid   pack?   design? allowed_aas
		// 1       1       0       ASP:NtermProteinFull,
		// 2       1       1       ALA,CYS,ASP,GLU,PHE,GLY,HIS,HIS_D,ILE,LYS,LEU,MET,ASN,PRO,GLN,ARG,SER,THR,VAL,TRP,TYR,
		// 3       1       1       ALA,CYS,ASP,GLU,PHE,GLY,HIS,HIS_D,ILE,LYS,LEU,MET,ASN,PRO,GLN,ARG,SER,THR,VAL,TRP,TYR,
		// 4       1       1       ALA,CYS,ASP,GLU,PHE,GLY,HIS,HIS_D,ILE,LYS,LEU,MET,ASN,PRO,GLN,ARG,SER,THR,VAL,TRP,TYR,
		// 5       1       0       ILE,

		TS_ASSERT_EQUALS( task->residue_task(1).being_packed(), true );
		TS_ASSERT_EQUALS( task->residue_task(1).being_designed(), false );

		TS_ASSERT_EQUALS( task->residue_task(2).being_packed(), true );
		TS_ASSERT_EQUALS( task->residue_task(2).being_designed(), true );

		TS_ASSERT_EQUALS( task->residue_task(3).being_packed(), true );
		TS_ASSERT_EQUALS( task->residue_task(3).being_designed(), true );

		TS_ASSERT_EQUALS( task->residue_task(4).being_packed(), true );
		TS_ASSERT_EQUALS( task->residue_task(4).being_designed(), true );

		TS_ASSERT_EQUALS( task->residue_task(5).being_packed(), true );
		TS_ASSERT_EQUALS( task->residue_task(5).being_designed(), false );

	}

};

