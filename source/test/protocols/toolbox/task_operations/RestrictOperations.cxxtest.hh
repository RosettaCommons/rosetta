// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/toolbox/task_operations/RestrictToInterfaceOperation.cxxtest.hh
/// @brief  test for RestrictToOperations
/// @author Steven Lewis
/// @author Jared Adolf-Bryfogle

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/UTracer.hh>

// Unit header
#include <protocols/toolbox/task_operations/RestrictToInterfaceOperation.hh>
#include <protocols/toolbox/task_operations/RestrictToNeighborhoodOperation.hh>
#include <protocols/toolbox/task_operations/RestrictByCalculatorsOperation.hh>
#include <protocols/toolbox/task_operations/RestrictToInterfaceVectorOperation.hh>
#include <protocols/toolbox/task_operations/RestrictInterGroupVectorOperation.hh>
#include <protocols/toolbox/task_operations/RestrictToMoveMapChiOperation.hh>

// project headers
#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <core/pose/metrics/CalculatorFactory.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/kinematics/MoveMap.hh>

#include <protocols/toolbox/pose_metric_calculators/NeighborhoodByDistanceCalculator.hh>
#include <core/pose/metrics/simple_calculators/InterfaceNeighborDefinitionCalculator.hh>

#include <basic/Tracer.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <utility/vector1.hh>


// C++ headers
#include <set>

static basic::Tracer TR("protocols.toolbox.task_operations.RestrictOperations.cxxtest");

// --------------- Test Class --------------- //

class RestrictOperationsTests : public CxxTest::TestSuite {

	core::pose::Pose pose;

public:

	// --------------- Fixtures --------------- //

	//ctor sets up the pose once
	RestrictOperationsTests(){
		core_init();
		//reuse for comparison with Interface class
		core::import_pose::centroid_pose_from_pdb( pose, "core/conformation/dock_in.pdb" ); //JAB - centroid for speed?
	}

	virtual ~RestrictOperationsTests() {}

	static RestrictOperationsTests* createSuite() {
		return new RestrictOperationsTests();
	}

	static void destroySuite( RestrictOperationsTests *suite ) {
		delete suite;
	}

	void setUp() {}

	void tearDown() {}

	// ------------- Helper Functions ------------- //


	// --------------- Test Cases --------------- //

	/// @brief
	void test_RestrictToInterfaceOperation() {

		/* if you are debugging this test, I suggest passing the flag -print_pymol_selection for ease of use.
		The flag -pose_metrics::interface_cutoff 8 will cause the results to match core/conformation/Interface.cxxtest
		although the flag actually defaults to 10 (thus a larger interface)
		*/

		//set up test
		using namespace core::pack::task;
		using protocols::toolbox::task_operations::RestrictToInterfaceOperation;
		using core::pack::task::operation::TaskOperationCOP;
		TaskFactory RTIO_factory;
		RTIO_factory.push_back( TaskOperationCOP( new RestrictToInterfaceOperation() ) ); //defaults to interface between chains 1 and 2

		//run
		test::UTracer UT_RTIO("protocols/toolbox/task_operations/RestrictToInterfaceOperation.u");
		//this call returns PackerTaskOP; we are dumping the ptask to utracer
		UT_RTIO << *(RTIO_factory.create_task_and_apply_taskoperations( pose )) << std::endl;

	}

	void test_RestrictToNeighborhoodOperation() {

		std::set< core::Size > crset;
		crset.insert(77); crset.insert(215), crset.insert(45); //surface, interface, buried residues to test

		using namespace core::pack::task;
		using protocols::toolbox::task_operations::RestrictToNeighborhoodOperation;
		using core::pack::task::operation::TaskOperationCOP;
		TaskFactory RTNO_factory;
		RTNO_factory.push_back( TaskOperationCOP( new RestrictToNeighborhoodOperation( crset ) ) );

		test::UTracer UT_RTNO("protocols/toolbox/task_operations/RestrictToNeighborhoodOperation.u");
		UT_RTNO << *(RTNO_factory.create_task_and_apply_taskoperations( pose )) << std::endl;
	}

	void test_RestrictByCalculatorsOperation() {
		using core::pose::metrics::PoseMetricCalculatorOP;

		//first we set up the calculators that the Operation will use
		std::string const interface_calc("interface"), neighborhood_calc("neighborhood");
		std::set< core::Size > crset_RBC;
		crset_RBC.insert(127); crset_RBC.insert(170), crset_RBC.insert(46);
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( interface_calc, PoseMetricCalculatorOP( new core::pose::metrics::simple_calculators::InterfaceNeighborDefinitionCalculator( core::Size(1), core::Size(2) ) ) );
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( neighborhood_calc, PoseMetricCalculatorOP( new protocols::toolbox::pose_metric_calculators::NeighborhoodByDistanceCalculator( crset_RBC ) ) );

		//this is the constructor parameter for the calculator - pairs of calculators and calculations to perform
		utility::vector1< std::pair< std::string, std::string> > calcs_and_calcns;
		calcs_and_calcns.push_back(std::make_pair(interface_calc, "interface_residues"));
		calcs_and_calcns.push_back(std::make_pair(neighborhood_calc, "neighbors"));

		using protocols::toolbox::task_operations::RestrictByCalculatorsOperation;
		using core::pack::task::operation::TaskOperationCOP;
		core::pack::task::TaskFactory RBC_factory;
		RBC_factory.push_back( TaskOperationCOP( new RestrictByCalculatorsOperation( calcs_and_calcns ) ) );

		test::UTracer UT_RBC("protocols/toolbox/task_operations/RestrictByCalculatorsOperation.u");
		UT_RBC << *(RBC_factory.create_task_and_apply_taskoperations( pose )) << std::endl;

	}//end test_RestrictByCalculatorsOperation

	void test_RestrictToInterfaceVectorOperation(){

		//set up test
		using namespace core::pack::task;
		using protocols::toolbox::task_operations::RestrictToInterfaceVectorOperation;
		using core::pack::task::operation::TaskOperationCOP;
		TaskFactory RTIVO_factory;
		//these are the default values but hard code anyway, test for chain #s
		RTIVO_factory.push_back( TaskOperationCOP( new RestrictToInterfaceVectorOperation(1,2,10,5.5,75,9.0) ) );
		//run
		test::UTracer UT_RTIVO("protocols/toolbox/task_operations/RestrictToInterfaceVectorOperation.u");
		//this call returns PackerTaskOP; we are dumping the ptask to utracer
		//std::cout << *(RTIVO_factory.create_task_and_apply_taskoperations( pose )) << std::endl;
		UT_RTIVO << *(RTIVO_factory.create_task_and_apply_taskoperations( pose )) << std::endl;

		//now test again for jump deffinition
		test::UTracer UT_RTIVO2("protocols/toolbox/task_operations/RestrictToInterfaceVectorOperation.u");
		RTIVO_factory.clear();
		utility::vector1< int > interface_jump;
		interface_jump.push_back(1);
		RTIVO_factory.push_back( TaskOperationCOP( new RestrictToInterfaceVectorOperation(interface_jump,10,5.5,75,9.0) ) );
		//std::cout <<"Interface Jump RestrictToInterfaceVectorOperation \n "
		//<< *(RTIVO_factory.create_task_and_apply_taskoperations( pose )) << std::endl;
		UT_RTIVO2 << *(RTIVO_factory.create_task_and_apply_taskoperations( pose )) << std::endl;

	}//end test_RestrictToInterfaceVectorOperation

	void test_RestrictInterGroupVectorOperation(){
		//set up test
		using namespace core::pack::task;
		using protocols::toolbox::task_operations::RestrictInterGroupVectorOperation;
		using core::pack::task::operation::TaskOperationCOP;
		TaskFactory RIGV_factory;
		test::UTracer UT_RIGV("protocols/toolbox/task_operations/RestrictInterGroupVector.u");
		//make a few groups
		std::set< core::Size > side1, side2, partA, partB;
		//std::pair< std::set< core::Size >, std::set< core::Size > > interface, otherparts;
		utility::vector1< std::pair< std::set< core::Size >, std::set< core::Size > > > full_vec;
		//for loops to add residues to sets
		for ( core::Size ii=1 ; ii<= 245; ++ii ) {
			side1.insert(ii);
		}
		for ( core::Size ii=246 ; ii<= 301; ++ii ) {
			side2.insert(ii);
		}
		//make other parts by hand
		partA.insert(105);
		partA.insert(107);
		partB.insert(245);
		partB.insert(166);
		//residues 105, 107, and 245 should be included in design, 166 should not.
		std::pair< std::set< core::Size >, std::set< core::Size > > interface( side1, side2 );
		std::pair< std::set< core::Size >, std::set< core::Size > > otherparts( partA, partB );
		full_vec.push_back(interface);
		full_vec.push_back(otherparts);
		//now make the task
		RIGV_factory.push_back( TaskOperationCOP( new RestrictInterGroupVectorOperation( full_vec,10,5.5,75,9.0 ) ) );
		//output
		//std::cout <<"RestrictInterGroupVectorOperation \n "
		//<< *(RIGV_factory.create_task_and_apply_taskoperations( pose )) << std::endl;
		UT_RIGV << *(RIGV_factory.create_task_and_apply_taskoperations( pose )) << std::endl;


	}//end RestrictInterGroupVectorOperation

	void test_RestrictToMoveMapChiOperation() {
		using namespace core::pack::task;
		using namespace core::kinematics;
		using namespace protocols::toolbox::task_operations;
		using core::pack::task::operation::TaskOperationCOP;
		using utility::vector1;

		MoveMapOP mm( new MoveMap() );
		for ( core::Size i = 1; i <=10; ++i ) {
			mm->set_chi(i, true);
		}
		mm->show(TR);
		TaskFactory tf;

		RestrictToMoveMapChiOperationOP mm_op( new RestrictToMoveMapChiOperation(mm) );
		tf.push_back(mm_op);
		PackerTaskOP task = tf.create_task_and_apply_taskoperations(pose);
		task->show(TR);

		vector1<bool> repacking_residues(pose.size(), false);
		for ( core::Size i = 1; i <=10; ++i ) {
			repacking_residues[i] = true;
		}
		TS_ASSERT(repacking_residues == task->repacking_residues());

		//Testing neighbor detection.

		tf.clear();
		mm_op->set_include_neighbors(true);
		mm_op->set_cutoff_distance(10.0);
		tf.push_back(mm_op);
		task = tf.create_task_and_apply_taskoperations(pose);
		task->show(TR);
		test::UTracer UT_MMNEI("protocols/toolbox/task_operations/RestrictToMoveMapChiWNeighbors.u");

		UT_MMNEI << *task << std::endl;

		//Testing design.
		tf.clear();
		mm_op->set_include_neighbors(false);
		mm_op->set_design(true);
		tf.push_back(mm_op);
		task = tf.create_task_and_apply_taskoperations(pose);
		TS_ASSERT(repacking_residues == task->designing_residues());

	}

};
