// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/RestrictToInterfaceOperation.cxxtest.hh
/// @brief  test for RestrictToInterfaceOperation
/// @author Steven Lewis

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/UTracer.hh>

// Unit header
#include <protocols/toolbox/task_operations/RestrictToInterfaceOperation.hh>
#include <protocols/toolbox/task_operations/RestrictToNeighborhoodOperation.hh>
#include <protocols/toolbox/task_operations/RestrictByCalculatorsOperation.hh>
#include <protocols/toolbox/task_operations/RestrictToInterfaceVectorOperation.hh>


// project headers
#include <core/types.hh>

// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
#include <core/pose/Pose.hh>
#include <core/pose/metrics/CalculatorFactory.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>

#include <protocols/toolbox/pose_metric_calculators/NeighborhoodByDistanceCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/InterfaceNeighborDefinitionCalculator.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <utility/vector1.hh>


// C++ headers
//#include <set>

// --------------- Test Class --------------- //

class RestrictOperationsTests : public CxxTest::TestSuite {

	core::pose::Pose pose;

public:

	// --------------- Fixtures --------------- //

	//ctor sets up the pose once
	RestrictOperationsTests(){
		core_init();
		//reuse for comparison with Interface class
		core::import_pose::centroid_pose_from_pdb( pose, "core/conformation/dock_in.pdb" );
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

	///@brief
	void test_RestrictToInterfaceOperation() {

		/* if you are debugging this test, I suggest passing the flag -print_pymol_selection for ease of use.
			 The flag -pose_metrics::interface_cutoff 8 will cause the results to match core/conformation/Interface.cxxtest
			 although the flag actually defaults to 10 (thus a larger interface)
		*/

		//set up test
		using namespace core::pack::task;
		using protocols::toolbox::task_operations::RestrictToInterfaceOperation;
		TaskFactory RTIO_factory;
		RTIO_factory.push_back( new RestrictToInterfaceOperation() ); //defaults to interface between chains 1 and 2

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
		TaskFactory RTNO_factory;
		RTNO_factory.push_back( new RestrictToNeighborhoodOperation( crset ) );

		test::UTracer UT_RTNO("protocols/toolbox/task_operations/RestrictToNeighborhoodOperation.u");
		UT_RTNO << *(RTNO_factory.create_task_and_apply_taskoperations( pose )) << std::endl;
	}

	void test_RestrictByCalculatorsOperation() {

		//first we set up the calculators that the Operation will use
		std::string const interface_calc("interface"), neighborhood_calc("neighborhood");
		std::set< core::Size > crset_RBC;
		crset_RBC.insert(127); crset_RBC.insert(170), crset_RBC.insert(46);
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( interface_calc, new protocols::toolbox::pose_metric_calculators::InterfaceNeighborDefinitionCalculator( core::Size(1), core::Size(2) ) );
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( neighborhood_calc, new protocols::toolbox::pose_metric_calculators::NeighborhoodByDistanceCalculator( crset_RBC ) );

		//this is the constructor parameter for the calculator - pairs of calculators and calculations to perform
		utility::vector1< std::pair< std::string, std::string> > calcs_and_calcns;
		calcs_and_calcns.push_back(std::make_pair(interface_calc, "interface_residues"));
		calcs_and_calcns.push_back(std::make_pair(neighborhood_calc, "neighbors"));

		using protocols::toolbox::task_operations::RestrictByCalculatorsOperation;
		core::pack::task::TaskFactory RBC_factory;
		RBC_factory.push_back( new RestrictByCalculatorsOperation( calcs_and_calcns ) );

		test::UTracer UT_RBC("protocols/toolbox/task_operations/RestrictByCalculatorsOperation.u");
		UT_RBC << *(RBC_factory.create_task_and_apply_taskoperations( pose )) << std::endl;

	}//end test_RestrictByCalculatorsOperation

	void test_RestrictToInterfaceVectorOperation(){

		//set up test
		using namespace core::pack::task;
		using protocols::toolbox::task_operations::RestrictToInterfaceVectorOperation;
		TaskFactory RTIVO_factory;
		//these are the default values but hard code anyway, test for chain #s
		RTIVO_factory.push_back( new RestrictToInterfaceVectorOperation(1,2,10,5.5,75,9.0) );
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
		RTIVO_factory.push_back( new RestrictToInterfaceVectorOperation(interface_jump,10,5.5,75,9.0) );
		//std::cout <<"Interface Jump RestrictToInterfaceVectorOperation \n "
		//<< *(RTIVO_factory.create_task_and_apply_taskoperations( pose )) << std::endl;
		UT_RTIVO2 << *(RTIVO_factory.create_task_and_apply_taskoperations( pose )) << std::endl;

	}//end test_RestrictToInterfaceVectorOperation

};//end class
