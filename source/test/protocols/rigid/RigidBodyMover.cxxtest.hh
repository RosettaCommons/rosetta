// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/rigid/RigidBodyMover.cxxtest.hh
/// @brief  test for RigidBodyMover
/// @author Sid Chaudhury

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/UTracer.hh>
#include <test/util/rosettascripts.hh>

// Unit header
#include <protocols/rigid/RigidBodyMover.hh>

// project headers
#include <core/types.hh>

// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Edge.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/RotamerSetBase.fwd.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
#include <core/conformation/symmetry/SymDof.fwd.hh>
#include <core/conformation/symmetry/SymDof.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/DOF_ID_Range.fwd.hh>
#include <core/id/DOF_ID_Range.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/TorsionID.hh>
#include <core/id/TorsionID_Range.fwd.hh>
#include <core/id/TorsionID_Range.hh>
#include <core/id/types.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/file_data.fwd.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/kinematics/Edge.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/EnergyGraph.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/LREnergyContainer.fwd.hh>
#include <core/scoring/MinimizationGraph.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionInfo.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/methods/ContextDependentLRTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextDependentOneBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextDependentTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/EnergyMethod.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/methods/LongRangeTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/TwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/WholeStructureEnergy.fwd.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/canonical_sampling/MetropolisHastingsMover.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/rigid/RigidBodyMover.fwd.hh>
#include <protocols/canonical_sampling/ThermodynamicMover.fwd.hh>
#include <protocols/canonical_sampling/ThermodynamicMover.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/signals/BufferedSignalHub.fwd.hh>
#include <utility/signals/BufferedSignalHub.hh>
#include <utility/signals/Link.fwd.hh>
#include <utility/signals/Link.hh>
#include <utility/signals/LinkUnit.fwd.hh>
#include <utility/signals/LinkUnit.hh>
#include <utility/signals/SignalHub.fwd.hh>
#include <utility/signals/SignalHub.hh>
#include <utility/tag/Tag.fwd.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/sphericalVector.fwd.hh>
#include <numeric/trig.functions.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/xyzVector.hh>
#include <ObjexxFCL/Dimension.fwd.hh>
#include <ObjexxFCL/Dimension.hh>
#include <ObjexxFCL/DimensionExpression.hh>
#include <ObjexxFCL/DynamicIndexRange.fwd.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>
#include <ObjexxFCL/FArray.fwd.hh>
#include <ObjexxFCL/FArray.hh>
#include <ObjexxFCL/FArray1.fwd.hh>
#include <ObjexxFCL/FArray1.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2.fwd.hh>
#include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArrayInitializer.fwd.hh>
#include <ObjexxFCL/FArrayInitializer.hh>
#include <ObjexxFCL/FArraySection.fwd.hh>
#include <ObjexxFCL/FArraySection.hh>
#include <ObjexxFCL/FArrayTraits.fwd.hh>
#include <ObjexxFCL/FArrayTraits.hh>
#include <ObjexxFCL/IndexRange.fwd.hh>
#include <ObjexxFCL/IndexRange.hh>
#include <ObjexxFCL/InitializerSentinel.hh>
#include <ObjexxFCL/Observer.fwd.hh>
#include <ObjexxFCL/Observer.hh>
#include <ObjexxFCL/ObserverMulti.hh>
#include <ObjexxFCL/ObserverSingle.hh>
#include <ObjexxFCL/ProxySentinel.hh>
#include <ObjexxFCL/SetWrapper.fwd.hh>
#include <ObjexxFCL/Star.fwd.hh>
#include <ObjexxFCL/Star.hh>
#include <ObjexxFCL/proxy_const_assert.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <execinfo.h>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>
#include <basic/MetricValue.fwd.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>


// --------------- Test Class --------------- //

class RigidBodyMoverTests : public CxxTest::TestSuite {

	core::pose::Pose pose; //RigidBodyMover.pdb
	core::Size rb_jump; //1

public:

	// --------------- Fixtures --------------- //

	void setUp() {
		core_init();
		core::import_pose::pose_from_pdb( pose, "protocols/rigid/RigidBodyMover.pdb" );
		rb_jump = 1;

		//setting up the fold tree as is used in docking
		core::kinematics::FoldTree fold_tree;

		core::Size jump_pos1 = 197;
		core::Size jump_pos2 = 282;
		core::Size cutpoint = 245;

    fold_tree.clear();
    fold_tree.add_edge( jump_pos1, jump_pos2, rb_jump );
    fold_tree.add_edge( 1, jump_pos1, core::kinematics::Edge::PEPTIDE );
    fold_tree.add_edge( jump_pos1, cutpoint, core::kinematics::Edge::PEPTIDE );
    fold_tree.add_edge( jump_pos2, cutpoint+1, core::kinematics::Edge::PEPTIDE );
    fold_tree.add_edge( jump_pos2, pose.total_residue(), core::kinematics::Edge::PEPTIDE );
    fold_tree.reorder( 1 );

		pose.fold_tree(fold_tree);
		//std::cout << "set_up() fold tree: " << pose.fold_tree() << std::endl;
	}

	void tearDown() {
		pose.clear();
	}

	// ------------- Helper Functions ------------- //


	// --------------- Test Cases --------------- //

	///@brief test a RigidBodyPerturbMover without rot_center calculated at the interface
	void test_RigidBodyPerturbMover_no_int() {

		////////////////////////RBmover///////////////////////////////////////////////
		using protocols::rigid::RigidBodyPerturbMoverOP;
		using protocols::rigid::RigidBodyPerturbMover;
		core::Real rot_mag(3.0);
		core::Real trans_mag(8.0);

		RigidBodyPerturbMoverOP RB_mover = new RigidBodyPerturbMover(rb_jump, rot_mag, trans_mag, protocols::rigid::partner_downstream, false);

		/////////////////////////run
		RB_mover->apply(pose);
		test::UTracer UT("protocols/rigid/RigidBodyPerturbMover_no_int.pdb");
		UT.abs_tolerance(0.003);
		UT << std::endl;
		pose.dump_pdb(UT);

	}//end test_RigidBodyPerturbMover_no_int

	///@brief test a RigidBodyPerturbMover with rot_center calculated at the interface
	void test_RigidBodyPerturbMover_int() {

		////////////////////////RBmover///////////////////////////////////////////////
		using protocols::rigid::RigidBodyPerturbMoverOP;
		using protocols::rigid::RigidBodyPerturbMover;
		core::Real rot_mag(3.0);
		core::Real trans_mag(8.0);

		core::scoring::ScoreFunctionOP scorefxn;
		scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( core::scoring::PRE_TALARIS_2013_STANDARD_WTS, core::scoring::DOCK_PATCH ) ;
		(*scorefxn)(pose);  //sets up EnergyGraph for interface calculation

		RigidBodyPerturbMoverOP RB_mover = new RigidBodyPerturbMover(rb_jump, rot_mag, trans_mag, protocols::rigid::partner_downstream, true);

		/////////////////////////run
		RB_mover->apply(pose);
		test::UTracer UT("protocols/rigid/RigidBodyPerturbMover_int.pdb");
		UT.abs_tolerance(0.003);
		UT << std::endl;
		pose.dump_pdb(UT);

	}//end test_RigidBodyPerturbMover_int

	///@brief test a RigidBodyPerturbNoCenterMover
	void test_RigidBodyPerturbNoCenterMover() {

		////////////////////////RBmover///////////////////////////////////////////////
		using protocols::rigid::RigidBodyPerturbNoCenterMoverOP;
		using protocols::rigid::RigidBodyPerturbNoCenterMover;
		core::Real rot_mag(3.0);
		core::Real trans_mag(8.0);

		RigidBodyPerturbNoCenterMoverOP RB_mover = new RigidBodyPerturbNoCenterMover(rb_jump, rot_mag, trans_mag);

		/////////////////////////run
		RB_mover->apply(pose);
		test::UTracer UT("protocols/rigid/RigidBodyPerturbNoCenterMover.pdb");
		UT.abs_tolerance(0.003);
		UT << std::endl;
		pose.dump_pdb(UT);

	}//end test_RigidBodyPerturbNoCenterMover

	///@brief test a RigidBodyRandomizeMover
	void test_RigidBodyRandomizeMover() {

		////////////////////////RBmover///////////////////////////////////////////////
		using protocols::rigid::RigidBodyRandomizeMoverOP;
		using protocols::rigid::RigidBodyRandomizeMover;

		RigidBodyRandomizeMoverOP RB_mover = new RigidBodyRandomizeMover(pose, rb_jump, protocols::rigid::partner_downstream);

		/////////////////////////run
		RB_mover->apply(pose);
		test::UTracer UT("protocols/rigid/RigidBodyRandomizeMover.pdb");
		UT.abs_tolerance(0.003);
		UT << std::endl;
		pose.dump_pdb(UT);

	}//end test_RigidBodyRandomizeMover

	///@brief test a RigidBodySpinMover
	void test_RigidBodySpinMover() {

		////////////////////////RBmover///////////////////////////////////////////////
		using protocols::rigid::RigidBodySpinMoverOP;
		using protocols::rigid::RigidBodySpinMover;

		RigidBodySpinMoverOP RB_mover = new RigidBodySpinMover(rb_jump);

		/////////////////////////run
		RB_mover->apply(pose);
		test::UTracer UT("protocols/rigid/RigidBodySpinMover.pdb");
		UT.abs_tolerance(0.003);
		UT << std::endl;
		pose.dump_pdb(UT);

	}//end test_RigidBodySpinMover

	///@brief test a RigidBodyTransMover
	void test_RigidBodyTransMover() {

		////////////////////////RBmover///////////////////////////////////////////////
		using protocols::rigid::RigidBodyTransMoverOP;
		using protocols::rigid::RigidBodyTransMover;

		RigidBodyTransMoverOP RB_mover = new RigidBodyTransMover(pose, rb_jump);

		/////////////////////////run
		RB_mover->apply(pose);
		test::UTracer UT("protocols/rigid/RigidBodyTransMover.pdb");
		UT.abs_tolerance(0.003);
		UT << std::endl;
		pose.dump_pdb(UT);

	}//end test_RigidBodyTransMover

	///@brief test parsing RigidBodyTransMover
	void test_RigidBodyTransMover_parse() {

		////////////////////////RBmover///////////////////////////////////////////////
		using protocols::rigid::RigidBodyTransMoverOP;
		using protocols::rigid::RigidBodyTransMover;

		RigidBodyTransMoverOP RB_mover = new RigidBodyTransMover();

		DataMap data;
		Filters_map filters;
		Movers_map movers;
		TagPtr tag = tagptr_from_string("<RigidBodyTransMover name=trans distance=1 jump=1/>\n");
		RB_mover->parse_my_tag( tag, data, filters, movers, pose );

		/////////////////////////run
		RB_mover->apply(pose);
		test::UTracer UT("protocols/rigid/RigidBodyTransMover.pdb");
		UT.abs_tolerance(0.003);
		UT << std::endl;
		pose.dump_pdb(UT);

	}//end test_RigidBodyTransMover_parse()



	///@brief test a RigidBodyTransMover
	//void test_RigidBodyUniformSphereTransMover() {

	//	////////////////////////RBmover///////////////////////////////////////////////
	//	using protocols::rigid::UniformSphereTransMoverOP;
	//	using protocols::rigid::UniformSphereTransMover;

	//	core::Real mag(10.0);
	//	UniformSphereTransMoverOP RB_mover = new UniformSphereTransMover(rb_jump, mag);

	//	/////////////////////////run
	//	RB_mover->apply(pose);
	//	test::UTracer UT("protocols/rigid/RigidBodyUniformSphereTransMover.pdb");
	//	UT.abs_tolerance(0.003);
	//	UT << std::endl;
	//	pose.dump_pdb(UT);

	//}//end test_RigidBodyUniformSphereTransMover

};//end class
