// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/rotamer_recovery/RotamerRecovery.cxxtest.hh
/// @brief  Test RotamerRecovery class
/// @author Matthew O'Meara (mattjomeara@gmail.com)


// Test Headers
#include <cxxtest/TestSuite.h>
#include <util/pose_funcs.hh>

// Unit Headers
#include <protocols/rotamer_recovery/RotamerRecovery.hh>
#include <protocols/rotamer_recovery/RotamerRecoveryFactory.hh>

// Project Headers
#include <test/core/init_util.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// C++ Headers
#include <iostream>

//Auto Headers
#include <platform/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/RotamerSetBase.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/pack/rotamer_set/RotamerCouplings.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.fwd.hh>
#include <core/pack/task/IGEdgeReweightContainer.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/ResfileReader.fwd.hh>
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/task/operation/TaskOperation.fwd.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/EnergyGraph.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/LREnergyContainer.fwd.hh>
#include <core/scoring/MinimizationGraph.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionInfo.fwd.hh>
#include <core/scoring/ScoreType.hh>
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
#include <protocols/rotamer_recovery/RRComparer.fwd.hh>
#include <protocols/rotamer_recovery/RRProtocol.fwd.hh>
#include <protocols/rotamer_recovery/RRReporter.fwd.hh>
#include <protocols/rotamer_recovery/RotamerRecovery.fwd.hh>
#include <protocols/rotamer_recovery/RotamerRecoveryCreator.fwd.hh>
#include <protocols/rotamer_recovery/RotamerRecoveryFactory.fwd.hh>
#include <utility/down_cast.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/factory/WidgetRegistrator.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/tag/Tag.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iosfwd>
#include <limits>
#include <list>
#include <map>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>
#include <basic/Tracer.fwd.hh>


static basic::Tracer TR("protocols.rotamer_recovery.RotamerRecovery.cxxtest");

class RotamerRecoveryTests : public CxxTest::TestSuite {

public:

	void
	setUp() {
		using core::scoring::get_score_function;
		using core::pack::task::TaskFactory;
		using core::pack::task::operation::RestrictToRepacking;
		using core::pack::task::operation::TaskOperationCOP;

		core_init_with_additional_options( "-patch_selectors VIRTUAL_BB" );
		pose_1ten_ = fullatom_pose_from_string( pdb_string_1ten() );
		score_function_ = get_score_function();

		TaskFactory task_factory;

		task_factory.push_back( TaskOperationCOP( new RestrictToRepacking ) );
		packer_task_1ten_ = task_factory.create_task_and_apply_taskoperations( pose_1ten_ );

	}

	void test_RotamerRecovery_main() {
		do_test_RotamerRecovery_automatic_construction();
	}

	void
	do_test_RotamerRecovery_automatic_construction() {

		using std::endl;
		using core::Real;
		using protocols::rotamer_recovery::RotamerRecoveryOP;
		using protocols::rotamer_recovery::RotamerRecoveryFactory;
		score_function_->setup_for_scoring(pose_1ten_);

		RotamerRecoveryFactory* factory(RotamerRecoveryFactory::get_instance());
		{
			RotamerRecoveryOP rr(
				factory->get_rotamer_recovery(
				"RRProtocolRTMin", "RRComparerRotBins", "RRReporterSimple"));

			rr->run( pose_1ten_,*score_function_,*packer_task_1ten_);
			rr->show( TR );
			TS_ASSERT_DELTA( rr->recovery_rate(), Real(/*57*/59)/Real(89) , .001 );
		}

		{
			RotamerRecoveryOP rr(
				factory->get_rotamer_recovery(
				"RRProtocolRTMin", "RRComparerAutomorphicRMSD", "RRReporterSimple"));
			rr->run( pose_1ten_,*score_function_,*packer_task_1ten_);
			rr->show( TR );
			TS_ASSERT_DELTA( rr->recovery_rate(), Real(12)/Real(89) , .001 );
		}


	}

private:
	core::pose::Pose pose_1ten_;
	core::scoring::ScoreFunctionOP score_function_;
	core::pack::task::PackerTaskOP packer_task_1ten_;

};
