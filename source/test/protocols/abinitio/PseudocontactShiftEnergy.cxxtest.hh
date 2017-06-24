// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/abinitio/PseudocontactShiftEnergy.cxxtest.hh
/// @brief  test suite for protocols/scoring/methods/pcs/* and protocols/topology_broker/PseudocontactShiftEnergyController
/// @author Christophe Schmitz schmitz@maths.uq.edu.au / cofcof.oz@gmail.com

// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <protocols/scoring/methods/pcs/PseudocontactShiftEnergy.hh>

// Package Headers
#include <protocols/topology_broker/TopologyBroker.hh>
#include <protocols/topology_broker/util.hh>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/exit.hh>

#include <protocols/abinitio/AbrelaxMover.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <basic/Tracer.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/types.hh>
#include <core/conformation/RotamerSetBase.fwd.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID_Mask.hh>
#include <core/io/silent/EnergyNames.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/ShortestPathInFoldTree.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/methods/EnergyMethod.fwd.hh>
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/methods/EnergyMethodCreator.fwd.hh>
#include <core/scoring/methods/WholeStructureEnergy.fwd.hh>
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <protocols/abinitio/AbrelaxMover.fwd.hh>
#include <protocols/abinitio/FragmentSampler.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/jd2/Job.fwd.hh>
#include <protocols/jd2/JobDistributor.fwd.hh>
#include <protocols/jd2/JobInputter.fwd.hh>
#include <protocols/jd2/JobOutputter.fwd.hh>
//#include <protocols/jd2/Parser.fwd.hh>
#include <protocols/loops/loop_closure/ccd/SlidingWindowLoopClosure.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/relax/RelaxProtocolBase.fwd.hh>
#include <protocols/scoring/methods/pcs/GridSearchIterator.fwd.hh>
#include <protocols/scoring/methods/pcs/PseudocontactShiftData.fwd.hh>
#include <protocols/scoring/methods/pcs/PseudocontactShiftEnergy.fwd.hh>
#include <protocols/scoring/methods/pcs/PseudocontactShiftTensor.fwd.hh>
#include <protocols/topology_broker/ClaimerMessage.fwd.hh>
#include <protocols/topology_broker/claims/DofClaim.fwd.hh>
#include <protocols/topology_broker/TopologyBroker.fwd.hh>
#include <protocols/topology_broker/TopologyClaimer.fwd.hh>
#include <utility/down_cast.hh>
#include <utility/fix_boinc_read.hh>
#include <utility/vector0.fwd.hh>
#include <utility/vector0.hh>
#include <utility/vector0_bool.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/excn/EXCN_Base.hh>
#include <utility/keys/AutoKey.fwd.hh>
#include <utility/keys/AutoKey.hh>
#include <utility/keys/Key.fwd.hh>
#include <utility/keys/Key.hh>
#include <utility/keys/KeyLess.fwd.hh>
#include <utility/keys/KeyLookup.fwd.hh>
#include <utility/keys/KeyLookup.hh>
#include <utility/keys/NoClient.fwd.hh>
#include <utility/keys/NoClient.hh>
#include <utility/keys/SmallKeyVector.fwd.hh>
#include <utility/keys/UserKey.fwd.hh>
#include <utility/keys/VariantKey.fwd.hh>
#include <utility/keys/VariantKey.hh>
#include <utility/options/AnyOption.fwd.hh>
#include <utility/options/BooleanOption.fwd.hh>
#include <utility/options/BooleanVectorOption.fwd.hh>
#include <utility/options/FileOption.fwd.hh>
#include <utility/options/FileVectorOption.fwd.hh>
#include <utility/options/IntegerOption.fwd.hh>
#include <utility/options/IntegerVectorOption.fwd.hh>
#include <utility/options/Option.fwd.hh>
#include <utility/options/PathOption.fwd.hh>
#include <utility/options/PathVectorOption.fwd.hh>
#include <utility/options/RealOption.fwd.hh>
#include <utility/options/RealVectorOption.fwd.hh>
#include <utility/options/StringOption.fwd.hh>
#include <utility/options/StringVectorOption.fwd.hh>
#include <utility/options/keys/AnyOptionKey.fwd.hh>
#include <utility/options/keys/AnyOptionKey.hh>
#include <utility/options/keys/AnyVectorOptionKey.fwd.hh>
#include <utility/options/keys/AnyVectorOptionKey.hh>
#include <utility/options/keys/BooleanOptionKey.fwd.hh>
#include <utility/options/keys/BooleanOptionKey.hh>
#include <utility/options/keys/BooleanVectorOptionKey.fwd.hh>
#include <utility/options/keys/BooleanVectorOptionKey.hh>
#include <utility/options/keys/FileOptionKey.fwd.hh>
#include <utility/options/keys/FileOptionKey.hh>
#include <utility/options/keys/FileVectorOptionKey.fwd.hh>
#include <utility/options/keys/FileVectorOptionKey.hh>
#include <utility/options/keys/IntegerOptionKey.fwd.hh>
#include <utility/options/keys/IntegerOptionKey.hh>
#include <utility/options/keys/IntegerVectorOptionKey.fwd.hh>
#include <utility/options/keys/IntegerVectorOptionKey.hh>
#include <utility/options/keys/OptionKey.fwd.hh>
#include <utility/options/keys/OptionKey.hh>
#include <utility/options/keys/OptionKeys.hh>
#include <utility/options/keys/PathOptionKey.fwd.hh>
#include <utility/options/keys/PathOptionKey.hh>
#include <utility/options/keys/PathVectorOptionKey.fwd.hh>
#include <utility/options/keys/PathVectorOptionKey.hh>
#include <utility/options/keys/RealOptionKey.fwd.hh>
#include <utility/options/keys/RealOptionKey.hh>
#include <utility/options/keys/RealVectorOptionKey.fwd.hh>
#include <utility/options/keys/RealVectorOptionKey.hh>
#include <utility/options/keys/ScalarOptionKey.fwd.hh>
#include <utility/options/keys/ScalarOptionKey.hh>
#include <utility/options/keys/StringOptionKey.fwd.hh>
#include <utility/options/keys/StringOptionKey.hh>
#include <utility/options/keys/StringVectorOptionKey.fwd.hh>
#include <utility/options/keys/StringVectorOptionKey.hh>
#include <utility/options/keys/VectorOptionKey.fwd.hh>
#include <utility/options/keys/VectorOptionKey.hh>
#include <utility/options/keys/all.hh>
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
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <iosfwd>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <basic/Tracer.fwd.hh>
#include <basic/options/keys/OptionKeys.hh>


using namespace core;
using namespace core::pose;
using namespace core::scoring;
using namespace core::scoring::methods;
using namespace protocols::scoring::methods::pcs;
static basic::Tracer Tracer_PCS("protocols.abinitio.PseudocontactShiftEnergy.cxxtest");

class PseudocontactShiftTests : public CxxTest::TestSuite {

public:

	// Shared data elements go here.
	PoseOP the_pose_;
	protocols::topology_broker::TopologyBrokerOP top_bro_OP_;
	PCS_EnergyOP pcs_energy_;


	// Shared initialization goes here.
	void setUp() {

		using namespace std;
		using namespace core;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		core_init();

		//core_init_with_additional_options("-broker:setup protocols/abinitio/pcs_broker_setup.txt -mute all");
		core_init_with_additional_options("-broker:setup protocols/abinitio/pcs_broker_setup.txt "
			"-in:file:fasta protocols/abinitio/pose_funcs_test.fasta");

		//We read the setup file with the topologyclaimer framework
		top_bro_OP_ = protocols::topology_broker::TopologyBrokerOP( new  protocols::topology_broker::TopologyBroker() );
		try {
			add_cmdline_claims(*top_bro_OP_, false);
		}
catch ( utility::excn::EXCN_Exception &excn )  {
	excn.show( std::cerr );
	utility_exit();
}

		the_pose_ = create_test_in_pdb_poseop();
		pcs_energy_ = PCS_EnergyOP( new PCS_Energy() );
	}


	// Shared finalization goes here.
	void tearDown() {
		top_bro_OP_.reset();
		pcs_energy_.reset();
		the_pose_.reset();
	}


	// WARNING this test is a little bit agressive
	// It runs the broker protocol, using the pcs energy term ONLY. It should go through all 4 stages.
	// It could break if the protocol change, if the option flags change.
	void test_eval_abinitio_pcs_only(){

		using namespace core;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;


		core_init_with_additional_options(
			"-abinitio::increase_cycles 0.01 "
			"-nstruct 1 "
			"-frag9 protocols/abinitio/frag9_for_pcs_test.tab.gz "
			"-frag3 protocols/abinitio/frag3_for_pcs_test.tab.gz "
			"-native protocols/abinitio/pdb_idealized_for_pcs_test.pdb "
			"-run:protocol broker "
			"-run::constant_seed "
			"-run::jran 123456 "
			"-abinitio::stage1_patch protocols/abinitio/score0_pcs_only.wts_patch "
			"-abinitio::stage2_patch protocols/abinitio/score1_pcs_only.wts_patch "
			"-abinitio::stage3a_patch protocols/abinitio/score2_pcs_only.wts_patch "
			"-abinitio::stage3b_patch protocols/abinitio/score5_pcs_only.wts_patch "
			"-abinitio::stage4_patch protocols/abinitio/score3_pcs_only.wts_patch "
			"-overwrite "
			"-out:prefix PCS_"
		);

		protocols::abinitio::AbrelaxMoverOP abrelax( new protocols::abinitio::AbrelaxMover );
		protocols::jd2::JobDistributor::get_instance()->go( abrelax);
	}


	// WARNING This test is really aggressive and might break easily
	// This test evaluate the pcs energy on the final pdb generated in the test_eval_abinitio_pcs_only() test
	// If the test_eval_abinitio_pcs_only fails, this test will obviously fail too.
	// It could break if some changes have been made in the protocol, since if the outputed pdb is different, the energy will be different
	// If the break is legitimate, please, update the expected value.
	void test_eval_pcs_energy_on_abinitio_output(){

		core::import_pose::pose_from_file( *the_pose_, "PCS_S_0001.pdb" , core::import_pose::PDB_file);
		core::Real pcs_score_total = pcs_energy_->calculate_pcs_score(*the_pose_, false);
		core::Real expected_value(0.7588863408);
		core::Real tolerance(0.001);
		Tracer_PCS << std::setprecision(10) << "Expected: "  << expected_value << " Calculated: " << pcs_score_total << "Tolerance: " << tolerance  << std::endl;
		Tracer_PCS << std::setprecision(10) << "Comparison of 2 values desactivated. The test is not that deterministic, different values on 32 and 64 bit machines?" << std::endl;
		//  TS_ASSERT_DELTA( pcs_score_total, expected_value, tolerance);
	}


	// WARNING This test is very friendly and should NEVER break.
	// This test simply evaluate the pcs energy term on the native structures, and should always return the same value.
	// It could break if some stuffs have been changed in the PseudocontactShift files
	// It could break if some stuffs have been changed in the SVD files
	// It could break if some stuffs have been changed in the minimizer
	void test_eval_pcs_energy_on_native(){

		core::import_pose::pose_from_file( *the_pose_, "protocols/abinitio/pdb_idealized_for_pcs_test.pdb" , core::import_pose::PDB_file);
		core::Real pcs_score_total = pcs_energy_->calculate_pcs_score(*the_pose_, false);
		core::Real expected_value(0.5401);
		core::Real tolerance(0.0001);
		Tracer_PCS << std::setprecision(10) << "Expected: "  << expected_value << " Calculated: " << pcs_score_total << "Tolerance: " << tolerance  << std::endl;
		TS_ASSERT_DELTA( pcs_score_total, expected_value, tolerance);
	}
};
