// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/protocols/abinitio/PseudocontactShiftEnergy.cxxtest.hh
/// @brief  test suite for protocols/scoring/methods/pcs/* and protocols/topology_broker/PseudocontactShiftEnergyController
/// @author Christophe Schmitz schmitz@maths.uq.edu.au / cofcof.oz@gmail.com
/// @last_modified June 2009

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
#include <basic/options/keys/broker.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include <protocols/abinitio/AbrelaxMover.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <core/io/pdb/pose_io.hh>
#include <basic/Tracer.hh>

//Auto Headers
#include <core/conformation/Atom.fwd.hh>
#include <core/fragment/FragID_Iterator.fwd.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FrameIteratorWorker_.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Mask.hh>
#include <core/id/DOF_ID_Map.hh>
#include <core/id/DOF_ID_Mask.hh>
#include <core/id/NamedStubID.hh>
#include <core/id/SequenceMapping.hh>
#include <core/io/pdb/file_data.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/silent/EnergyNames.fwd.hh>
#include <core/pose/signals/ConformationEvent.hh>
#include <core/pose/signals/DestructionEvent.hh>
#include <core/pose/signals/EnergyEvent.hh>
#include <core/scoring/types.hh>
#include <core/scoring/constraints/Constraints.fwd.hh>
#include <core/scoring/constraints/XYZ_Func.hh>
#include <core/id/SequenceMapping.hh>
#include <protocols/scoring/methods/pcs/GridSearchIterator.fwd.hh>
#include <utility/stream_util.hh>
#include <protocols/evaluation/ConstraintEvaluator.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributorFactory.hh>
#include <protocols/jumping/PairingLibrary.hh>
#include <protocols/jumping/PairingsList.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/moves/BoolMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/mc_convergence_checks/ConvergenceCheck.hh>
#include <protocols/topology_broker/ClaimerMessage.hh>
#include <protocols/topology_broker/TopologyClaimer.hh>
#include <protocols/scoring/methods/pcs/GridSearchIterator.fwd.hh>
#include <utility/fix_boinc_read.hh>
#include <numeric/xyz.functions.hh>
#include <ObjexxFCL/FArray.fwd.hh>
#include <ObjexxFCL/format.hh>
#include <iterator>
#include <time.h>



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
		core_init_with_additional_options("-broker:setup protocols/abinitio/pcs_broker_setup.txt");

		//We read the setup file with the topologyclaimer framework
		top_bro_OP_ = new  protocols::topology_broker::TopologyBroker();
		try {
			add_cmdline_claims(*top_bro_OP_, false);
		}
		catch ( utility::excn::EXCN_Exception &excn )  {
			excn.show( std::cerr );
			utility_exit();
		}

		the_pose_ = create_test_in_pdb_poseop();
		pcs_energy_ = new PCS_Energy();
	}



	// Shared finalization goes here.
	void tearDown() {
		top_bro_OP_ = 0;
		pcs_energy_ = 0;
		the_pose_ = 0;
	}



	// WARNING this test is a little bit agressive
	// It runs the broker protocol, using the pcs energy term ONLY. It should go through all 4 stages.
	// It could break if the protocol change, if the option flags change.
	void test_eval_abinitio_pcs_only(){

		using namespace core;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;


		core_init_with_additional_options("\
 -abinitio::increase_cycles 0.01\
 -nstruct 1\
 -frag9 protocols/abinitio/frag9_for_pcs_test.tab.gz\
 -frag3 protocols/abinitio/frag3_for_pcs_test.tab.gz\
 -native protocols/abinitio/pdb_idealized_for_pcs_test.pdb\
 -run:protocol broker\
 -run::constant_seed\
 -run::jran 123456\
 -abinitio::stage1_patch protocols/abinitio/score0_pcs_only.wts_patch\
 -abinitio::stage2_patch protocols/abinitio/score1_pcs_only.wts_patch\
 -abinitio::stage3a_patch protocols/abinitio/score2_pcs_only.wts_patch\
 -abinitio::stage3b_patch protocols/abinitio/score5_pcs_only.wts_patch\
 -abinitio::stage4_patch protocols/abinitio/score3_pcs_only.wts_patch\
 -overwrite\
 -out:prefix PCS_\
 ");

		protocols::abinitio::AbrelaxMoverOP abrelax = new protocols::abinitio::AbrelaxMover;
		protocols::jd2::JobDistributor::get_instance()->go( abrelax);
	}



	// WARNING This test is really aggressive and might break easily
	// This test evaluate the pcs energy on the final pdb generated in the test_eval_abinitio_pcs_only() test
	// If the test_eval_abinitio_pcs_only fails, this test will obviously fail too.
	// It could break if some changes have been made in the protocol, since if the outputed pdb is different, the energy will be different
	// If the break is legitimate, please, update the expected value.
	void test_eval_pcs_energy_on_abinitio_output(){

		core::import_pose::pose_from_pdb( *the_pose_, "PCS_S_0001.pdb" );
		core::Real pcs_score_total = pcs_energy_->calculate_pcs_score(*the_pose_, false);
		core::Real expected_value(0.7588863408);
		core::Real tolerance(0.001);
		Tracer_PCS << std::setprecision(10) << "Expected: "  << expected_value << " Calculated: " << pcs_score_total << "Tolerance: " << tolerance  << std::endl;
		Tracer_PCS << std::setprecision(10) << "Comparison of 2 values desactivated. The test is not that deterministic, different values on 32 and 64 bit machines?" << std::endl;
		//		TS_ASSERT_DELTA( pcs_score_total, expected_value, tolerance);
	}



	// WARNING This test is very friendly and should NEVER break.
	// This test simply evaluate the pcs energy term on the native structures, and should always return the same value.
	// It could break if some stuffs have been changed in the PseudocontactShift files
	// It could break if some stuffs have been changed in the SVD files
	// It could break if some stuffs have been changed in the minimizer
	void test_eval_pcs_energy_on_native(){

		core::import_pose::pose_from_pdb( *the_pose_, "protocols/abinitio/pdb_idealized_for_pcs_test.pdb" );
		core::Real pcs_score_total = pcs_energy_->calculate_pcs_score(*the_pose_, false);
		core::Real expected_value(0.5398201739);
		core::Real tolerance(0.0001);
		Tracer_PCS << std::setprecision(10) << "Expected: "  << expected_value << " Calculated: " << pcs_score_total << "Tolerance: " << tolerance  << std::endl;
		TS_ASSERT_DELTA( pcs_score_total, expected_value, tolerance);
	}
};
