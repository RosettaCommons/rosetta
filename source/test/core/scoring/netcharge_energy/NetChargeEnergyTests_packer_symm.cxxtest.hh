// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/core/scoring/netcharge_energy/NetChargeEnergy.cxxtest.hh
/// @brief  Test suite for core::scoring::netcharge_energy::NetChargeEnergy, an energy term for controlling
/// net charge during design.
/// @details See also the core::conformation::symmetry::MirrorSymmetricConformation unit tests.  These have
/// another example of NetCharge being set up from code (with constraints attached to the pose).
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <core/scoring/netcharge_energy/NetChargeEnergySetup.hh>
#include <core/scoring/netcharge_energy/NetChargeEnergy.hh>

// Unit headers

#include <platform/types.hh>

// Package Headers
#include <test/util/deriv_funcs.hh>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>

#include <core/pose/annotated_sequence.hh>

#include <core/pack/packer_neighbors.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/make_symmetric_task.hh>
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSets.hh>

#include <core/pack/interaction_graph/ResidueArrayAnnealingEvaluator.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/pose/util.hh>
#include <basic/Tracer.hh>

// Convenience for setting up symmetry and packing:
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>

//Auto Headers
#include <utility/vector1.hh>


static basic::Tracer TR("core.scoring.netcharge_energy.NetChargeEnergyTests_packer_symm.cxxtest");

// --------------- Test Class --------------- //

// using declarations
using namespace core;
using namespace core::pose;
using namespace core::scoring;
using namespace core::scoring::methods;
using namespace core::scoring::annealing;
using namespace core::scoring::symmetry;

using namespace core::pack;
using namespace core::pack::task;
using namespace core::pack::rotamer_set;

using namespace protocols::simple_moves::symmetry;

class NetChargeEnergyTests_packer_symm : public CxxTest::TestSuite {

public:

	// --------------- Fixtures --------------- //

	// Shared initialization goes here.
	void setUp() {
		core_init_with_additional_options("-score:netcharge_setup_file core/scoring/netcharge_energy/negative_three.charge" );
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	/// @brief Test the energy calculation with the packer and symmetry.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).  Modified from Alex Ford's tests for the aa_composition energy.
	void test_energy_annealing_symmetric() {
		SymmetricScoreFunctionOP scorefxn( new SymmetricScoreFunction );
		scorefxn->set_weight( netcharge, 1 );

		Pose pose;
		make_pose_from_sequence( pose, "GGGGEKGG", "fa_standard" );
		SetupForSymmetryMover symm( "core/scoring/symmetry/C3.symm" );
		symm.apply(pose); //Symmetrize the pose.
		pose.update_residue_neighbors();

		core::Real const initial_energy( (*scorefxn)(pose) );
		TS_ASSERT_DELTA( initial_energy, 90.0 , 1e-6);

		// Setup packer task and packer objects
		PackerTaskOP task( TaskFactory::create_packer_task( pose ) );

		utility::vector1< bool > keep_aas( core::chemical::num_canonical_aas, false );
		keep_aas[ core::chemical::aa_lys ] = true;
		keep_aas[ core::chemical::aa_asp ] = true;
		keep_aas[ core::chemical::aa_gly ] = true;
		keep_aas[ core::chemical::aa_glu ] = true;

		for ( core::Size i = 1; i <= pose.total_residue(); i++ ) {
			task->nonconst_residue_task( i ).restrict_absent_canonical_aas( keep_aas );
		}

		task->request_symmetrize_by_intersection();
		task = core::pack::make_new_symmetric_PackerTask_by_requested_method( pose, task );

		core::pack::interaction_graph::AnnealableGraphBaseOP symmetric_ig;
		core::pack::rotamer_set::symmetry::SymmetricRotamerSetsOP sym_rotamer_sets( new core::pack::rotamer_set::symmetry::SymmetricRotamerSets );
		utility::vector0< int > rot_to_pack;
		core::pack::pack_rotamers_setup( pose, *scorefxn, task, sym_rotamer_sets, symmetric_ig );

		//protocols::minimization_packing::PackRotamersMover packer( scorefxn, task, 1 );

		for ( core::Size i(1); i<=5; ++i ) {
			// Test via pack_rotamers run five times.
			core::pack::pack_rotamers_run( pose, task, sym_rotamer_sets, symmetric_ig, rot_to_pack );

			//packer.apply(pose);
			core::Real const postpack_score( (*scorefxn)(pose) );
			TS_ASSERT_DELTA(postpack_score, 0, 1e-6);

			TR << "Seq" << i << ":\t" << pose.sequence() << std::endl;

			std::map< std::string, int > aa_count;
			aa_count["ASP"] = 0;
			aa_count["GLU"] = 0;
			aa_count["LYS"] = 0;
			aa_count["GLY"] = 0;

			for ( core::Size r = 1; r<= pose.size(); ++r ) {
				aa_count[ pose.residue(r).name3() ] += 1;
			}

			TS_ASSERT_EQUALS( aa_count["ASP"] + aa_count["GLU"] - aa_count["LYS"], 3);
		}
	}

};
