// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/hbonds/HBondEnergy.cxxtest.hh
/// @brief  Test the score and derivative evaluation for the HBondEnergy class
/// @author Andrew Leaver-Fay


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/UTracer.hh>
#include <test/util/deriv_funcs.hh>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Package headers
// AUTO-REMOVED #include <core/scoring/hbonds/constants.hh>
// AUTO-REMOVED #include <core/scoring/hbonds/hbonds_geom.hh>
// AUTO-REMOVED #include <core/scoring/hbonds/types.hh>
// AUTO-REMOVED #include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
// AUTO-REMOVED #include <core/scoring/hbonds/HBondTypeManager.hh>

#include <core/scoring/methods/EnergyMethodOptions.hh>

// Project headers
// AUTO-REMOVED #include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/kinematics/MoveMap.hh>
// AUTO-REMOVED #include <core/optimization/MinimizerMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>

// AUTO-REMOVED #include <basic/Tracer.hh>

// AUTO-REMOVED #include <core/id/DOF_ID.hh>
// AUTO-REMOVED #include <core/id/TorsionID.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/scoring/Energies.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
// AUTO-REMOVED #include <basic/options/option.hh>


// Utility headers
#include <utility/vector1.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/types.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/optimization/AtomTreeMinimizer.fwd.hh>
#include <core/optimization/MinimizerOptions.fwd.hh>
#include <core/optimization/NumericalDerivCheckResult.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/SecondaryStructureWeights.hh>
#include <core/scoring/hbonds/HBondOptions.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/mm/MMBondAngleResidueTypeParamSet.fwd.hh>
#include <utility/down_cast.hh>
#include <utility/vector1.fwd.hh>
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
#include <numeric/xyzVector.fwd.hh>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iosfwd>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

//Auto using namespaces
namespace std { } using namespace std; // AUTO USING NS
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
namespace ObjexxFCL { namespace fmt { } } using namespace ObjexxFCL::fmt; // AUTO USING NS
//Auto using namespaces end


using namespace core;
  using namespace conformation;
  using namespace scoring;
    using namespace hbonds;


class HBondEnergyTests : public CxxTest::TestSuite {

public:
  void setUp() {
		core_init();
	}

	void tearDown(){}

	/// @brief Setup for minimization using a move map that says "minimize all bb and sc torsions"
	/// Make sure that start_score matches start_func.
	void test_hbonds_deriv_check_w_full_torsional_flexibility()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn;
		//sfxn.set_weight( fa_atr, 0.5 ); /// just test hbonds; ignore etable
		//sfxn.set_weight( fa_rep, 0.5 );
		sfxn.set_weight( hbond_sr_bb, 0.75 );
		sfxn.set_weight( hbond_lr_bb, 0.5  );
		sfxn.set_weight( hbond_bb_sc, 1.25 );
		sfxn.set_weight( hbond_sc, 1.125 );

		kinematics::MoveMap movemap;
		movemap.set_bb( true );
		movemap.set_chi( true );

		AtomDerivValidator adv( pose, sfxn, movemap );

		/// This call runs a numeric deriv check on all the free dofs in the system and makes sure
		/// that the analytic norm matches the numeric norm to 1e-4.
		adv.simple_deriv_check( true, 1e-6 );

	}

	void test_hbonds_w_sp2term_deriv_check_w_full_torsional_flexibility()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn;
		EnergyMethodOptionsOP emopts( new EnergyMethodOptions( sfxn.energy_method_options() ));
		emopts->hbond_options().use_sp2_chi_penalty( true );
		sfxn.set_energy_method_options( *emopts );
		//sfxn.set_weight( fa_atr, 0.5 ); /// just test hbonds; ignore etable
		//sfxn.set_weight( fa_rep, 0.5 );
		sfxn.set_weight( hbond_sr_bb, 0.75 );
		sfxn.set_weight( hbond_lr_bb, 0.5  );
		sfxn.set_weight( hbond_bb_sc, 1.25 );
		sfxn.set_weight( hbond_sc, 1.125 );

		kinematics::MoveMap movemap;
		movemap.set_bb( true );
		movemap.set_chi( true );

		AtomDerivValidator adv( pose, sfxn, movemap );

		/// This call runs a numeric deriv check on all the free dofs in the system and makes sure
		/// that the analytic norm matches the numeric norm to 1e-4.
		adv.simple_deriv_check( true, 1e-6 );

	}

	/// @brief This cannot be made a unit test, because it is very unlikely that multiple machines
	/// would minimize to the same structure and score.  There is simply too much numerical noise
	/// in the process of minimization
	void dont_test_atom_tree_minimize_with_hbonds_and_etable()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn;
		sfxn.set_weight( hbond_sr_bb, 0.75 );
		sfxn.set_weight( hbond_lr_bb, 0.5  );
		sfxn.set_weight( hbond_bb_sc, 1.25 );
		sfxn.set_weight( hbond_sc, 1.125 );

		AtomTreeMinimizer minimizer;
		Real start_score = sfxn(pose);
		std::cout.precision( 16 );
		std::cout << "start score: " << sfxn(pose) << std::endl;
		TS_ASSERT_DELTA( -33.67974148760877, start_score, 1e-12 );


		kinematics::MoveMap movemap;
		movemap.set_bb( true );
		movemap.set_chi( true );

		MinimizerOptions min_options( "dfpmin_armijo", 0.01, true, true, false );

		minimizer.run( pose, movemap, sfxn, min_options );

		Real end_score = sfxn(pose);
		std::cout << "end score: " << sfxn(pose) << std::endl;
		TS_ASSERT_DELTA( -36.48274046926093, end_score, 1e-12 );
	}

	/// @brief Create a move map that has the central residue changing, while everything
	/// else is held fixed; then the domain map will have color 1 for residues 1-9, and
	/// color 2 for residues 11-20.
	void test_hbonds_deriv_check_w_partial_torsional_flexibility()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn;
		EnergyMethodOptionsOP emopts( new EnergyMethodOptions( sfxn.energy_method_options() ));
		emopts->hbond_options().decompose_bb_hb_into_pair_energies( true );
		sfxn.set_energy_method_options( *emopts );

		sfxn.set_weight( hbond_sr_bb, 0.75 );
		sfxn.set_weight( hbond_lr_bb, 0.5  );
		sfxn.set_weight( hbond_bb_sc, 1.25 );
		sfxn.set_weight( hbond_sc, 1.125 );

		kinematics::MoveMap movemap;
		movemap.set_bb( 10, true );
		AtomDerivValidator adv( pose, sfxn, movemap );
		adv.simple_deriv_check( true, 1e-7 );
	}

	/// @brief Create a move map that has the central residue changing, while everything
	/// else is held fixed; then the domain map will have color 1 for residues 1-9, and
	/// color 2 for residues 11-20.
	void test_hbonds_w_sp2term_deriv_check_w_partial_torsional_flexibility()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn;
		EnergyMethodOptionsOP emopts( new EnergyMethodOptions( sfxn.energy_method_options() ));
		emopts->hbond_options().decompose_bb_hb_into_pair_energies( true );
		emopts->hbond_options().use_sp2_chi_penalty( true );
		sfxn.set_energy_method_options( *emopts );

		sfxn.set_weight( hbond_sr_bb, 0.75 );
		sfxn.set_weight( hbond_lr_bb, 0.5  );
		sfxn.set_weight( hbond_bb_sc, 1.25 );
		sfxn.set_weight( hbond_sc, 1.125 );

		kinematics::MoveMap movemap;
		movemap.set_bb( 10, true );
		AtomDerivValidator adv( pose, sfxn, movemap );
		adv.simple_deriv_check( true, 1e-7 );
	}

	void dont_test_atom_tree_minimize_with_hbonds_and_etable2()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn;
		sfxn.set_weight( fa_atr, 0.5 );
		sfxn.set_weight( fa_rep, 0.5 );
		sfxn.set_weight( hbond_sr_bb, 0.75 );
		sfxn.set_weight( hbond_lr_bb, 0.5  );
		sfxn.set_weight( hbond_bb_sc, 1.25 );
		sfxn.set_weight( hbond_sc, 1.125 );

		AtomTreeMinimizer minimizer;
		Real start_score = sfxn(pose);
		std::cout.precision( 16 );
		std::cout << "start score: " << sfxn(pose) << std::endl;
		TS_ASSERT_DELTA( -33.67974148760877, start_score, 1e-12 );


		kinematics::MoveMap movemap;
		movemap.set_bb( 10, true );

		MinimizerOptions min_options( "dfpmin_armijo", 0.01, true, false, false );

		minimizer.run( pose, movemap, sfxn, min_options );

		Real end_score = sfxn(pose);
		std::cout << "end score: " << sfxn(pose) << std::endl;
		TS_ASSERT_DELTA( -33.88566264995315, end_score, 1e-12 );
	}


};
