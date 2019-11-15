// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/elec/FA_ElecEnergy_pre_talaris.cxxtest.hh
/// @brief  Test the score and derivative evaluation for the FA_ElecEnergy class with the pre-Talaris
/// settings and scorefunction.
/// @author Andrew Leaver-Fay
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org) -- Split the pre-Talaris test
/// into its own test suite, since it requires different initialization which was messing up
/// the initialization of the tests of the current scorefunction.

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/UTracer.hh>
#include <test/util/deriv_funcs.hh>
#include <test/util/pose_funcs.hh>
#include <test/util/pdb1ubq.hh>
#include <test/core/init_util.hh>

// Package headers
#include <core/scoring/elec/FA_ElecEnergy.hh>
#include <core/scoring/etable/coulomb/Coulomb.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/EnergyMap.hh>

// Utility headers
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR("core.scoring.elec.FA_ElecEnergy_pre_talaris.cxxtest");

using namespace core;

class FA_ElecEnergyTests_pre_talaris : public CxxTest::TestSuite {

public:
	void setUp() {
		core_init_with_additional_options( "-restore_pre_talaris_2013_behavior -override_rsd_type_limit" );
	}

	void tearDown(){}

	void test_elec_pre_talaris_2013_settings()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;

		core::pose::Pose pose = create_twores_1ubq_pose();

		core::Real const TOL(1e-5);

		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( fa_elec, 1 );

		methods::EnergyMethodOptions options; // default is what we want

		core::scoring::etable::coulomb::Coulomb coulomb( options );
		//options.show(TR);

		EnergyMap emap;

		TS_ASSERT_DELTA( coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("CB") ), -0.18,
			pose.residue(2).xyz(pose.residue(1).atom_index("CB") ), -0.18 ), 0, TOL ); // CB-CB sits at 5.56174
		TS_ASSERT_DELTA( coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("CB") ), -0.18,
			pose.residue(2).xyz(pose.residue(1).atom_index("CA") ), 0.07 ), -0.006643886, TOL ); // CB-CA at 4.49784 (no count pair)
		TS_ASSERT_DELTA( coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("C") ), 0.51,
			pose.residue(2).xyz(pose.residue(1).atom_index("N") ), -0.47 ), -3.175849739, TOL ); // C-N at 1.29914 (no count pair)

		core::scoring::elec::FA_ElecEnergy elec( options );
		elec.residue_pair_energy( pose.residue(1), pose.residue(2), pose, sfxn, emap );
		// Value calculated for 1.5/5.5 min/max, die=10r,
		// with no fa_elec on atoms 1,2,or 3 bonds apart, and a scaling of 0.2 on those 4 bonds apart
		TS_ASSERT_DELTA( emap[ fa_elec ] , 0.337019896208238, TOL);
		TS_ASSERT_DELTA( sfxn(pose), 0.337019896208238, TOL);

	}

};
