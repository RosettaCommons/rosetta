// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/core/scoring/util.cxxtest.hh
/// @brief  Tests for core/scoring/util.hh
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

#include <core/scoring/util.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/hbonds/HBondOptions.hh>

#include <core/types.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>

#include <map>

class ScoringUtilTests : public CxxTest::TestSuite {
public:
	static constexpr core::Real delta_ = 1e-6;

	void setUp() {
		core_init();
	}

	void test_atomistic_energies() {
		using namespace core::scoring;

		core::pose::PoseOP pose = create_test_in_pdb_poseop();
		core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function("ref2015");
		core::scoring::methods::EnergyMethodOptionsOP emopts( new core::scoring::methods::EnergyMethodOptions( scorefxn->energy_method_options() ) );
		emopts->hbond_options().decompose_bb_hb_into_pair_energies( true );
		scorefxn->set_energy_method_options( *emopts );

		// All the energy terms from ref2015 where we know the atomistic decomposability is present,
		// and for which the atom sum should match the residue-level values.
		core::scoring::ScoreTypes known_types = { fa_atr, fa_rep, fa_sol, fa_intra_rep, fa_intra_sol_xover4, lk_ball_wtd, fa_elec, hbond_sr_bb, hbond_lr_bb, hbond_bb_sc, hbond_sc };
		core::Size n_known_types = known_types.size();

		std::map< core::id::AtomID,  utility::vector1< core::Real > > id_energy_map = get_atomistic_energies( *pose, *scorefxn, known_types ); // Over all the pose.

		utility::vector1< core::Real > blank_vec( n_known_types, 0.0 );
		utility::vector1< utility::vector1< core::Real > > per_residue_sum( pose->size(), blank_vec );
		for ( auto const & entry: id_energy_map ) {
			core::Size rsd = entry.first.rsd();
			utility::vector1< core::Real > const & energies = entry.second;
			for ( core::Size ii(1); ii <= known_types.size(); ++ii ) {
				per_residue_sum[rsd][ii] += energies[ii];
			}
		}

		scorefxn->score( *pose );
		for ( core::Size rsd(1); rsd <= pose->size(); ++rsd ) {
			EnergyMap reference = pose->energies().residue_total_energies(rsd) * scorefxn->weights();
			for ( core::Size ii(1); ii <= known_types.size(); ++ii ) {
				ScoreType st = known_types[ii];
				TSM_ASSERT_DELTA( "Residue " + std::to_string(rsd) + " term " + name_from_score_type(st), per_residue_sum[rsd][ii], reference[st], delta_ );
			}
		}
	}

};
