// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/mm/MMLJEnergy.cxxtest.hh
/// @brief  Unit tests for the mm lennard-jones energy
/// @author Andy Watkins (andy.watkins2@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/util/deriv_funcs.hh>

#include <core/types.hh>

// Unit headers
#include <core/scoring/methods/MMLJEnergyIntra.hh>

// Package headers
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

// Project headers
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>

#include <utility/graph/Graph.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/conformation/AbstractRotamerTrie.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/vector1.hh>


using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.scoring.mm.MMLJEnergy.cxxtest");

// using declarations
using namespace std;
using namespace core;
using namespace scoring;
using namespace etable;

///////////////////////////////////////////////////////////////////////////
/// @name MMLJEnergyTest
/// @brief: Test the functionality of the MMLJEnergy class
///////////////////////////////////////////////////////////////////////////
class MMLJEnergyTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void tearDown() {}

	void test_etab_numeric_deriv_check_analytic_version()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( mm_lj_intra_atr, 0.5 );
		sfxn.set_weight( mm_lj_intra_rep, 0.25 );
		kinematics::MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv( pose, sfxn, movemap );

		adv.simple_deriv_check( true, 1e-6 );
	}

};
/*
To be added back in once e.g. setup_for_minimizing_for_residue is enabled for
this energy method in master!

void eval_intra_residue_energy_w_minimization_data_analytic_version()
{
using namespace core::pose;
using namespace core::scoring;
using namespace core::scoring::etable;
using namespace core::scoring::methods;

Pose pose = create_trpcage_ideal_pose();

EnergyMap emap;
ScoreFunction sfxn;
MMLJEnergy mmljenergy;

sfxn.set_weight( mm_lj_intra_atr, 0.5 ); // enable intra-residue energy evaluation; skipped if all intra-residue weights are 0
sfxn.set_weight( mm_lj_intra_rep, 0.25 ); // enable intra-residue energy evaluation; skipped if all intra-residue weights are 0
//mmljenergy.eval_intrares_energy( pose.residue( 4 ), pose, sfxn, emap );

optimization::MinimizerMap minmap;
EnergyMap emap2;
ResSingleMinimizationData min_data;
mmljenergy.setup_for_minimizing_for_residue( pose.residue( 4 ), pose, sfxn,  minmap, min_data );
mmljenergy.eval_intrares_energy_ext( pose.residue( 4 ), min_data, pose, sfxn, emap2 );

TS_ASSERT_DELTA( emap[ fa_intra_atr ], emap2[ fa_intra_atr ], 1e-12 );
TS_ASSERT_DELTA( emap[ fa_intra_rep ], emap2[ fa_intra_rep ], 1e-12 );
TS_ASSERT_DELTA( emap[ fa_intra_sol ], emap2[ fa_intra_sol ], 1e-12 );

}*/

