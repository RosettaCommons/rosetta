// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/core/scoring/methods/ProClosureEnergy.cxxtest.hh
/// @brief  test suite for core::scoring::ProClosureEnergy.cc
/// @author Andrew Leaver-Fay

// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <core/scoring/methods/CartesianBondedEnergy.hh>

#include <platform/types.hh>

// Package Headers
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <test/util/pose_funcs.hh>
#include <test/util/deriv_funcs.hh>
#include <test/core/init_util.hh>
#include <core/kinematics/DomainMap.hh>

//Auto Headers
#include <utility/vector1.hh>

#include <basic/Tracer.hh>

// --------------- Test Class --------------- //

// using declarations
using namespace core;
using namespace core::id;
using namespace core::kinematics;
using namespace core::pose;
using namespace core::scoring;
using namespace core::scoring::methods;

static basic::Tracer TR("core.scoring.methods.CartesianBondedEnergy.cxxtest");

class CartesianBondedEnergyTests : public CxxTest::TestSuite {

	public:

	// Shared initialization goes here.
	void setUp() {
		using namespace core;
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {}

	// --------------- Test Cases --------------- //
	void test_eval_energy()
	{
		Pose trpcage( create_trpcage_ideal_pose() );

		EnergyMethodOptions opts;
		CartesianBondedEnergy cartbond_energy( opts );
		ScoreFunction sfxn;
		sfxn.set_weight( cart_bonded, 0.5 );


		/*Size before_precision = std::cout.precision();
		std::cout.precision( 16 );
		for ( Size ii = 2; ii <= trpcage.total_residue(); ++ii ) {

			EnergyMap emap;
			std::cout << "{\n"
			"EnergyMap emap;\n";

			cartbond_energy.residue_pair_energy( trpcage.residue( ii-1 ), trpcage.residue( ii ), trpcage, sfxn, emap );
			std::cout << "cartbond_energy.residue_pair_energy( trpcage.residue( " << ii-1 << " ), trpcage.residue( " << ii << " ), trpcage, sfxn, emap );\n";
			std::cout << "TS_ASSERT_DELTA( emap[ cart_bonded ], " << emap[ cart_bonded ] << ", 1e-12 );\n";
			std::cout << "}\n";
		}
		std::cout.precision( before_precision );*/

		{
			EnergyMap emap;
			cartbond_energy.residue_pair_energy( trpcage.residue( 1 ), trpcage.residue( 2 ), trpcage, sfxn, emap );
			TS_ASSERT_DELTA( emap[ cart_bonded ], 6.969547693227086, 1e-12 );
		}
		{
			EnergyMap emap;
			cartbond_energy.residue_pair_energy( trpcage.residue( 2 ), trpcage.residue( 3 ), trpcage, sfxn, emap );
			TS_ASSERT_DELTA( emap[ cart_bonded ], 0.8857209065751711, 1e-12 );
		}
		{
			EnergyMap emap;
			cartbond_energy.residue_pair_energy( trpcage.residue( 3 ), trpcage.residue( 4 ), trpcage, sfxn, emap );
			TS_ASSERT_DELTA( emap[ cart_bonded ], 3.105841396648269, 1e-12 );
		}
		{
			EnergyMap emap;
			cartbond_energy.residue_pair_energy( trpcage.residue( 4 ), trpcage.residue( 5 ), trpcage, sfxn, emap );
			TS_ASSERT_DELTA( emap[ cart_bonded ], 0.2924107074751637, 1e-12 );
		}
		{
			EnergyMap emap;
			cartbond_energy.residue_pair_energy( trpcage.residue( 5 ), trpcage.residue( 6 ), trpcage, sfxn, emap );
			TS_ASSERT_DELTA( emap[ cart_bonded ], 1.766149336774214, 1e-12 );
		}
		{
			EnergyMap emap;
			cartbond_energy.residue_pair_energy( trpcage.residue( 6 ), trpcage.residue( 7 ), trpcage, sfxn, emap );
			TS_ASSERT_DELTA( emap[ cart_bonded ], 0.9622667976521625, 1e-12 );
		}
		{
			EnergyMap emap;
			cartbond_energy.residue_pair_energy( trpcage.residue( 7 ), trpcage.residue( 8 ), trpcage, sfxn, emap );
			TS_ASSERT_DELTA( emap[ cart_bonded ], 0.8804912704560846, 1e-12 );
		}
		{
			EnergyMap emap;
			cartbond_energy.residue_pair_energy( trpcage.residue( 8 ), trpcage.residue( 9 ), trpcage, sfxn, emap );
			TS_ASSERT_DELTA( emap[ cart_bonded ], 0.1803429928439102, 1e-12 );
		}
		{
			EnergyMap emap;
			cartbond_energy.residue_pair_energy( trpcage.residue( 9 ), trpcage.residue( 10 ), trpcage, sfxn, emap );
			TS_ASSERT_DELTA( emap[ cart_bonded ], 1.10972814373643, 1e-12 );
		}
		{
			EnergyMap emap;
			cartbond_energy.residue_pair_energy( trpcage.residue( 10 ), trpcage.residue( 11 ), trpcage, sfxn, emap );
			TS_ASSERT_DELTA( emap[ cart_bonded ], 0.03877637220795598, 1e-12 );
		}
		{
			EnergyMap emap;
			cartbond_energy.residue_pair_energy( trpcage.residue( 11 ), trpcage.residue( 12 ), trpcage, sfxn, emap );
			TS_ASSERT_DELTA( emap[ cart_bonded ], 1.890381532828944, 1e-12 );
		}
		{
			EnergyMap emap;
			cartbond_energy.residue_pair_energy( trpcage.residue( 12 ), trpcage.residue( 13 ), trpcage, sfxn, emap );
			TS_ASSERT_DELTA( emap[ cart_bonded ], 1.44461258712037, 1e-12 );
		}
		{
			EnergyMap emap;
			cartbond_energy.residue_pair_energy( trpcage.residue( 13 ), trpcage.residue( 14 ), trpcage, sfxn, emap );
			TS_ASSERT_DELTA( emap[ cart_bonded ], 0.9956672637671763, 1e-12 );
		}
		{
			EnergyMap emap;
			cartbond_energy.residue_pair_energy( trpcage.residue( 14 ), trpcage.residue( 15 ), trpcage, sfxn, emap );
			TS_ASSERT_DELTA( emap[ cart_bonded ], 15.63800101544166, 1e-12 );
		}
		{
			EnergyMap emap;
			cartbond_energy.residue_pair_energy( trpcage.residue( 15 ), trpcage.residue( 16 ), trpcage, sfxn, emap );
			TS_ASSERT_DELTA( emap[ cart_bonded ], 0.03798202544109427, 1e-12 );
		}
		{
			EnergyMap emap;
			cartbond_energy.residue_pair_energy( trpcage.residue( 16 ), trpcage.residue( 17 ), trpcage, sfxn, emap );
			TS_ASSERT_DELTA( emap[ cart_bonded ], 4.91119354198335, 1e-12 );
		}
		{
			EnergyMap emap;
			cartbond_energy.residue_pair_energy( trpcage.residue( 17 ), trpcage.residue( 18 ), trpcage, sfxn, emap );
			TS_ASSERT_DELTA( emap[ cart_bonded ], 4.143634791534636, 1e-12 );
		}
		{
			EnergyMap emap;
			cartbond_energy.residue_pair_energy( trpcage.residue( 18 ), trpcage.residue( 19 ), trpcage, sfxn, emap );
			TS_ASSERT_DELTA( emap[ cart_bonded ], 3.878703682394249, 1e-12 );
		}
		{
			EnergyMap emap;
			cartbond_energy.residue_pair_energy( trpcage.residue( 19 ), trpcage.residue( 20 ), trpcage, sfxn, emap );
			TS_ASSERT_DELTA( emap[ cart_bonded ], 4.616681252917811, 1e-12 );
		}
	}

	void test_cartbonded_start_score_start_func_match_w_total_flexibility()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( cart_bonded, 0.5 );
		kinematics::MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv;
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );
		adv.validate_start_func_matches_start_score( 26.87406665551287, true );

	}

	void test_cartbonded_deriv_check_w_total_flexibility()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( cart_bonded, 0.5 );
		kinematics::MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv;
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );
		adv.simple_deriv_check( false, 1e-6 );

	}

	void fail_test_create_parameters_for_restype()
	{
		core::chemical::ResidueTypeSetCAP rs;
		std::string rss = core::chemical::FA_STANDARD;

		EnergyMethodOptions opts;
		core::Real cartbonded_len;
		core::Real cartbonded_ang;
		core::Real cartbonded_tors;
		core::Real cartbonded_proton;
		core::Real cartbonded_improper;
		opts.get_cartesian_bonded_parameters(
			cartbonded_len,
			cartbonded_ang,
			cartbonded_tors,
			cartbonded_proton ,
			cartbonded_improper );
		IdealParametersDatabase ipd(
			cartbonded_len,
			cartbonded_ang,
			cartbonded_tors,
			cartbonded_proton ,
			cartbonded_improper );

		rs = core::chemical::ChemicalManager::get_instance()->residue_type_set(rss);
		for(
			core::chemical::ResidueTypeSet::const_residue_iterator
				t=rs->all_residues_begin(), te=rs->all_residues_end();
			t != te; ++t
		) {
			TR << t->second->name() << std::endl;
			{
				bool prepro(true);
				core::scoring::methods::ResidueCartBondedParameters const & rcbp(
					ipd.parameters_for_restype(*(t->second), prepro));
				rcbp.bb_N_index();
				rcbp.bb_CA_index();
				rcbp.bb_O_index();
				rcbp.bb_H_index();
			}

			{
				bool prepro(false);
				core::scoring::methods::ResidueCartBondedParameters const & rcbp(
					ipd.parameters_for_restype(*(t->second), prepro));
				rcbp.bb_N_index();
				rcbp.bb_CA_index();
				rcbp.bb_O_index();
				rcbp.bb_H_index();
			}
		}
	}

};


