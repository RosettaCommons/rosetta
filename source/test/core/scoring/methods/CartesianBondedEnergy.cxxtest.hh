// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
#include <test/util/pose_funcs.hh>
#include <test/util/deriv_funcs.hh>
#include <test/core/init_util.hh>
#include <core/kinematics/DomainMap.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>


//Auto Headers
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>

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
		for ( Size ii = 2; ii <= trpcage.size(); ++ii ) {

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
			TS_ASSERT_DELTA( emap[ cart_bonded ], 13.58958408336582, 1e-12 );
		}
		{
			EnergyMap emap;
			cartbond_energy.residue_pair_energy( trpcage.residue( 2 ), trpcage.residue( 3 ), trpcage, sfxn, emap );
			TS_ASSERT_DELTA( emap[ cart_bonded ], 1.678549902433499, 1e-12 );
		}
		{
			EnergyMap emap;
			cartbond_energy.residue_pair_energy( trpcage.residue( 3 ), trpcage.residue( 4 ), trpcage, sfxn, emap );
			TS_ASSERT_DELTA( emap[ cart_bonded ], 1.042759035471679, 1e-12 );
		}
		{
			EnergyMap emap;
			cartbond_energy.residue_pair_energy( trpcage.residue( 4 ), trpcage.residue( 5 ), trpcage, sfxn, emap );
			TS_ASSERT_DELTA( emap[ cart_bonded ], 0.5222409639108707, 1e-12 );
		}
		{
			EnergyMap emap;
			cartbond_energy.residue_pair_energy( trpcage.residue( 5 ), trpcage.residue( 6 ), trpcage, sfxn, emap );
			TS_ASSERT_DELTA( emap[ cart_bonded ], 3.592030906461361, 1e-12 );
		}
		{
			EnergyMap emap;
			cartbond_energy.residue_pair_energy( trpcage.residue( 6 ), trpcage.residue( 7 ), trpcage, sfxn, emap );
			TS_ASSERT_DELTA( emap[ cart_bonded ], 1.930873905641382, 1e-12 );
		}
		{
			EnergyMap emap;
			cartbond_energy.residue_pair_energy( trpcage.residue( 7 ), trpcage.residue( 8 ), trpcage, sfxn, emap );
			TS_ASSERT_DELTA( emap[ cart_bonded ], 1.652105725187637, 1e-12 );
		}
		{
			EnergyMap emap;
			cartbond_energy.residue_pair_energy( trpcage.residue( 8 ), trpcage.residue( 9 ), trpcage, sfxn, emap );
			TS_ASSERT_DELTA( emap[ cart_bonded ], 0.3838832036725081, 1e-12 );
		}
		{
			EnergyMap emap;
			cartbond_energy.residue_pair_energy( trpcage.residue( 9 ), trpcage.residue( 10 ), trpcage, sfxn, emap );
			TS_ASSERT_DELTA( emap[ cart_bonded ], 2.020811577596775, 1e-12 );
		}
		{
			EnergyMap emap;
			cartbond_energy.residue_pair_energy( trpcage.residue( 10 ), trpcage.residue( 11 ), trpcage, sfxn, emap );
			TS_ASSERT_DELTA( emap[ cart_bonded ], 0.1839301975150479, 1e-12 );
		}
		{
			EnergyMap emap;
			cartbond_energy.residue_pair_energy( trpcage.residue( 11 ), trpcage.residue( 12 ), trpcage, sfxn, emap );
			TS_ASSERT_DELTA( emap[ cart_bonded ], 3.981797963933229, 1e-12 );
		}
		{
			EnergyMap emap;
			cartbond_energy.residue_pair_energy( trpcage.residue( 12 ), trpcage.residue( 13 ), trpcage, sfxn, emap );
			TS_ASSERT_DELTA( emap[ cart_bonded ], 2.392589574554818, 1e-12 );
		}
		{
			EnergyMap emap;
			cartbond_energy.residue_pair_energy( trpcage.residue( 13 ), trpcage.residue( 14 ), trpcage, sfxn, emap );
			TS_ASSERT_DELTA( emap[ cart_bonded ], 1.854623831539971, 1e-12 );
		}
		{
			EnergyMap emap;
			cartbond_energy.residue_pair_energy( trpcage.residue( 14 ), trpcage.residue( 15 ), trpcage, sfxn, emap );
			TS_ASSERT_DELTA( emap[ cart_bonded ], 3.733038058903135, 1e-12 );
		}
		{
			EnergyMap emap;
			cartbond_energy.residue_pair_energy( trpcage.residue( 15 ), trpcage.residue( 16 ), trpcage, sfxn, emap );
			TS_ASSERT_DELTA( emap[ cart_bonded ], 0.236528075657657, 1e-12 );
		}
		{
			EnergyMap emap;
			cartbond_energy.residue_pair_energy( trpcage.residue( 16 ), trpcage.residue( 17 ), trpcage, sfxn, emap );
			TS_ASSERT_DELTA( emap[ cart_bonded ], 9.743150769702101, 1e-12 );
		}
		{
			EnergyMap emap;
			cartbond_energy.residue_pair_energy( trpcage.residue( 17 ), trpcage.residue( 18 ), trpcage, sfxn, emap );
			TS_ASSERT_DELTA( emap[ cart_bonded ], 6.114907368286072, 1e-12 );
		}
		{
			EnergyMap emap;
			cartbond_energy.residue_pair_energy( trpcage.residue( 18 ), trpcage.residue( 19 ), trpcage, sfxn, emap );
			TS_ASSERT_DELTA( emap[ cart_bonded ], 6.315488021858255, 1e-12 );
		}
		{
			EnergyMap emap;
			cartbond_energy.residue_pair_energy( trpcage.residue( 19 ), trpcage.residue( 20 ), trpcage, sfxn, emap );
			cartbond_energy.eval_intrares_energy( trpcage.residue( 20 ), trpcage, sfxn, emap );
			TS_ASSERT_DELTA( emap[ cart_bonded ], 4.431813105970174, 1e-12 );
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
		adv.validate_start_func_matches_start_score( 32.700353135831, false );
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
		core::chemical::ResidueTypeSetCOP rs;
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
		core::chemical::ResidueTypeCOPs base_residue_types = rs->base_residue_types();
		for (
				Size q = 1; q <= base_residue_types.size(); q++ )  {
			core::chemical::ResidueTypeCOP restype = base_residue_types[ q ];
			TR << restype->name() << std::endl;
			{
				bool prepro(true);
				core::scoring::methods::ResidueCartBondedParameters const & rcbp(
					ipd.parameters_for_restype(*(restype), prepro));
				rcbp.bb_N_index();
				rcbp.bb_CA_index();
				rcbp.bb_O_index();
				rcbp.bb_H_index();
			}

			{
				bool prepro(false);
				core::scoring::methods::ResidueCartBondedParameters const & rcbp(
					ipd.parameters_for_restype(*(restype), prepro));
				rcbp.bb_N_index();
				rcbp.bb_CA_index();
				rcbp.bb_O_index();
				rcbp.bb_H_index();
			}
		}
	}

};

class IdealParametersDatabaseTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init_with_additional_options("-extra_res_fa core/scoring/methods/HEM.fa.params");
	}

	void tearDown() {}

	// Test to make sure that the autogenerated torsion parameters are what we expect
	void test_torsion_generation() {
		using namespace core::scoring::methods;
		using namespace core::chemical;

		ResidueTypeSetCOP restypeset( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
		ResidueType const & restype( restypeset->name_map("HEM") );

		// Load the reference items
		utility::io::izstream input_file( "core/scoring/methods/HEM.fa.tors" );
		utility::vector1< std::pair< utility::vector1< core::Size >, utility::vector1< core::Real > > > refdata;
		while ( input_file ) {
			std::string res, a1, a2, a3, a4;
			core::Real mu, K, period;
			input_file >> res >> a1 >> a2 >> a3 >> a4 >> mu >> K >> period;
			if ( ! input_file ) break;
			TR << "READING " << res << " " << a1 << " " << a2 << " " << a3 << " " << a4 << std::endl;
			utility::vector1< core::Size > atoms;
			atoms.push_back( restype.atom_index( a1 ) );
			atoms.push_back( restype.atom_index( a2 ) );
			atoms.push_back( restype.atom_index( a3 ) );
			atoms.push_back( restype.atom_index( a4 ) );
			utility::vector1< core::Real > params;
			params.push_back( mu );
			params.push_back( K );
			params.push_back( period );
			refdata.push_back( std::make_pair( atoms, params ) );
		}

		EnergyMethodOptions options; //go with defaults
		core::Real cartbonded_len, cartbonded_ang, cartbonded_tors, cartbonded_proton , cartbonded_improper;
		options.get_cartesian_bonded_parameters( cartbonded_len, cartbonded_ang, cartbonded_tors, cartbonded_proton , cartbonded_improper );
		IdealParametersDatabase param_db(cartbonded_len, cartbonded_ang, cartbonded_tors, cartbonded_proton , cartbonded_improper);

		ResidueCartBondedParameters const & params( param_db.parameters_for_restype( restype, false ) );

		utility::vector1< ResidueCartBondedParameters::torsion_parameter > const & torsions( params.improper_parameters() );

		TS_ASSERT_EQUALS(torsions.size(), refdata.size() );

		for ( core::Size ii(1); ii <= torsions.size(); ++ii ) {
			ResidueCartBondedParameters::Size4 atoms( torsions[ii].first );
			CartBondedParametersCOP params( torsions[ii].second );

			bool torsion_found( false );
			for ( core::Size jj(1); jj <= refdata.size(); ++jj ) {
				utility::vector1< core::Size > const & refatoms( refdata[jj].first );
				if ( atoms[1] != refatoms[1] && atoms[1] != refatoms[4] ) continue;
				if ( atoms[2] != refatoms[2] && atoms[2] != refatoms[3] ) continue;
				if ( atoms[3] != refatoms[2] && atoms[3] != refatoms[3] ) continue;
				if ( atoms[4] != refatoms[1] && atoms[4] != refatoms[4] ) continue;
				torsion_found = true;
				utility::vector1< core::Real > refparams( refdata[jj].second );
				TS_ASSERT_DELTA( refparams[1], params->mu(0,0), 0.0001 );
				TS_ASSERT_DELTA( refparams[2], params->K(0,0), 0.0001 );
				TS_ASSERT_DELTA( refparams[3], params->period(), 0.0001 );
				// Funky syntax to delete an item by index
				// We delete so we don't double compare
				refdata.erase( refdata.begin() + ( jj - 1 ) );
				break;
			}
			TS_ASSERT( torsion_found );
			if ( ! torsion_found ) {
				TR << "Did not find in reference: "
					<< restype.atom_name(atoms[1]) << " - "
					<< restype.atom_name(atoms[2]) << " - "
					<< restype.atom_name(atoms[3]) << " - "
					<< restype.atom_name(atoms[4]) << std::endl;
			}
		}

		for ( core::Size tt(1); tt <= refdata.size(); ++tt ) {
			utility::vector1< core::Size > const & refatoms( refdata[tt].first );
			TR << "Was not generated: "
				<< restype.atom_name(refatoms[1]) << " - "
				<< restype.atom_name(refatoms[2]) << " - "
				<< restype.atom_name(refatoms[3]) << " - "
				<< restype.atom_name(refatoms[4]) << std::endl;
		}
	}

};

