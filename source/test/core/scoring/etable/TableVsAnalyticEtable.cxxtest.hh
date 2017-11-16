// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/etable/EtableEnergy.cxxtest.hh
/// @brief  Unit tests for the lennard-jones and EEF1 solvation model.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>
#include <test/util/deriv_funcs.hh>

#include <core/types.hh>

// Unit headers
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/etable/EtableEnergy.hh>

// Package headers
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/AtomType.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

// Project headers
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>

#include <test/UTracer.hh>

//Auto Headers
#include <utility/vector1.hh>


using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.scoring.etable.EtableEnergy.cxxtest");

// using declarations
using namespace std;
using namespace core;
using namespace scoring;

///////////////////////////////////////////////////////////////////////////
/// @name EtableEnergyTest
/// @brief: Test the functionality of the EtableEnergy class
///////////////////////////////////////////////////////////////////////////
class TableVsAnalyticEtableTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void tearDown() {}

	Real
	normalized_difference(
		Real val1,
		Real val2
	)
	{
		if ( val1 == 0.0 && val2 == 0.0 ) return 0.0;
		return std::abs( val1 - val2 ) / std::max( std::abs(val1), std::abs(val2));
	}

	void
	etable_numeric_deriv(
		core::scoring::etable::Etable const & e,
		core::conformation::Atom const & at1,
		core::conformation::Atom const & at2,
		core::Real & dlj_atrE,
		core::Real & dlj_repE,
		core::Real & dfa_solE
	)
	{
		Real const delta = 1e-10;
		Real ljatr1, ljrep1, fasol1, dummy;
		Real ljatr2, ljrep2, fasol2;
		core::conformation::Atom at2moved( at2 );
		at2moved.xyz( Vector( at2.xyz().x() - delta, 0.0, 0.0 ) );
		e.analytic_etable_evaluation( at1, at2moved, ljatr1, ljrep1, fasol1, dummy );
		at2moved.xyz( Vector( at2.xyz().x() + delta, 0.0, 0.0 ) );
		e.analytic_etable_evaluation( at1, at2moved, ljatr2, ljrep2, fasol2, dummy );
		dlj_atrE = (ljatr2-ljatr1)/(2*delta);
		dlj_repE = (ljrep2-ljrep1)/(2*delta);
		dfa_solE = (fasol2-fasol1)/(2*delta);
	}

	void test_etable_analytic_derivatives()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::etable;
		using namespace core::scoring::methods;

		Pose pose = create_trpcage_ideal_pose();
		EnergyMethodOptions options; // default is fine

		Etable const & etable( *( ScoringManager::get_instance()->etable( options.etable_type() ).lock()) );
		AnalyticEtableEnergy ana_lj_energy( *( ScoringManager::get_instance()->etable( options.etable_type() ).lock()), options );
		TableLookupEtableEnergy tab_lj_energy( *( ScoringManager::get_instance()->etable( options.etable_type() ).lock()), options );

		conformation::Atom at1, at2;
		at1.type(1); at2.type(2);
		at1.xyz( Vector(0,0,0) ); at2.xyz( Vector(0,0,1) );

		Real step = 1e-2;
		Real range = ana_lj_energy.atomic_interaction_cutoff();
		Size nsteps = Size( range / step ) + 1;
		//Real d2;

		Real offset = 0; //1e-5;
		for ( Size ii = 1; ii <= static_cast< Size >(etable.n_atomtypes()); ++ii ) {
			at1.type(ii);
			for ( Size jj = 1; jj <= static_cast< Size >(etable.n_atomtypes()); ++jj ) {
				at2.type(jj);
				//std::cout << "looking at " << (*etable.atom_set().lock())[ii].name() << " " << (*etable.atom_set())[jj].name() << std::endl;


				for ( Size kk = 2; kk <= nsteps; ++kk ) { // skip d = 0.01; this is etable.min_dis and messes up derivative calculations
					Real dljatr_ana, dljrep_ana, dfasol_ana, dummy;
					Real dljatr_num, dljrep_num, dfasol_num;
					at2.xyz( Vector(step*kk+offset,0,0) );
					etable.analytic_etable_derivatives( at1, at2, dljatr_ana, dljrep_ana, dfasol_ana, dummy );
					etable_numeric_deriv( etable, at1, at2, dljatr_num, dljrep_num, dfasol_num );

					Real tolerance = 4e-4;
					TS_ASSERT( normalized_difference( dljatr_ana, dljatr_num ) < tolerance || dljatr_ana - dljatr_num < tolerance );
					TS_ASSERT( normalized_difference( dljrep_ana, dljrep_num ) < tolerance || dljrep_ana - dljrep_num < tolerance );
					TS_ASSERT( normalized_difference( dfasol_ana, dfasol_num ) < tolerance || dfasol_ana - dfasol_num < tolerance );
					if ( normalized_difference( dljatr_ana, dljatr_num ) > tolerance && dljatr_ana - dljatr_num > tolerance ) {
						std::cout << (*etable.atom_set().lock())[ii].name() << " " << (*etable.atom_set().lock())[jj].name() << " " << step * kk <<
							"   Atr: " << dljatr_ana << " " << dljatr_num << " diff: " << dljatr_ana - dljatr_num << std::endl;
					}
					if ( normalized_difference( dljrep_ana, dljrep_num  ) > tolerance && dljrep_ana - dljrep_num > tolerance ) {
						std::cout << (*etable.atom_set().lock())[ii].name() << " " << (*etable.atom_set().lock())[jj].name() << " " << step * kk <<
							"   Rep: " << dljrep_ana << " " << dljrep_num << " diff: " << dljrep_ana - dljrep_num << std::endl;
						std::cout << "atr: " << dljatr_ana << " rep: " << dljrep_ana << " sol: " << dfasol_ana << std::endl;
					}
					if ( normalized_difference( dfasol_ana, dfasol_num  ) > tolerance && dfasol_ana - dfasol_num > tolerance ) {
						std::cout << (*etable.atom_set().lock())[ii].name() << " " << (*etable.atom_set().lock())[jj].name() << " " << step * kk <<
							"   Sol: " << dfasol_ana << " " << dfasol_num << " diff: " << dfasol_ana - dfasol_num << std::endl;
					}

				}
			}
		}
	}

	void test_etable_analytic_lk_individual()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::etable;
		using namespace core::scoring::methods;

		Pose pose = create_trpcage_ideal_pose();
		EnergyMethodOptions options; // default is fine

		Etable const & etable( *( ScoringManager::get_instance()->etable( options.etable_type() ).lock()) );

		conformation::Atom at1, at2;
		at1.type(1); at2.type(2);
		at1.xyz( Vector(0,0,0) ); at2.xyz( Vector(0,0,1) );

		Real step = 1e-2;
		Real range = etable.max_dis();
		Size nsteps = Size( range / step ) + 1;
		//Real d2;

		Real offset = 0; //1e-5;
		for ( Size ii = 1; ii <= static_cast< Size >(etable.n_atomtypes()); ++ii ) {
			at1.type(ii);
			for ( Size jj = 1; jj <= static_cast< Size >(etable.n_atomtypes()); ++jj ) {
				at2.type(jj);
				//std::cout << "looking at " << (*etable.atom_set().lock())[ii].name() << " " << (*etable.atom_set())[jj].name() << std::endl;


				for ( Size kk = 2; kk <= nsteps; ++kk ) { // skip d = 0.01; this is etable.min_dis and messes up derivative calculations
					Real dummy_ljatr, dummy_ljrep, fasol_regular, dummy_d2;
					Real fasol_individual_1, fasol_individual_2;

					at2.xyz( Vector(step*kk+offset,0,0) );

					etable.analytic_etable_evaluation( at1, at2, dummy_ljatr, dummy_ljrep, fasol_regular, dummy_d2 );
					etable.analytic_lk_energy( at1, at2, fasol_individual_1, fasol_individual_2 );

					/// There's surprising disagreement in some atom pairs (22 vs 22) as high as 5e-4
					/// hpark: even higher with beta_nov15, over 1e-3 b/w 44-50
					TS_ASSERT_DELTA( fasol_regular, fasol_individual_1 + fasol_individual_2, 1.1e-3/*1e-3*/  );
					if ( std::abs( fasol_regular - ( fasol_individual_1 + fasol_individual_2 )) > 1.1e-3/*1e-3*/  ) {
						std::cout << "Error: " << ii << " " << jj << " " << step*kk+offset << " " << fasol_regular << " vs " << fasol_individual_1 + fasol_individual_2 << " ( " << fasol_individual_1 << " + " << fasol_individual_2 << "); diff= " << fasol_regular - ( fasol_individual_1 + fasol_individual_2 ) << std::endl;
					}

				}
			}
		}
	}

	void dont_test_etable_analytic_evaluation()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::etable;
		using namespace core::scoring::methods;

		Pose pose = create_trpcage_ideal_pose();
		EnergyMethodOptions options; // default is fine

		Etable const & etable( *( ScoringManager::get_instance()->etable( options.etable_type() ).lock()) );
		TableLookupEtableEnergy tab_lj_energy( *( ScoringManager::get_instance()->etable( options.etable_type() ).lock()), options );
		AnalyticEtableEnergy    ana_lj_energy( *( ScoringManager::get_instance()->etable( options.etable_type() ).lock()), options );
		EnergyMap emap, emap2;

		conformation::Atom at1, at2;
		at1.type(1); at2.type(2);
		at1.xyz( Vector(0,0,0) ); at2.xyz( Vector(0,0,1) );

		Real step = 1e-2;
		Real range = tab_lj_energy.atomic_interaction_cutoff();
		Size nsteps = Size( range / step ) + 1;
		Real d2;
		//int const OCbb_idx = etable.atom_set()->atom_type_index("OCbb");
		int const Hha_idx = etable.atom_set().lock()->atom_type_index("Hha" );
		int const HREPS_idx = etable.atom_set().lock()->atom_type_index("HREPS" );

		Size count_failures = 0;
		ifstream infile( "save_etable_values.txt"  );
		for ( Size ii = 1; ii <= static_cast< Size >(etable.n_atomtypes()); ++ii ) {
			at1.type(ii);
			for ( Size jj = ii; jj <= static_cast< Size >(etable.n_atomtypes()); ++jj ) {
				at2.type(jj);
				//std::cout << "looking at " << (*etable.atom_set().lock())[ii].name() << " " << (*etable.atom_set())[jj].name() << std::endl;

				for ( Size kk = 1; kk <= nsteps; ++kk ) {
					bool failed = false;
					Real dummy;

					// 1. Make sure the table lookup still works.
					emap.zero();
					emap2.zero();
					at2.xyz( Vector(step*kk,0,0) );

					tab_lj_energy.atom_pair_energy( at1, at2, 1.0, emap, d2 );
					tab_lj_energy.atom_pair_energy( at2, at1, 1.0, emap2, d2 );
					TS_ASSERT_DELTA( emap[ fa_atr ], emap2[ fa_atr ], 1e-16 );
					TS_ASSERT_DELTA( emap[ fa_rep ], emap2[ fa_rep ], 1e-16 );
					TS_ASSERT_DELTA( emap[ fa_sol ], emap2[ fa_sol ], 1e-16 );

					//std::cout.precision(16);
					//std::cout << ii << " " << jj << " " << step*kk << " " << emap[fa_atr] << " " << emap[fa_rep] << " " << emap[ fa_sol ] << std::endl;
					Real inii, injj, instep, infaatr, infarep, infasol;
					infile >> inii >> injj >> instep >> infaatr >> infarep >> infasol;
					TS_ASSERT( inii == ii );
					TS_ASSERT( injj == jj );
					TS_ASSERT_DELTA( instep, step*kk, 1e-6 );
					TS_ASSERT_DELTA( emap[fa_atr], infaatr, 1e-3 );
					TS_ASSERT_DELTA( emap[fa_rep], infarep, 1e-3 );
					TS_ASSERT_DELTA( emap[fa_sol], infasol, 1e-3 );

					EnergyMap emap_ana;
					Real ana_vs_table_percent_diff_tolerance = ( ii==Size(Hha_idx) &&  ( jj==Size(Hha_idx) || jj == Size(HREPS_idx) ) ) ? 1 : ( d2 < 1 ) ? 1e-1 : 5e-2;
					ana_lj_energy.atom_pair_energy( at1, at2, 1.0, emap_ana, d2 );
					TS_ASSERT( ( emap[ fa_atr ] < 0.01 && emap_ana[ fa_atr ] < 0.01 ) ||
						std::abs( emap[ fa_atr ] - emap_ana[ fa_atr ] ) / std::max( std::abs( emap[ fa_atr ]), std::abs(emap_ana[ fa_atr ])) < ana_vs_table_percent_diff_tolerance );
					TS_ASSERT( ( emap[ fa_rep ] < 0.01 && emap_ana[ fa_rep ] < 0.01 ) ||
						std::abs( emap[ fa_rep ] - emap_ana[ fa_rep ] ) / std::max( std::abs( emap[ fa_rep ]), std::abs(emap_ana[ fa_rep ])) < ana_vs_table_percent_diff_tolerance );
					TS_ASSERT( ( emap[ fa_sol ] < 0.01 && emap_ana[ fa_sol ] < 0.01 ) ||
						std::abs( emap[ fa_sol ] - emap_ana[ fa_sol ] ) / std::max( std::abs( emap[ fa_sol ]), std::abs(emap_ana[ fa_sol ])) < ana_vs_table_percent_diff_tolerance );

					if ( ! (( emap[ fa_atr ] < 0.01 && emap_ana[ fa_atr ] < 0.01 ) ||
							( std::abs( emap[ fa_atr ] - emap_ana[ fa_atr ] ) / std::max( std::abs( emap[ fa_atr ]), std::abs(emap_ana[ fa_atr ])) < ana_vs_table_percent_diff_tolerance )) ) {
						std::cout <<  (*etable.atom_set().lock())[ii].name() << " " << (*etable.atom_set().lock())[jj].name() << " " << step * kk << " fa_atr " <<
							emap[ fa_atr ] << " " << emap_ana[ fa_atr ] << " " <<
							std::abs( emap[ fa_atr ] - emap_ana[ fa_atr ] ) / std::max( std::abs( emap[ fa_atr ]), std::abs(emap_ana[ fa_atr ])) << std::endl;
					}

					if ( ! ( ( emap[ fa_rep ] < 0.01 && emap_ana[ fa_rep ] < 0.01 ) ||
							std::abs( emap[ fa_rep ] - emap_ana[ fa_rep ] ) / std::max( std::abs( emap[ fa_rep ]), std::abs(emap_ana[ fa_rep ])) < ana_vs_table_percent_diff_tolerance ) ) {
						std::cout <<  (*etable.atom_set().lock())[ii].name() << " " << (*etable.atom_set().lock())[jj].name() << " " << step * kk << " fa_rep " <<
							emap[ fa_rep ] << " " << emap_ana[ fa_rep ] << " " <<
							std::abs( emap[ fa_rep ] - emap_ana[ fa_rep ] ) / std::max( std::abs( emap[ fa_rep ]), std::abs(emap_ana[ fa_rep ])) << std::endl;
					}
					if ( ! ( ( emap[ fa_sol ] < 0.01 && emap_ana[ fa_sol ] < 0.01 ) ||
							std::abs( emap[ fa_sol ] - emap_ana[ fa_sol ] ) / std::max( std::abs( emap[ fa_sol ]), std::abs(emap_ana[ fa_sol ])) < ana_vs_table_percent_diff_tolerance ) ) {
						std::cout <<  (*etable.atom_set().lock())[ii].name() << " " << (*etable.atom_set().lock())[jj].name() << " " << step * kk << " fa_sol " <<
							emap[ fa_sol ] << " " << emap_ana[ fa_sol ] << " " <<
							std::abs( emap[ fa_sol ] - emap_ana[ fa_sol ] ) / std::max( std::abs( emap[ fa_sol ]), std::abs(emap_ana[ fa_sol ])) << std::endl;
					}

					//if ( ii==OCbb_idx && jj==OCbb_idx && step*kk >= 2.1 && step*kk <=3.6 ) {
					// std::cout << "OCbb vs OCbb " << step*kk << " " << emap[fa_atr] << " " << emap[fa_rep] << " " << emap[ fa_sol ] << std::endl;
					//}
					Real an_ljatrE, an_ljrepE, an_fasolE;
					etable.interpolated_analytic_etable_evaluation( at1, at2, an_ljatrE, an_ljrepE, an_fasolE, dummy );
					if ( std::abs( an_ljatrE - emap[ fa_atr ]) > 1e-6 ) {
						failed = true;
						std::cout << (*etable.atom_set().lock())[ii].name() << " " << (*etable.atom_set().lock())[jj].name() << " " << step * kk <<
							"   Atr: " << an_ljatrE << " " << emap[ fa_atr ] << " diff: " << an_ljatrE - emap[ fa_atr ] << std::endl;
						//std::cout << (*etable.atom_set())[jj].name() << " " <<  (*etable.atom_set())[jj].name() << " e.ljatr_final_weight " << etable.ljatr_final_weight(ii,jj) <<  std::endl;
					}
					if ( std::abs( an_ljrepE - emap[ fa_rep ]) > 1e-6 ) {
						failed = true;
						std::cout << (*etable.atom_set().lock())[ii].name() << " " << (*etable.atom_set().lock())[jj].name() << " " << step * kk <<
							"   Rep: " << an_ljrepE << " " << emap[ fa_rep ] << " diff: " << an_ljrepE - emap[ fa_rep ] << std::endl;
					}
					if ( std::abs( an_fasolE - emap[ fa_sol ]) > 1e-4 ) {
						failed = true;
						std::cout << (*etable.atom_set().lock())[ii].name() << " " << (*etable.atom_set().lock())[jj].name() << " " << step * kk <<
							"   Sol: " << an_fasolE << " " << emap[ fa_sol ] << " diff: " << an_fasolE - emap[ fa_sol ] << std::endl;
						//std::cout << (*etable.atom_set().lock())[ii].name() << " " <<  (*etable.atom_set())[jj].name() << " e.fasol_final_weight " << etable.fasol_final_weight(ii,jj) << std::endl;
					}
					if ( failed ) ++count_failures;
				}
			}
		}
		std::cout << "nfailed: " << count_failures << std::endl;
	}

};
