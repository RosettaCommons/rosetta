// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
// AUTO-REMOVED #include <core/scoring/etable/BaseEtableEnergy.tmpl.hh>

// Package headers
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/AtomType.hh>
// AUTO-REMOVED #include <core/scoring/EnergyGraph.hh>
// AUTO-REMOVED #include <core/scoring/Energies.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/scoring/MinimizationData.hh>
// AUTO-REMOVED #include <core/scoring/MinimizationGraph.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

// Project headers
// AUTO-REMOVED #include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/kinematics/MoveMap.hh>
// AUTO-REMOVED #include <core/optimization/MinimizerMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
// AUTO-REMOVED #include <core/optimization/NumericalDerivCheckResult.hh>

// AUTO-REMOVED #include <basic/Tracer.hh>
#include <test/UTracer.hh>

//Auto Headers
#include <utility/vector1.hh>


using basic::T;
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

	/*void interpolated_analytic_etable_evaluation(
		core::scoring::etable::Etable const & e,
		core::conformation::Atom const & at1,
		core::conformation::Atom const & at2,
		core::Real & lj_atrE,
		core::Real & lj_repE,
		core::Real & fa_solE
	)
	{
		Real ljatrE_lo, ljrepE_lo, fasolE_lo;
		Real ljatrE_hi, ljrepE_hi, fasolE_hi;
		core::conformation::Atom at2_lo( at2 ), at2_hi( at2 );
		Real dis2 = at1.xyz().distance_squared( at2.xyz() );
		int dis2bin = (int) ( dis2 * e.bins_per_A2 );

		Real dis2lo = Real(dis2bin) / e.bins_per_A2;
		Real dis2hi = Real(dis2bin  + 1 )/e.bins_per_A2;
		at2_lo.xyz( Vector( std::sqrt( dis2lo ), 0, 0 ));
		at2_hi.xyz( Vector( std::sqrt( dis2hi ), 0, 0 ));

		analytic_etable_evaluation( e, at1, at2_lo, ljatrE_lo, ljrepE_lo, fasolE_lo );
		analytic_etable_evaluation( e, at1, at2_hi, ljatrE_hi, ljrepE_hi, fasolE_hi );

		Real alpha = (dis2*e.bins_per_A2 - dis2bin);
		//std::cout << "debug " << dis2 << " " << dis2lo << " " << dis2hi << " " << alpha << " " << 1-alpha << " " << ljatrE_lo << " " << ljatrE_hi << std::endl;

		lj_atrE = alpha * ljatrE_hi + (1-alpha) * ljatrE_lo;
		lj_repE = alpha * ljrepE_hi + (1-alpha) * ljrepE_lo;
		fa_solE = alpha * fasolE_hi + (1-alpha) * fasolE_lo;
	}

	void
	analytic_etable_evaluation(
		core::scoring::etable::Etable const & e,
		core::conformation::Atom const & at1,
		core::conformation::Atom const & at2,
		core::Real & lj_atrE,
		core::Real & lj_repE,
		core::Real & fa_solE
	)
	{
		using namespace core;
		using namespace core::scoring::etable;
		Real dis2 = at1.xyz().distance_squared( at2.xyz() );

		int atype1 = at1.type();
		int atype2 = at2.type();
		lj_atrE = 0; lj_repE = 0; fa_solE = 0;

		// Pasted directly from Etable.cc:calc_etable_value
		// locals
		Real ljE,d_ljE,x1,x2;
		Real dis;
		Real inv_dis,inv_dis2,inv_dis6,inv_dis7,inv_dis12,inv_dis13;
		Real dis2sigma;
		int xtra_atype1,xtra_atype2;
		//int const atype_sulfur = { 16 };

		// include after local variables to allow data statements to initialize
		Real atrE = 0.;
		Real d_atrE = 0.;
		Real repE = 0.;
		Real d_repE = 0.;
		Real solvE1 = 0.;
		Real solvE2 = 0.;
		Real dsolvE1 = 0.;
		Real dsolvE2 = 0.;

		//  ctsa - epsilon allows final bin value to be calculated
		if ( dis2 > e.max_dis2 + e.epsilon ) return;

		if ( dis2 < e.min_dis2 ) dis2 = e.min_dis2;

		dis = std::sqrt(dis2);
		inv_dis = 1.0/dis;
		inv_dis2 = inv_dis * inv_dis;



		//  ctsa - switch to disulfide bonded atom types
		//    when conditions are met
		// if ( ( atype1 == atype_sulfur && atype2 == atype_sulfur ) &&
		//  dis < disulfide_dis_thresh ) {
		// xtra_atype1 = n_atomtypes + 1;
		// xtra_atype2 = n_atomtypes + 1;
		// } else {
		xtra_atype1 = atype1;
		xtra_atype2 = atype2;
		//}


		dis2sigma = dis / e.lj_sigma(xtra_atype1,xtra_atype2);

		std::pair< Real, Real > const & ljatr_spline_xlo_xhi = e.ljatr_spline_xlo_xhi(atype1,atype2);
		if ( dis2sigma < e.lj_switch_dis2sigma ) {
			//  ctsa - use linear ramp instead of lj when the dis/sigma
			//    ratio drops below theshold
			d_ljE = e.lj_switch_slope(xtra_atype1,xtra_atype2);
			ljE = dis*d_ljE + e.lj_switch_intercept(xtra_atype1,xtra_atype2);
		} else if ( dis < ljatr_spline_xlo_xhi.first ) {
			//  ctsa - calc regular lennard-jones
			inv_dis6  = inv_dis2 * inv_dis2 * inv_dis2;
			inv_dis7  = inv_dis6 * inv_dis;
			inv_dis12 = inv_dis6 * inv_dis6;
			inv_dis13 = inv_dis12 * inv_dis;

			ljE = e.lj_r12_coeff(xtra_atype1,xtra_atype2) * inv_dis12 +
				e.lj_r6_coeff(xtra_atype1,xtra_atype2) * inv_dis6;

			d_ljE = -12.*e.lj_r12_coeff(xtra_atype1,xtra_atype2) * inv_dis13-6. *
				e.lj_r6_coeff(xtra_atype1,xtra_atype2) * inv_dis7;
		} else if ( dis < ljatr_spline_xlo_xhi.second ) {
			Real width = (ljatr_spline_xlo_xhi.second - ljatr_spline_xlo_xhi.first);
			Real invwidth = 1/width;
			Real a = ( ljatr_spline_xlo_xhi.second - dis ) * invwidth;
			Real b = ( dis - ljatr_spline_xlo_xhi.first ) * invwidth;
			SplineParameters sp = e.ljatr_spline_parameters( at1.type(), at2.type() );
			ljE = a*sp.ylo + b*sp.yhi + ((a*a*a-a)*sp.y2lo + (b*b*b-b)*sp.y2hi)*width*width / 6;
		} else {
			/// assuming e.ljatr_spline_xhi == LK distance cutoff
			return;
		}

		if ( e.ljrep_from_negcrossing(at1.type(), at2.type() )) {
			if (ljE < 0 ) {
				atrE = ljE;
			} else {
				repE = ljE;
			}
		} else if ( dis < e.lj_minima( at1.type(), at2.type() ) ) {
			atrE = e.lj_vals_at_minima( at1.type(), at2.type() );
			repE = ljE - atrE;
		} else {
			atrE = ljE;
		}

		ExtraQuadraticRepulsion const & exrep = e.ljrep_extra_repulsion(atype1,atype2);
		if ( dis < exrep.xhi ) {
			if ( dis < exrep.xlo ) {
				repE += ( dis - exrep.xlo ) * exrep.extrapolated_slope + exrep.ylo;
			} else {
				repE += ( exrep.xhi - dis ) * ( exrep.xhi - dis ) * exrep.slope;
			}
		}

		//if ( ljE < 0. ) {
		//	atrE = ljE;
		//	d_atrE = d_ljE;
		//} else {
		//	repE = ljE;
		//	d_repE = d_ljE;
		//}


		std::pair< Real, Real > fasol_spline_close_knots( e.fasol_spline_close_start_end(at1.type(),at2.type()) );
		if ( dis < fasol_spline_close_knots.first ) {
			fa_solE = e.fasol_spline_close(at1.type(),at2.type() ).ylo;
		} else if ( dis < fasol_spline_close_knots.second ) {
			Real fasol_spline_knots_diff = fasol_spline_close_knots.second - fasol_spline_close_knots.first;
			Real fasol_spline_knots_diff_inv = 1/fasol_spline_knots_diff;
      Real a = ( fasol_spline_close_knots.second - dis ) * fasol_spline_knots_diff_inv;
      Real b = ( dis - fasol_spline_close_knots.first ) * fasol_spline_knots_diff_inv;
      SplineParameters const & sp = e.fasol_spline_close( at1.type(), at2.type() );
      fa_solE = a*sp.ylo + b*sp.yhi + ((a*a*a-a)*sp.y2lo + (b*b*b-b)*sp.y2hi)*fasol_spline_knots_diff*fasol_spline_knots_diff / 6;
		} else if ( dis < e.fasol_spline_far_xlo ) {
			/// exponential evaluation
      Real dis_rad = dis - e.lj_radius(atype1);
      //Real dis_rad = dis - atom_type(atype1).lj_radius();
      x1 = ( dis_rad * dis_rad ) * e.lk_inv_lambda2(atype1);
      dis_rad = dis - e.lj_radius(atype2);
      //dis_rad = dis - atom_type(atype2).lj_radius();
      x2 = ( dis_rad * dis_rad ) * e.lk_inv_lambda2(atype2);

      fa_solE =  std::exp(-x1) * e.lk_coeff(atype1,atype2) * inv_dis2;
      fa_solE += std::exp(-x2) * e.lk_coeff(atype2,atype1) * inv_dis2;
		} else if ( dis < e.fasol_spline_far_xhi ) {
      Real a = ( e.fasol_spline_far_xhi - dis ) * e.fasol_spline_far_diff_xhi_xlo_inv;
      Real b = ( dis - e.fasol_spline_far_xlo ) * e.fasol_spline_far_diff_xhi_xlo_inv;
      SplineParameters const & sp = e.fasol_spline_far( at1.type(), at2.type() );
      fa_solE = a*sp.ylo + b*sp.yhi + ((a*a*a-a)*sp.y2lo + (b*b*b-b)*sp.y2hi)*e.fasol_spline_far_diff_xhi_xlo*e.fasol_spline_far_diff_xhi_xlo / 6;
		} else {
			fa_solE = 0;
		}

		lj_atrE = atrE * e.ljatr_final_weight(at1.type(),at2.type());
		lj_repE = repE;
		fa_solE *= e.fasol_final_weight(at1.type(),at2.type());

		//fa_solE = solvE1 + solvE2;
	}*/

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

		Etable const & etable( *( ScoringManager::get_instance()->etable( options.etable_type() )) );
		AnalyticEtableEnergy ana_lj_energy( *( ScoringManager::get_instance()->etable( options.etable_type() )), options );
		TableLookupEtableEnergy tab_lj_energy( *( ScoringManager::get_instance()->etable( options.etable_type() )), options );

		conformation::Atom at1, at2;
		at1.type(1); at2.type(2);
		at1.xyz( Vector(0,0,0) ); at2.xyz( Vector(0,0,1) );

		Real step = 1e-2;
		Real range = ana_lj_energy.atomic_interaction_cutoff();
		Size nsteps = Size( range / step ) + 1;
		//Real d2;

		Real offset = 0; //1e-5;
		for ( Size ii = 1; ii <= etable.n_atomtypes(); ++ii ) {
			at1.type(ii);
			for ( Size jj = 1; jj <= etable.n_atomtypes(); ++jj ) {
				at2.type(jj);
				//std::cout << "looking at " << (*etable.atom_set())[ii].name() << " " << (*etable.atom_set())[jj].name() << std::endl;


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
						std::cout << (*etable.atom_set())[ii].name() << " " << (*etable.atom_set())[jj].name() << " " << step * kk <<
							"   Atr: " << dljatr_ana << " " << dljatr_num << " diff: " << dljatr_ana - dljatr_num << std::endl;
					}
					if ( normalized_difference( dljrep_ana, dljrep_num  ) > tolerance && dljrep_ana - dljrep_num > tolerance ) {
						std::cout << (*etable.atom_set())[ii].name() << " " << (*etable.atom_set())[jj].name() << " " << step * kk <<
							"   Rep: " << dljrep_ana << " " << dljrep_num << " diff: " << dljrep_ana - dljrep_num << std::endl;
						std::cout << "atr: " << dljatr_ana << " rep: " << dljrep_ana << " sol: " << dfasol_ana << std::endl;
					}
					if ( normalized_difference( dfasol_ana, dfasol_num  ) > tolerance && dfasol_ana - dfasol_num > tolerance ) {
						std::cout << (*etable.atom_set())[ii].name() << " " << (*etable.atom_set())[jj].name() << " " << step * kk <<
							"   Sol: " << dfasol_ana << " " << dfasol_num << " diff: " << dfasol_ana - dfasol_num << std::endl;
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

		Etable const & etable( *( ScoringManager::get_instance()->etable( options.etable_type() )) );
		TableLookupEtableEnergy tab_lj_energy( *( ScoringManager::get_instance()->etable( options.etable_type() )), options );
		AnalyticEtableEnergy    ana_lj_energy( *( ScoringManager::get_instance()->etable( options.etable_type() )), options );
		EnergyMap emap, emap2;

		conformation::Atom at1, at2;
		at1.type(1); at2.type(2);
		at1.xyz( Vector(0,0,0) ); at2.xyz( Vector(0,0,1) );

		Real step = 1e-2;
		Real range = tab_lj_energy.atomic_interaction_cutoff();
		Size nsteps = Size( range / step ) + 1;
		Real d2;
		//int const OCbb_idx = etable.atom_set()->atom_type_index("OCbb");
		int const Hha_idx = etable.atom_set()->atom_type_index("Hha" );
		int const HREPS_idx = etable.atom_set()->atom_type_index("HREPS" );

		Size count_failures = 0;
		ifstream infile( "save_etable_values.txt"  );
		for ( Size ii = 1; ii <= etable.n_atomtypes(); ++ii ) {
			at1.type(ii);
			for ( Size jj = ii; jj <= etable.n_atomtypes(); ++jj ) {
				at2.type(jj);
				//std::cout << "looking at " << (*etable.atom_set())[ii].name() << " " << (*etable.atom_set())[jj].name() << std::endl;

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
					Real ana_vs_table_percent_diff_tolerance = ( ii==Hha_idx &&  ( jj==Hha_idx || jj == HREPS_idx ) ) ? 1 : ( d2 < 1 ) ? 1e-1 : 5e-2;
					ana_lj_energy.atom_pair_energy( at1, at2, 1.0, emap_ana, d2 );
					TS_ASSERT( ( emap[ fa_atr ] < 0.01 && emap_ana[ fa_atr ] < 0.01 ) ||
						std::abs( emap[ fa_atr ] - emap_ana[ fa_atr ] ) / std::max( std::abs( emap[ fa_atr ]), std::abs(emap_ana[ fa_atr ])) < ana_vs_table_percent_diff_tolerance );
					TS_ASSERT( ( emap[ fa_rep ] < 0.01 && emap_ana[ fa_rep ] < 0.01 ) ||
						std::abs( emap[ fa_rep ] - emap_ana[ fa_rep ] ) / std::max( std::abs( emap[ fa_rep ]), std::abs(emap_ana[ fa_rep ])) < ana_vs_table_percent_diff_tolerance );
					TS_ASSERT( ( emap[ fa_sol ] < 0.01 && emap_ana[ fa_sol ] < 0.01 ) ||
						std::abs( emap[ fa_sol ] - emap_ana[ fa_sol ] ) / std::max( std::abs( emap[ fa_sol ]), std::abs(emap_ana[ fa_sol ])) < ana_vs_table_percent_diff_tolerance );

					if ( ! (( emap[ fa_atr ] < 0.01 && emap_ana[ fa_atr ] < 0.01 ) ||
							( std::abs( emap[ fa_atr ] - emap_ana[ fa_atr ] ) / std::max( std::abs( emap[ fa_atr ]), std::abs(emap_ana[ fa_atr ])) < ana_vs_table_percent_diff_tolerance )) ) {
						std::cout <<  (*etable.atom_set())[ii].name() << " " << (*etable.atom_set())[jj].name() << " " << step * kk << " fa_atr " <<
							emap[ fa_atr ] << " " << emap_ana[ fa_atr ] << " " <<
							std::abs( emap[ fa_atr ] - emap_ana[ fa_atr ] ) / std::max( std::abs( emap[ fa_atr ]), std::abs(emap_ana[ fa_atr ])) << std::endl;
					}

					if ( ! ( ( emap[ fa_rep ] < 0.01 && emap_ana[ fa_rep ] < 0.01 ) ||
							std::abs( emap[ fa_rep ] - emap_ana[ fa_rep ] ) / std::max( std::abs( emap[ fa_rep ]), std::abs(emap_ana[ fa_rep ])) < ana_vs_table_percent_diff_tolerance )) {
						std::cout <<  (*etable.atom_set())[ii].name() << " " << (*etable.atom_set())[jj].name() << " " << step * kk << " fa_rep " <<
							emap[ fa_rep ] << " " << emap_ana[ fa_rep ] << " " <<
							std::abs( emap[ fa_rep ] - emap_ana[ fa_rep ] ) / std::max( std::abs( emap[ fa_rep ]), std::abs(emap_ana[ fa_rep ])) << std::endl;
          }
					if ( ! ( ( emap[ fa_sol ] < 0.01 && emap_ana[ fa_sol ] < 0.01 ) ||
							std::abs( emap[ fa_sol ] - emap_ana[ fa_sol ] ) / std::max( std::abs( emap[ fa_sol ]), std::abs(emap_ana[ fa_sol ])) < ana_vs_table_percent_diff_tolerance ) ) {
						std::cout <<  (*etable.atom_set())[ii].name() << " " << (*etable.atom_set())[jj].name() << " " << step * kk << " fa_sol " <<
							emap[ fa_sol ] << " " << emap_ana[ fa_sol ] << " " <<
							std::abs( emap[ fa_sol ] - emap_ana[ fa_sol ] ) / std::max( std::abs( emap[ fa_sol ]), std::abs(emap_ana[ fa_sol ])) << std::endl;
					}

					//if ( ii==OCbb_idx && jj==OCbb_idx && step*kk >= 2.1 && step*kk <=3.6 ) {
					//	std::cout << "OCbb vs OCbb " << step*kk << " " << emap[fa_atr] << " " << emap[fa_rep] << " " << emap[ fa_sol ] << std::endl;
					//}
					Real an_ljatrE, an_ljrepE, an_fasolE;
					etable.interpolated_analytic_etable_evaluation( at1, at2, an_ljatrE, an_ljrepE, an_fasolE, dummy );
					if ( std::abs( an_ljatrE - emap[ fa_atr ]) > 1e-6 ) {
						failed = true;
						std::cout << (*etable.atom_set())[ii].name() << " " << (*etable.atom_set())[jj].name() << " " << step * kk <<
							"   Atr: " << an_ljatrE << " " << emap[ fa_atr ] << " diff: " << an_ljatrE - emap[ fa_atr ] << std::endl;
						//std::cout << (*etable.atom_set())[jj].name() << " " <<  (*etable.atom_set())[jj].name() << " e.ljatr_final_weight " << etable.ljatr_final_weight(ii,jj) <<  std::endl;
					}
					if ( std::abs( an_ljrepE - emap[ fa_rep ]) > 1e-6 ) {
						failed = true;
						std::cout << (*etable.atom_set())[ii].name() << " " << (*etable.atom_set())[jj].name() << " " << step * kk <<
							"   Rep: " << an_ljrepE << " " << emap[ fa_rep ] << " diff: " << an_ljrepE - emap[ fa_rep ] << std::endl;
					}
					if ( std::abs( an_fasolE - emap[ fa_sol ]) > 1e-4 ) {
						failed = true;
						std::cout << (*etable.atom_set())[ii].name() << " " << (*etable.atom_set())[jj].name() << " " << step * kk <<
							"   Sol: " << an_fasolE << " " << emap[ fa_sol ] << " diff: " << an_fasolE - emap[ fa_sol ] << std::endl;
						//std::cout << (*etable.atom_set())[ii].name() << " " <<  (*etable.atom_set())[jj].name() << " e.fasol_final_weight " << etable.fasol_final_weight(ii,jj) << std::endl;
					}
					if ( failed ) ++count_failures;
				}
			}
		}
		std::cout << "nfailed: " << count_failures << std::endl;
	}

};
