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
#include <core/scoring/EnergyGraph.hh>

// Project headers
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>

// Utility Headers
#include <utility/vector1.hh>
//#include <utility/alignedfixedsizearray1.hh>
#include <utility/fixedsizearray1.hh>
#include <utility/FixedSizeLexicographicalIterator.hh>
#include <utility/FixedSizeLexicographicalIterator.tmpl.hh>

// Basic headers
#include <basic/Tracer.hh>

#include <ctime>

static basic::Tracer TR("core.scoring.etable.VectorizedEtable.cxxtest");

// using declarations
using namespace std;
using namespace core;
using namespace core::scoring;
using namespace core::scoring::etable;
using namespace core::scoring::methods;
using namespace core::pose;

typedef float VReal;

class EtableScratchSpace {
public:
	//typedef utility::vector0< Real > RealArray;
	//typedef utility::vector0< int > IntArray;


	typedef VReal * VRealArray;
	typedef int * IntArray;

	static int const N = 64;

	//typedef VReal[ N ] __attribute__((aligned(64))) VRealArray;
	//typedef int[ N ] __attribute__((aligned(64))) IntArray;

public:
	VRealArray x1, y1, z1;
	VRealArray x2, y2, z2;
	VRealArray lj_r6_coeff, lj_r12_coeff, ljrep_linramp_d2_cut, ljrep_linramp_slope, ljrep_linramp_intercept;
	VRealArray lj_exrep_xlo, lj_exrep_xhi, lj_exrep_slope, lj_exrep_extrap_slope, lj_exrep_ylo;
	VRealArray ljatr_final_weight, lj_minimum_x, lj_val_at_minimum;

	VRealArray ljatr_poly3_xhi, ljatr_poly3_xlo;
	VRealArray ljatr_poly3_c0, ljatr_poly3_c1, ljatr_poly3_c2, ljatr_poly3_c3;

	VRealArray lk_poly3_close_xlo, lk_poly3_close_xhi;
	VRealArray lk_poly3_close_flat; // value of the LK potential for distances below lk_poly3_close_xlo
	VRealArray lk_poly3_close_c0, lk_poly3_close_c1, lk_poly3_close_c2, lk_poly3_close_c3;

	VRealArray lj_radius_1, lj_radius_2, lk_inv_lambda2_1, lk_inv_lambda2_2, lk_coeff1, lk_coeff2;
	VRealArray lk_final_weight;

	VRealArray lk_poly3_far_xlo, lk_poly3_far_xhi;
	VRealArray lk_poly3_far_c0, lk_poly3_far_c1, lk_poly3_far_c2, lk_poly3_far_c3;

	// values that are computed along the way
	VRealArray d2, d, invd2, invd6;
	VRealArray lj_regular, ljrep_linramp, ljatr_ramp_to_zero, lj_extra_rep, ljE;
	VRealArray ljatr, ljrep;
	VRealArray lk_disrad1, lk_disrad2, lk_x1, lk_x2;
	VRealArray lk_close_poly3, lk_sol_regular, lk_far_poly3, lksol;

	// "boolean" arrays used to mask certain computations that should be discarded
	IntArray ljrep_from_negcrossing, lj_closer_than_minima, lj_below_linramp_start, lj_in_mid_region, lj_ramp_to_zero;
	IntArray lj_dbin; // 0 = beyond famaxdis, 1 = linear repulsion, 2 = regular, 3 = taper to zero

	IntArray lk_beneath_close_poly3_xlo, lk_in_close_poly3_range, lk_in_mid_range, lk_in_far_poly3_range;
	IntArray lk_dbin; // 0 = beyond famaxdis, 1 = beneath close xlo, 2 = in close sline range, 3 = mid, 4 = in far cubic_polynomial regionx


	EtableScratchSpace() :
		x1( new VReal[ N ] ),
		y1( new VReal[ N ] ),
		z1( new VReal[ N ] ),
		x2( new VReal[ N ] ),
		y2( new VReal[ N ] ),
		z2( new VReal[ N ] ),
		lj_r6_coeff( new VReal[ N ] ),
		lj_r12_coeff( new VReal[ N ] ),
		ljrep_linramp_d2_cut( new VReal[ N ] ),
		ljrep_linramp_slope( new VReal[ N ] ),
		ljrep_linramp_intercept( new VReal[ N ] ),
		lj_exrep_xlo( new VReal[ N ] ),
		lj_exrep_xhi( new VReal[ N ] ),
		lj_exrep_slope( new VReal[ N ] ),
		lj_exrep_extrap_slope( new VReal[ N ] ),
		lj_exrep_ylo( new VReal[ N ] ),
		ljatr_final_weight( new VReal[ N ] ),
		lj_minimum_x( new VReal[ N ] ),
		lj_val_at_minimum( new VReal[ N ] ),
		ljatr_poly3_xhi( new VReal[ N ] ),
		ljatr_poly3_xlo( new VReal[ N ] ),
		ljatr_poly3_c0( new VReal[ N ] ),
		ljatr_poly3_c1( new VReal[ N ] ),
		ljatr_poly3_c2( new VReal[ N ] ),
		ljatr_poly3_c3( new VReal[ N ] ),

		lk_poly3_close_xlo( new VReal[ N ] ),
		lk_poly3_close_xhi( new VReal[ N ] ),
		lk_poly3_close_flat( new VReal[ N ] ),

		lk_poly3_close_c0( new VReal[ N ] ),
		lk_poly3_close_c1( new VReal[ N ] ),
		lk_poly3_close_c2( new VReal[ N ] ),
		lk_poly3_close_c3( new VReal[ N ] ),

		lj_radius_1( new VReal[ N ] ),
		lj_radius_2( new VReal[ N ] ),
		lk_inv_lambda2_1( new VReal[ N ] ),
		lk_inv_lambda2_2( new VReal[ N ] ),
		lk_coeff1( new VReal[ N ] ),
		lk_coeff2( new VReal[ N ] ),
		lk_final_weight( new VReal[ N ] ),

		lk_poly3_far_xlo( new VReal[ N ] ),
		lk_poly3_far_xhi( new VReal[ N ] ),
		lk_poly3_far_c0( new VReal[ N ] ),
		lk_poly3_far_c1( new VReal[ N ] ),
		lk_poly3_far_c2( new VReal[ N ] ),
		lk_poly3_far_c3( new VReal[ N ] ),

		d2( new VReal[ N ] ),
		d( new VReal[ N ] ),
		invd2( new VReal[ N ] ),
		invd6( new VReal[ N ] ),
		lj_regular( new VReal[ N ] ),
		ljrep_linramp( new VReal[ N ] ),
		ljatr_ramp_to_zero( new VReal[ N ] ),
		lj_extra_rep( new VReal[ N ] ),
		ljE( new VReal[ N ] ),
		ljatr( new VReal[ N ] ),
		ljrep( new VReal[ N ] ),
		lk_disrad1( new VReal[ N ] ),
		lk_disrad2( new VReal[ N ] ),
		lk_x1( new VReal[ N ] ),
		lk_x2( new VReal[ N ] ),
		lk_close_poly3( new VReal[ N ] ),
		lk_sol_regular( new VReal[ N ] ),
		lk_far_poly3( new VReal[ N ] ),
		lksol( new VReal[ N ] ),
		ljrep_from_negcrossing( new int[ N ] ),
		lj_closer_than_minima( new int[ N ] ),
		lj_below_linramp_start( new int[ N ] ),
		lj_in_mid_region( new int[ N ] ),
		lj_ramp_to_zero( new int[ N ] ),
		lj_dbin( new int[ N ] ),

		lk_beneath_close_poly3_xlo( new int[ N ] ),
		lk_in_close_poly3_range( new int[ N ] ),
		lk_in_mid_range( new int[ N ] ),
		lk_in_far_poly3_range( new int[ N ]),
		lk_dbin( new int[ N ] )
	{}


	/*EtableScratchSpace() :
		x1( N ),
		y1( N ),
		z1( N ),
		x2( N ),
		y2( N ),
		z2( N ),
		lj_r6_coeff( N ),
		lj_r12_coeff( N ),
		ljrep_linramp_d2_cut( N ),
		ljrep_linramp_slope( N ),
		ljrep_linramp_intercept( N ),
		lj_exrep_xlo( N ),
		lj_exrep_xhi( N ),
		lj_exrep_slope( N ),
		lj_exrep_extrap_slope( N ),
		lj_exrep_ylo( N ),
		ljatr_final_weight( N ),
		lj_minimum_x( N ),
		lj_val_at_minimum( N ),
		ljatr_poly3_xhi( N ),
		ljatr_poly3_xlo( N ),
		ljatr_poly3_width( N ),
		ljatr_poly3_invwidth( N ),
		ljatr_poly3_ylo( N ),
		ljatr_poly3_yhi( N ),
		ljatr_poly3_y2lo( N ),
		ljatr_poly3_y2hi( N ),
		a( N ),
		b( N ),
		lk_poly3_close_xlo( N ),lk_poly3_close_xhi( N ),lk_poly3_close_width( N ),lk_poly3_close_invwidth( N ),
		lk_poly3_close_ylo( N ),lk_poly3_close_yhi( N ),lk_poly3_close_y2lo( N ),lk_poly3_close_y2hi( N ),

		lj_radius_1( N ),lj_radius_2( N ),lk_inv_lambda2_1( N ),lk_inv_lambda2_2( N ),lk_coeff1( N ),lk_coeff2( N ),
		lk_final_weight( N ),

		lk_poly3_far_xlo( N ),lk_poly3_far_xhi( N ),lk_poly3_far_width( N ),lk_poly3_far_invwidth( N ),
		lk_poly3_far_ylo( N ),lk_poly3_far_yhi( N ),lk_poly3_far_y2lo( N ),lk_poly3_far_y2hi( N ),

		d2( N ),
		d( N ),
		invd2( N ),
		invd6( N ),
		lj_regular( N ),
		ljrep_linramp( N ),
		ljatr_ramp_to_zero( N ),
		lj_extra_rep( N ),
		ljE( N ),
		ljatr( N ),
		ljrep( N ),
		lk_disrad1( N ),lk_disrad2( N ),lk_x1( N ),lk_x2( N ),
		lk_close_poly3( N ),lk_sol_regular( N ),lk_far_poly3( N ),lksol( N ),
		ljrep_from_negcrossing( N ),
		lj_closer_than_minima( N ),
		lj_below_linramp_start( N ),
		lj_in_mid_region( N ),
		lj_ramp_to_zero( N ),
		lk_beneath_close_poly3_xlo( N ),lk_in_close_poly3_range( N ),lk_in_mid_range( N ),lk_in_far_poly3_range(N)
		{}*/

	/*void resize( core::Size newsize ) {
		x1.resize( newsize );
		y1.resize( newsize );
		z1.resize( newsize );
		x2.resize( newsize );
		y2.resize( newsize );
		z2.resize( newsize );
		lj_r6_coeff.resize( newsize );
		lj_r12_coeff.resize( newsize );
		ljrep_linramp_d2_cut.resize( newsize );
		ljrep_linramp_slope.resize( newsize );
		ljrep_linramp_intercept.resize( newsize );
		lj_exrep_xlo.resize( newsize );
		lj_exrep_xhi.resize( newsize );
		lj_exrep_slope.resize( newsize );
		lj_exrep_extrap_slope.resize( newsize );
		lj_exrep_ylo.resize( newsize );
		ljatr_final_weight.resize( newsize );
		lj_minimum_x.resize( newsize );
		lj_val_at_minimum.resize( newsize );
		ljatr_poly3_xhi.resize( newsize );
		ljatr_poly3_xlo.resize( newsize );
		ljatr_poly3_width.resize( newsize );
		ljatr_poly3_invwidth.resize( newsize );
		ljatr_poly3_ylo.resize( newsize );
		ljatr_poly3_yhi.resize( newsize );
		ljatr_poly3_y2lo.resize( newsize );
		ljatr_poly3_y2hi.resize( newsize );
		a.resize( newsize );
		b.resize( newsize );
		lk_poly3_close_xlo.resize( newsize );
		lk_poly3_close_xhi.resize( newsize );
		lk_poly3_close_width.resize( newsize );
		lk_poly3_close_invwidth.resize( newsize );
		lk_poly3_close_ylo.resize( newsize );
		lk_poly3_close_yhi.resize( newsize );
		lk_poly3_close_y2lo.resize( newsize );
		lk_poly3_close_y2hi.resize( newsize );

		lj_radius_1.resize( newsize );
		lj_radius_2.resize( newsize );
		lk_inv_lambda2_1.resize( newsize );
		lk_inv_lambda2_2.resize( newsize );
		lk_coeff1.resize( newsize );
		lk_coeff2.resize( newsize );
		lk_final_weight.resize( newsize );

		lk_poly3_far_xlo.resize( newsize );
		lk_poly3_far_xhi.resize( newsize );
		lk_poly3_far_width.resize( newsize );
		lk_poly3_far_invwidth.resize( newsize );
		lk_poly3_far_ylo.resize( newsize );
		lk_poly3_far_yhi.resize( newsize );
		lk_poly3_far_y2lo.resize( newsize );
		lk_poly3_far_y2hi.resize( newsize );

		d2.resize( newsize );
		d.resize( newsize );
		invd2.resize( newsize );
		invd6.resize( newsize );
		lj_regular.resize( newsize );
		ljrep_linramp.resize( newsize );
		ljatr_ramp_to_zero.resize( newsize );
		lj_extra_rep.resize( newsize );
		ljE.resize( newsize );
		ljatr.resize( newsize );
		ljrep.resize( newsize );
		lk_disrad1.resize( newsize );
		lk_disrad2.resize( newsize );
		lk_x1.resize( newsize );
		lk_x2.resize( newsize );
		lk_close_poly3.resize( newsize );
		lk_sol_regular.resize( newsize );
		lk_far_poly3.resize( newsize );
		lksol.resize( newsize );
		ljrep_from_negcrossing.resize( newsize );
		lj_closer_than_minima.resize( newsize );
		lj_below_linramp_start.resize( newsize );
		lj_in_mid_region.resize( newsize );
		lj_ramp_to_zero.resize( newsize );
		lk_beneath_close_poly3_xlo.resize( newsize );
		lk_in_close_poly3_range.resize( newsize );
		lk_in_mid_range.resize( newsize );
		lk_in_far_poly3_range.resize( newsize );
	}*/
};

///////////////////////////////////////////////////////////////////////////
/// @name VectorizedEtableTests
/// @brief: Test the functionality of the EtableEnergy class
///////////////////////////////////////////////////////////////////////////
class VectorizedEtableTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void tearDown() {}

	VReal
	normalized_difference(
		VReal val1,
		VReal val2
	)
	{
		if ( val1 == 0.0 && val2 == 0.0 ) return 0.0;
		return std::abs( val1 - val2 ) / std::max( std::abs(val1), std::abs(val2));
	}

	void
	load_atom_pair_batch(
		core::scoring::etable::Etable const & etable,
		core::conformation::Residue const & r1,
		core::conformation::Residue const & r2,
		utility::FixedSizeLexicographicalIterator< 2 > & lex,
		int start_index,
		int jjlimit,
		EtableScratchSpace & s
	)
	{
		//for ( core::Size jj = 1; jj <= jjlimit; ++jj ) {
		for ( int jj = 0; jj < jjlimit; ++jj ) {
			//lex.set_position_from_index( jj + start_index - 1 );
			lex.set_position_from_index( jj + start_index );
			Vector const & a1coord = r1.xyz(lex[1]);
			Vector const & a2coord = r2.xyz(lex[2]);

			s.x1[ jj ] = a1coord.x();
			s.y1[ jj ] = a1coord.y();
			s.z1[ jj ] = a1coord.z();

			s.x2[ jj ] = a2coord.x();
			s.y2[ jj ] = a2coord.y();
			s.z2[ jj ] = a2coord.z();

			core::Size const a1 = r1.atom(lex[1]).type();
			core::Size const a2 = r2.atom(lex[2]).type();
			core::scoring::etable::EtableParamsOnePair const & params12 = etable.analytic_params_for_pair( a1, a2 );

			s.lj_r6_coeff[  jj ]            = params12.lj_r6_coeff;
			s.lj_r12_coeff[ jj ]            = params12.lj_r12_coeff;
			s.ljrep_linramp_d2_cut[ jj ]    = params12.ljrep_linear_ramp_d2_cutoff;
			s.ljrep_linramp_slope[ jj ]     = params12.lj_switch_slope;
			s.ljrep_linramp_intercept[ jj ] = params12.lj_switch_intercept;

			s.lj_exrep_xlo[ jj ]          = params12.ljrep_extra_repulsion.xlo;
			s.lj_exrep_xhi[ jj ]          = params12.ljrep_extra_repulsion.xhi;
			s.lj_exrep_slope[ jj ]        = params12.ljrep_extra_repulsion.slope;
			s.lj_exrep_extrap_slope[ jj ] = params12.ljrep_extra_repulsion.extrapolated_slope;
			s.lj_exrep_ylo[ jj ]          = params12.ljrep_extra_repulsion.ylo;

			s.ljatr_final_weight[ jj ] = params12.ljatr_final_weight;
			s.lj_minimum_x[ jj ]       = params12.lj_minimum;
			s.lj_val_at_minimum[ jj ]  = params12.lj_val_at_minimum;

			s.ljrep_from_negcrossing[ jj ] = params12.ljrep_from_negcrossing ? 1:0;

			s.ljatr_poly3_xlo[ jj ] = params12.ljatr_cubic_poly_xlo;
			s.ljatr_poly3_xhi[ jj ] = params12.ljatr_cubic_poly_xhi;
			s.ljatr_poly3_c0[ jj ]  = params12.ljatr_cubic_poly_parameters.c0;
			s.ljatr_poly3_c1[ jj ]  = params12.ljatr_cubic_poly_parameters.c1;
			s.ljatr_poly3_c2[ jj ]  = params12.ljatr_cubic_poly_parameters.c2;
			s.ljatr_poly3_c3[ jj ]  = params12.ljatr_cubic_poly_parameters.c3;

			s.lk_poly3_close_xlo[ jj ]  = params12.fasol_cubic_poly_close_start;
			s.lk_poly3_close_xhi[ jj ]  = params12.fasol_cubic_poly_close_end;
			s.lk_poly3_close_flat[ jj ] = params12.fasol_cubic_poly_close_flat;

			s.lk_poly3_close_c0[ jj ]  = params12.fasol_cubic_poly_close.c0;
			s.lk_poly3_close_c1[ jj ]  = params12.fasol_cubic_poly_close.c1;
			s.lk_poly3_close_c2[ jj ]  = params12.fasol_cubic_poly_close.c2;
			s.lk_poly3_close_c3[ jj ]  = params12.fasol_cubic_poly_close.c3;

			s.lk_poly3_far_xlo[ jj ] = etable.fasol_cubic_poly_far_xlo();
			s.lk_poly3_far_xhi[ jj ] = etable.fasol_cubic_poly_far_xhi();
			s.lk_poly3_far_c0[ jj ]  = params12.fasol_cubic_poly_far.c0;
			s.lk_poly3_far_c1[ jj ]  = params12.fasol_cubic_poly_far.c1;
			s.lk_poly3_far_c2[ jj ]  = params12.fasol_cubic_poly_far.c2;
			s.lk_poly3_far_c3[ jj ]  = params12.fasol_cubic_poly_far.c3;

			// the rule is: params12 is arranged such that atom-type-1 is smaller than or equal to atom-type-2.
			// and to calculate lk correctly, some terms need to be correcty paired with each other
			s.lj_radius_1[ jj ] = etable.lj_radius( a1 <= a2 ? a1 : a2 );
			s.lj_radius_2[ jj ] = etable.lj_radius( a1 <= a2 ? a2 : a1 );
			s.lk_inv_lambda2_1[ jj ] = etable.lk_inv_lambda2( a1 <= a2 ? a1 : a2 );
			s.lk_inv_lambda2_2[ jj ] = etable.lk_inv_lambda2( a1 <= a2 ? a2 : a1 );
			s.lk_coeff1[ jj ] = params12.lk_coeff1;
			s.lk_coeff2[ jj ] = params12.lk_coeff2;
			s.lk_final_weight[ jj ] = params12.fasol_final_weight;
		}

	}

	void width(
		VReal * __restrict__ width,
		VReal * __restrict__ xhi,
		VReal * __restrict__ xlo,
		int limit
	)
	{
		for ( int ii = 0; ii < limit; ++ii ) {
			width[ ii ] = xhi[ ii ] - xlo[ ii ];
		}
	}

	void calc_inv(
		VReal * __restrict__ invx,
		VReal * __restrict__ x,
		int limit
	)
	{
		for ( int ii = 0; ii < limit; ++ii ) {
			invx[ ii ] = 1.0 / x[ ii ];
		}
	}

	void
	calc_d2(
		VReal * __restrict__ d2,
		VReal * __restrict__ x1,
		VReal * __restrict__ x2,
		VReal * __restrict__ y1,
		VReal * __restrict__ y2,
		VReal * __restrict__ z1,
		VReal * __restrict__ z2,
		int limit
	)
	{
		// VECTORIZES! 2/27/2014 3:15PM
		for ( int ii = 0; ii < limit; ++ii ) {
			d2[ ii ] = ( x1[ii]-x2[ii] ) * ( x1[ii]-x2[ii] ) +
				( y1[ii]-y2[ii] ) * ( y1[ii]-y2[ii] ) +
				( z1[ii]-z2[ii] ) * ( z1[ii]-z2[ii] );
		}
	}

	void
	calc_d(
		VReal * __restrict__ d,
		VReal * __restrict__ d2,
		int limit
	)
	{
		for ( int ii = 0; ii < limit; ++ii ) {
			d[ ii ] = std::sqrt( d2[ ii ] );
		}
	}

	void
	calc_invd6(
		VReal * __restrict__ invd6,
		VReal * __restrict__ invd2,
		int limit
	)
	{
		for ( int ii = 0; ii < limit; ++ii ) {
			invd6[ ii ] = invd2[ii]*invd2[ii]*invd2[ii];
		}
	}

	void
	calc_lk_dbin(
		int   * __restrict__ lk_dbin,
		VReal * __restrict__ d,
		VReal * __restrict__ lk_poly3_close_xlo,
		VReal * __restrict__ lk_poly3_close_xhi,
		VReal * __restrict__ lk_poly3_far_xlo,
		VReal * __restrict__ lk_poly3_far_xhi,
		int limit
	)
	{
		for ( int ii = 0; ii < limit; ++ii ) {
			lk_dbin[ ii ] = (
				( d[ ii ] > lk_poly3_close_xlo[ ii ] ? 1 : 0 ) +
				( d[ ii ] > lk_poly3_close_xhi[ ii ] ? 1 : 0 ) +
				( d[ ii ] > lk_poly3_far_xlo[ ii ] ? 1 : 0 ) +
				1 ) *
				( d[ ii ] > lk_poly3_far_xhi[ ii ] ? 0 : 1 );
		}
	}

	void
	eval_lk_regular(
		VReal * __restrict__ lk_sol_regular,
		VReal * __restrict__ lk_x1,
		VReal * __restrict__ lk_x2,
		VReal * __restrict__ d,
		VReal * __restrict__ lj_radius_1,
		VReal * __restrict__ lj_radius_2,
		VReal * __restrict__ lk_inv_lambda2_1,
		VReal * __restrict__ lk_inv_lambda2_2,
		int   * __restrict__ lk_dbin,
		VReal * __restrict__ lk_coeff1,
		VReal * __restrict__ lk_coeff2,
		VReal * __restrict__ invd2,
		int limit
	)
	{
		for ( int ii = 0; ii < limit; ++ii ) {
			lk_x1[ii] = ( d[ ii ] - lj_radius_1[ ii ] ) * ( d[ ii ] - lj_radius_1[ ii ] ) * lk_inv_lambda2_1[ii];
			lk_x2[ii] = ( d[ ii ] - lj_radius_2[ ii ] ) * ( d[ ii ] - lj_radius_2[ ii ] ) * lk_inv_lambda2_2[ii];
		}

		for ( int ii = 0; ii < limit; ++ii ) {
			if ( lk_dbin[ ii ] == 3 ) {
			//if ( lk_in_mid_range[ ii ] ) {
				lk_sol_regular[ ii ] = std::exp( -1*lk_x1[ii] ) * lk_coeff1[ii] + std::exp( -1*lk_x2[ii] ) * lk_coeff2[ii];
			}
		}

		for ( int ii = 0; ii < limit; ++ii ) {
			lk_sol_regular[ ii ] *= invd2[ii];
		}
	}

	void
	evaluate_atom_pair_energies(
		core::scoring::etable::Etable const & etable,
		int jjlimit,
		EtableScratchSpace & s
	)
	{
		// here beginneth the vectoriazble code
		calc_d2( s.d2, s.x1, s.x2, s.y1, s.y2, s.z1, s.z2, jjlimit );

		calc_d( s.d, s.d2, jjlimit );
		calc_inv( s.invd2, s.d2, jjlimit );
		calc_invd6( s.invd6, s.invd2, jjlimit );

		vectorized_poly3_evaluation(
			jjlimit, s.d, s.ljatr_poly3_c0, s.ljatr_poly3_c1, s.ljatr_poly3_c2, s.ljatr_poly3_c3, s.ljatr_ramp_to_zero
		);

		// VECTORIZES! 2/27/2014 3:05PM
		VReal const famd = etable.max_dis();

		for ( int ii = 0; ii <= jjlimit; ++ii ) {
			s.lj_dbin[ ii ] =
				( 1 +
					( s.d[ ii ] > s.ljatr_poly3_xlo[ ii ] ? 1 : 0 ) +
					( s.d2[ ii ] > s.ljrep_linramp_d2_cut[ ii ] ? 1 : 0 )) *
				( s.d[ ii ] > famd ? 0 : 1 );
		}

		for ( int jj = 0; jj < jjlimit; ++jj ) {
			s.lj_closer_than_minima[ jj ] = s.d[ jj ] <= s.lj_minimum_x[ jj ];
		}

		// compute the extra repulsion (only applies between OCbb/OCbb pairs!)
		// VECTORIZES! 2/27/2014 3:10PM
		for ( int jj = 0; jj < jjlimit; ++jj ) {
			s.lj_extra_rep[ jj ] = (s.d[jj]<s.lj_exrep_xhi[jj]?1:0) * (
				(s.d[jj]<s.lj_exrep_xlo[jj]?1:0)*( (s.d[jj]-s.lj_exrep_xlo[jj] )*s.lj_exrep_extrap_slope[jj] + s.lj_exrep_ylo[jj] ) +
				(s.d[jj]<s.lj_exrep_xlo[jj]?0:1)*( (s.lj_exrep_xhi[jj]-s.d[jj])*(s.lj_exrep_xhi[jj]-s.d[jj]) * s.lj_exrep_slope[jj] ) );
		}

		// now we're ready to calculate the full lennard-jones energy, which hasn't yet been
		// split into its attractive and repulsive components
		// VECTORIZES! 2/27/2014 2:16PM
		//for ( int jj = 0; jj < jjlimit; ++jj ) {
		//	s.ljE[ jj ] = s.lj_below_linramp_start[ jj ] * s.ljrep_linramp[ jj ] +
		//		s.lj_in_mid_region[ jj ] * s.lj_regular[ jj ] +
		//		s.lj_ramp_to_zero[ jj ] * s.ljatr_ramp_to_zero[ jj ];
		//}

		for ( int ii = 0; ii < jjlimit; ++ii ) {
			s.ljE[ ii ] =
				( s.lj_dbin[ii] == 1 ? 1:0 ) * ( s.d[ ii ] * s.ljrep_linramp_slope[ ii ] + s.ljrep_linramp_intercept[ ii ] ) +
				( s.lj_dbin[ii] == 2 ? 1:0 ) * (( s.lj_r12_coeff[ ii ] * s.invd6[ii] + s.lj_r6_coeff[ii] ) * s.invd6[ii] ) +
				( s.lj_dbin[ii] == 3 ? 1:0 ) * s.ljatr_ramp_to_zero[ii];
		}

		// Next divvy it up between the attractive and repulsive components.
		// Hydrogen interactions get a "0" final weight for their attractive
		// component (unless the -Hatr flag is on the command line)
		// VECTORIZES! 2/27/2014 4:27PM
		for ( int jj = 0; jj < jjlimit; ++jj ) {
			s.ljatr[ jj ] =	s.ljatr_final_weight[ jj ] * (
				( s.ljrep_from_negcrossing[jj]?0:1 ) * (
					( s.lj_closer_than_minima[jj]?0:1 ) * s.ljE[jj] +
					( s.lj_closer_than_minima[jj]?1:0 ) * s.lj_val_at_minimum[jj] ) +
				( s.ljrep_from_negcrossing[jj]     ) * ( s.ljE[jj]<0?1:0 ) * s.ljE[jj] );
		}

		// VECTORIZES! 2/27/2014 4:27PM
		for ( int jj = 0; jj < jjlimit; ++jj ) {
			s.ljrep[ jj ] =
				(s.ljrep_from_negcrossing[jj]?0:1)*( s.lj_closer_than_minima[jj]*(s.ljE[jj] - s.lj_val_at_minimum[jj]) ) +
				(s.ljrep_from_negcrossing[jj]    )*( (s.ljE[jj]<0?0:1)*s.ljE[jj] ) +
				s.lj_extra_rep[ jj ];
		}

		// LK computations!
		calc_lk_dbin( s.lk_dbin, s.d, s.lk_poly3_close_xlo, s.lk_poly3_close_xhi, s.lk_poly3_far_xlo, s.lk_poly3_far_xhi, jjlimit );

		// 1. close poly3
		vectorized_poly3_evaluation(
			jjlimit, s.d, s.lk_poly3_close_c0, s.lk_poly3_close_c1, s.lk_poly3_close_c2, s.lk_poly3_close_c3, s.lk_close_poly3
		);

		// 2. mid-range values
		eval_lk_regular( s.lk_sol_regular, s.lk_x1, s.lk_x2,
			s.d, s.lj_radius_1, s.lj_radius_2, s.lk_inv_lambda2_1, s.lk_inv_lambda2_2,
			s.lk_dbin, s.lk_coeff1, s.lk_coeff2, s.invd2, jjlimit );

		// 3. far poly3
		vectorized_poly3_evaluation(
			jjlimit, s.d, s.lk_poly3_far_c0, s.lk_poly3_far_c1, s.lk_poly3_far_c2, s.lk_poly3_far_c3, s.lk_far_poly3
		);

		// Now put it all together!
		for ( int jj = 0; jj < jjlimit; ++jj ) {
			s.lksol[ jj ] = s.lk_final_weight[ jj ] * (
				( s.lk_dbin[jj] == 1 ? 1 : 0 ) * s.lk_poly3_close_flat[ jj ] +
				( s.lk_dbin[jj] == 2 ? 1 : 0 ) * s.lk_close_poly3[ jj ] +
				( s.lk_dbin[jj] == 3 ? 1 : 0 ) * s.lk_sol_regular[ jj ] +
				( s.lk_dbin[jj] == 4 ? 1 : 0 ) * s.lk_far_poly3[ jj ] );
		}

	}

	void
	vectorized_poly3_evaluation(
		int limit,
		VReal *__restrict__ x,
		VReal *__restrict__ c0,
		VReal *__restrict__ c1,
		VReal *__restrict__ c2,
		VReal *__restrict__ c3,
		VReal *__restrict__ y  // where the results go
	)
	{
		for ( int ii = 0; ii < limit; ++ii ) {
			y[ii] = ((x[ii]*c3[ii] + c2[ii] )*x[ii] + c1[ii] ) * x[ii] + c0[ii];
		}
	}

	void test_vectorization() {
		EnergyMethodOptions options; // default is fine
		core::scoring::etable::Etable const & etable( *( ScoringManager::get_instance()->etable( options.etable_type() )) );

		EtableScratchSpace s;


		Pose pose = create_trpcage_ideal_pose();
		core::conformation::Residue const & r1 = pose.residue(1);
		core::conformation::Residue const & r2 = pose.residue(2);

		typedef utility::fixedsizearray1< Size, 2 > Size2;
		utility::FixedSizeLexicographicalIterator< 2 > lex;

		Size2 natoms_for_respair;

		natoms_for_respair[ 1 ] = r1.natoms(); natoms_for_respair[ 2 ] = r2.natoms();
		lex.set_dimension_sizes( natoms_for_respair );

		int npairs = natoms_for_respair[ 1 ] * natoms_for_respair[ 2 ];
		for ( int ii = 1; ii <= npairs; ii += s.N ) {
			int jjlimit = std::min( ii+s.N - 1, npairs ) + 1 - ii;

			load_atom_pair_batch( etable, r1, r2, lex, ii, jjlimit, s );
			evaluate_atom_pair_energies( etable, jjlimit, s );

			// Now compare the computed ljatr and ljrep energies against those computed by the etable
			for ( int jj = 0; jj < jjlimit; ++jj ) {
				Real ljatr, ljrep, fasol, d2;
				lex.set_position_from_index( jj + ii ); // 0-based indexing

				Vector const & a1coord = r1.xyz(lex[1]);
				Vector const & a2coord = r2.xyz(lex[2]);
				d2 = a1coord.distance_squared( a2coord );
				etable.analytic_etable_evaluation( r1.atom( lex[1] ), r2.atom( lex[2] ), ljatr, ljrep, fasol, d2 );

				//std::cout << "Atoms " << lex[1] << " " << lex[2] << " " << r1.atom_name( lex[1] ) << " " << r2.atom_name( lex[2] );
				//std::cout << " correct: (" << ljatr << ", " << ljrep << ", " << fasol <<  "), actual: (";
				//std::cout << s.ljatr[jj] << ", " << s.ljrep[ jj ] << ", " << s.lksol[ jj ] << ")"  << std::endl;
				//std::cout << "   " << "lj_reg: " << s.lj_regular[jj] << " ljrep_linramp: " << s.ljrep_linramp[ jj ] << std::endl;
				//std::cout << "   " << "ljatr_ramp_to_zero: " << s.ljatr_ramp_to_zero[ jj ] << std::endl;
				//std::cout << "   " << "below_linramp_start: " << s.lj_below_linramp_start[ jj ] <<  " lj_in_mid_region: " << s.lj_in_mid_region[ jj ] << std::endl;
				//std::cout << "   " << "s.lj_ramp_to_zero: " << s.lj_ramp_to_zero[ jj ] << " s.lj_closer_than_minima[ jj ] " << s.lj_closer_than_minima[ jj ] << std::endl;
				//std::cout << "   " << "s.lj_dbin " << s.lj_dbin[ jj ] << std::endl;
				//std::cout << "   " << "ljatr_final_weight: " << s.ljatr_final_weight[ jj ] << " ljrep_from_negcrossing: " << s.ljrep_from_negcrossing[ jj ] << std::endl;
				//std::cout << "   " << "s.lj_val_at_minimum: " << s.lj_val_at_minimum[ jj ] << std::endl;
				//std::cout << "   " << "s.lj_minimum_x: " << s.lj_minimum_x[ jj ] << " d: " << s.d[ jj ] << std::endl;
				//std::cout << "   " << "s.lj_exrep_xhi: " << s.lj_exrep_xhi[ jj ] << " s.lj_exrep_xlo: " << s.lj_exrep_xlo[ jj ] << std::endl;
				//std::cout << "   " << "s.lj_extra_rep: " << s.lj_extra_rep[ jj ] << std::endl;
				//std::cout << "    lk_final_weight: " << s.lk_final_weight[ jj ] << " lk_beneath_close_poly3_xlo: " << s.lk_beneath_close_poly3_xlo[ jj ] << std::endl;
				//std::cout << "    lk_in_close_poly3_range: " << s.lk_in_close_poly3_range[ jj ] << std::endl;
				//std::cout << "    lk_in_mid_range: " <<  s.lk_in_mid_range[ jj ] << std::endl;
				//std::cout << "    lk_in_far_poly3_range: " << s.lk_in_far_poly3_range[ jj ] << std::endl;
				//std::cout << "    s.lk_poly3_close_ylo[ jj ]: " << s.lk_poly3_close_ylo[ jj ] << std::endl;
				//std::cout << "    s.lk_close_poly3[ jj ]: " << s.lk_close_poly3[ jj ] << std::endl;
				//std::cout << "    s.lk_sol_regular[ jj ]: " << s.lk_sol_regular[ jj ] << std::endl;
				//std::cout << "    s.lk_far_poly3[ jj ]: " << s.lk_far_poly3[ jj ] << std::endl;
				TS_ASSERT_DELTA( ljatr, s.ljatr[ jj ], 1e-5 );
				TS_ASSERT_DELTA( ljrep, s.ljrep[ jj ], 1e-4 );
				TS_ASSERT_DELTA( fasol, s.lksol[ jj ], 1e-5 );
			}
		}

		TS_ASSERT( true );
		//TS_ASSERT( false );
	}

	void test_vectorizer_performance() {
		using namespace core::graph;

		EnergyMethodOptions options; // default is fine
		core::scoring::etable::Etable const & etable( *( ScoringManager::get_instance()->etable( options.etable_type() )) );
		AnalyticEtableEnergy ana_lj_energy( *( ScoringManager::get_instance()->etable( options.etable_type() )), options );
		EtableScratchSpace s;

		Pose pose = create_test_in_pdb_pose(); //create_trpcage_ideal_pose();

		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( fa_atr, 0.8 );
		sfxn.set_weight( fa_rep, 0.44 );
		sfxn.set_weight( fa_sol, 0.6 );

		sfxn( pose );

		typedef utility::fixedsizearray1< Size, 2 > Size2;
		utility::FixedSizeLexicographicalIterator< 2 > lex;

		Size2 natoms_for_respair;

		core::Size const niters = 100;

		core::Size count_work_1( 0 ), count_work_2( 0 );

		clock_t start_time1 = clock();
		time_t start = time(NULL);
		for ( core::Size ii = 1; ii <= niters; ++ii ) {
			for ( core::Size jj = 1; jj <= pose.total_residue(); ++jj ) {
				core::conformation::Residue const & r1 = pose.residue( jj );
				for ( Node::EdgeListConstIter
						iter     = pose.energies().energy_graph().get_node( jj )->const_upper_edge_list_begin(),
						iter_end = pose.energies().energy_graph().get_node( jj )->const_upper_edge_list_end();
						iter != iter_end; ++iter ) {
					Size kk = (*iter)->get_second_node_ind();
					core::conformation::Residue const & r2 = pose.residue( kk );

					natoms_for_respair[ 1 ] = r1.natoms(); natoms_for_respair[ 2 ] = r2.natoms();
					lex.set_dimension_sizes( natoms_for_respair );

					int npairs = natoms_for_respair[ 1 ] * natoms_for_respair[ 2 ];
					//s.resize( npairs );

					for ( int ll = 1; ll <= npairs; ll += s.N ) {
						int mmlimit = std::min( ll+s.N - 1, npairs ) + 1 - ll;

						load_atom_pair_batch( etable, r1, r2, lex, ll, mmlimit, s );
						evaluate_atom_pair_energies( etable, mmlimit, s );

						count_work_1 += mmlimit;

					}
				}
			}
		}
		clock_t stop_time1 = clock();
		time_t stop = time(NULL);


		core::scoring::EnergyMap emap;

		Real ljatr_tot = 0, ljrep_tot = 0, fasol_tot = 0;
		clock_t start_time2 = clock();
		for ( core::Size ii = 1; ii <= niters; ++ii ) {
			for ( core::Size jj = 1; jj <= pose.total_residue(); ++jj ) {
				core::conformation::Residue const & r1 = pose.residue( jj );
				for ( Node::EdgeListConstIter
						iter     = pose.energies().energy_graph().get_node( jj )->const_upper_edge_list_begin(),
						iter_end = pose.energies().energy_graph().get_node( jj )->const_upper_edge_list_end();
						iter != iter_end; ++iter ) {
					Size kk = (*iter)->get_second_node_ind();
					core::conformation::Residue const & r2 = pose.residue( kk );
					//ana_lj_energy.residue_pair_energy( r1, r2, pose, sfxn, emap );
					for ( core::Size ll = 1; ll <= r1.natoms(); ++ll ) {
						core::conformation::Atom const & llat = r1.atom(ll);
						for ( core::Size mm = 1; mm <= r2.natoms(); ++mm ) {
							core::conformation::Atom const & mmat = r2.atom(mm);
							Real d2 = llat.xyz().distance_squared( mmat.xyz() );
							Real ljatr, ljrep, fasol;
							etable.analytic_etable_evaluation( llat, mmat, ljatr, ljrep, fasol, d2 );
							ljatr_tot += ljatr;
							ljrep_tot += ljrep;
							fasol_tot += fasol;
							++count_work_2;
						}
					}
				}
			}
		}
		clock_t stop_time2 = clock();

		TR << "Vectorized: " << ((double) stop_time1 - start_time1 )/CLOCKS_PER_SEC << " (wall time: " << stop-start << ")";
		TR << " regular: " << ((double) stop_time2 - start_time2 )/CLOCKS_PER_SEC << std::endl;

		TR << " atr: " << ljatr_tot << " rep: " << ljrep_tot << " sol: " << fasol_tot << std::endl;
		TR << " count work 1: " << count_work_1 << " count work 2: " << count_work_2 << std::endl;
	}

	VReal eval_spline_poly(
		VReal xlo, VReal xhi,
		VReal ylo, VReal yhi,
		VReal y2lo, VReal y2hi,
		VReal x
	)
	{
		VReal a = ( xhi - x ) / ( xhi - xlo );
		VReal b = ( x - xlo ) / ( xhi - xlo );
		return a*ylo + b*yhi + ( (a*a*a-a)*y2lo + (b*b*b-b)*y2hi ) * (xhi-xlo)*(xhi-xlo) / 6;
	}

	VReal
	eval_poly4(
		VReal x,
		VReal c0,
		VReal c1,
		VReal c2,
		VReal c3
	) {
		return ((c3*x+c2)*x+c1)*x + c0;
	}


	void test_convert_spline_polynomial_to_regular() {
		VReal xlo( 1.25 ), xhi( 1.75 );
		VReal ylo( 3.0 ), yhi( 1.8 ), y2lo( -0.25 ), y2hi( 0 );

		VReal a( xlo ), b( xhi ), c( yhi ), d( ylo ), e(y2hi), f(y2lo);

		VReal c0 = ( (b*b*b*f - a*a*a*e)/(b-a) + (a*e - b*f) * (b-a) )  / 6 + (  b*d - a*c ) / ( b-a );
		VReal c1 = ( 3*a*a*e/(b-a) - e*(b-a) + f*(b-a) - 3*b*b*f / (b-a) ) / 6 + ( c - d ) / (b-a);
		VReal c2 = ( 3*b*f - 3*a*e ) / ( 6 * (b-a) );
		VReal c3 = ( e-f ) / ( 6 * (b-a) );

		TR << "c0 " << c0 << std::endl;
		TR << "(  b*d - a*c ) / ( b-a ): " << (  b*d - a*c ) / ( b-a ) << std::endl;
		TR << " b*d/(b-a) " << b*d/(b-a) << std::endl;
		TR << " a*c/(b-a) " << a*c/(b-a) << std::endl;

		VReal step = 0.01;
		for ( VReal x = xlo; x < xhi; x += step ) {
			VReal ygold = eval_spline_poly( xlo, xhi, ylo, yhi, y2lo, y2hi, x );
			VReal ynew  = eval_poly4( x, c0, c1, c2, c3 );
			TS_ASSERT_DELTA( ygold, ynew, 1e-5 );
			//std::cout << "x: " << x << " ygold: " << ygold << " ynew: " << ynew << " diff: " << ygold - ynew << std::endl;
		}

	}

	void test_vectorized_square_distance_performance() {
		using namespace core::graph;

		EnergyMethodOptions options; // default is fine
		core::scoring::etable::Etable const & etable( *( ScoringManager::get_instance()->etable( options.etable_type() )) );
		AnalyticEtableEnergy ana_lj_energy( *( ScoringManager::get_instance()->etable( options.etable_type() )), options );
		EtableScratchSpace s;

		Pose pose = create_test_in_pdb_pose(); //create_trpcage_ideal_pose();

		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( fa_atr, 0.8 );
		sfxn.set_weight( fa_rep, 0.44 );
		sfxn.set_weight( fa_sol, 0.6 );

		sfxn( pose );

		EtableScratchSpace::VRealArray x1, x2, y1, y2, z1, z2, d2;
		x1 = new VReal[ EtableScratchSpace::N ];
		x2 = new VReal[ EtableScratchSpace::N ];
		y1 = new VReal[ EtableScratchSpace::N ];
		y2 = new VReal[ EtableScratchSpace::N ];
		z1 = new VReal[ EtableScratchSpace::N ];
		z2 = new VReal[ EtableScratchSpace::N ];
		d2 = new VReal[ EtableScratchSpace::N ];

		for ( int ii = 0; ii < EtableScratchSpace::N; ++ii ) {
			x1[ ii ] = x2[ ii ] = y1[ ii ] = y2[ ii ] = z1[ ii ] = z2[ ii ] = 0;
		}

		core::Size const niters = 100;

		typedef utility::fixedsizearray1< Size, 2 > Size2;
		utility::FixedSizeLexicographicalIterator< 2 > lex;

		Size2 natoms_for_respair;

		core::Size count_work_1( 0 ), count_work_2( 0 );
		double d2_1( 0 ), d2_2( 0 );

		clock_t start_time1 = clock();
		time_t start = time(NULL);
		for ( core::Size ii = 1; ii <= niters; ++ii ) {
			for ( core::Size jj = 1; jj <= pose.total_residue(); ++jj ) {
				core::conformation::Residue const & r1 = pose.residue( jj );
				for ( Node::EdgeListConstIter
						iter     = pose.energies().energy_graph().get_node( jj )->const_upper_edge_list_begin(),
						iter_end = pose.energies().energy_graph().get_node( jj )->const_upper_edge_list_end();
						iter != iter_end; ++iter ) {
					Size kk = (*iter)->get_second_node_ind();
					core::conformation::Residue const & r2 = pose.residue( kk );

					natoms_for_respair[ 1 ] = r1.natoms(); natoms_for_respair[ 2 ] = r2.natoms();
					lex.set_dimension_sizes( natoms_for_respair );

					int npairs = natoms_for_respair[ 1 ] * natoms_for_respair[ 2 ];

					for ( int ll = 1; ll <= npairs; ll += EtableScratchSpace::N ) {
						int mmlimit = std::min( ll+EtableScratchSpace::N - 1, npairs ) + 1 - ll;

						//for ( int mm = 0; mm <= mmlimit; ++mm ) {
						//	lex.set_position_from_index( mm + ll );
						//	Vector const & a1coord = r1.xyz(lex[1]);
						//	Vector const & a2coord = r2.xyz(lex[2]);
						//
						//	x1[ mm ] = a1coord.x();
						//	y1[ mm ] = a1coord.y();
						//	z1[ mm ] = a1coord.z();
						//	x2[ mm ] = a2coord.x();
						//	y2[ mm ] = a2coord.y();
						//	z2[ mm ] = a2coord.z();
						//}

						for ( int mm = 0; mm < mmlimit; ++mm ) {
							d2[ mm ] = ( x1[ mm ] - x2[ mm ] ) * ( x1[ mm ] - x2[ mm ] ) +
								( y1[ mm ] - y2[ mm ] ) * ( y1[ mm ] - y2[ mm ] ) +
								( z1[ mm ] - z2[ mm ] ) * ( z1[ mm ] - z2[ mm ] );
						}

						for ( int mm = 0; mm < mmlimit; ++mm ) d2_1 += d2[ mm ];

						count_work_1 += mmlimit;

					}
				}
			}
		}
		clock_t stop_time1 = clock();
		time_t stop = time(NULL);


		core::scoring::EnergyMap emap;

		Real ljatr_tot = 0, ljrep_tot = 0, fasol_tot = 0;
		clock_t start_time2 = clock();
		for ( core::Size ii = 1; ii <= niters; ++ii ) {
			for ( core::Size jj = 1; jj <= pose.total_residue(); ++jj ) {
				core::conformation::Residue const & r1 = pose.residue( jj );
				for ( Node::EdgeListConstIter
						iter     = pose.energies().energy_graph().get_node( jj )->const_upper_edge_list_begin(),
						iter_end = pose.energies().energy_graph().get_node( jj )->const_upper_edge_list_end();
						iter != iter_end; ++iter ) {
					Size kk = (*iter)->get_second_node_ind();
					core::conformation::Residue const & r2 = pose.residue( kk );
					//ana_lj_energy.residue_pair_energy( r1, r2, pose, sfxn, emap );
					for ( core::Size ll = 1; ll <= r1.natoms(); ++ll ) {
						core::conformation::Atom const & llat = r1.atom(ll);
						for ( core::Size mm = 1; mm <= r2.natoms(); ++mm ) {
							core::conformation::Atom const & mmat = r2.atom(mm);
							Real d2 = llat.xyz().distance_squared( mmat.xyz() );
							d2_2 += d2;
							++count_work_2;
						}
					}
				}
			}
		}
		clock_t stop_time2 = clock();

		TR << "Square distance vectorized: " << ((double) stop_time1 - start_time1 )/CLOCKS_PER_SEC << " (wall time: " << stop-start << ")";
		TR << " regular: " << ((double) stop_time2 - start_time2 )/CLOCKS_PER_SEC << std::endl;

		TR << " count work 1: " << count_work_1 << " count work 2: " << count_work_2 << std::endl;
		TR << " d2_1: " << d2_1 << " d2_2: " << d2_2 << std::endl;
	}


};
