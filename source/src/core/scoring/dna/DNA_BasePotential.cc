// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author

#include <core/scoring/dna/DNA_BasePotential.hh>

#include <core/scoring/dna/base_geometry.hh>

#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>
#include <core/kinematics/Stub.hh>
#include <basic/database/open.hh>

#include <utility/io/izstream.hh>
#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <numeric/conversions.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyz.functions.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/format.hh>


using namespace ObjexxFCL;
using namespace ObjexxFCL::format;

namespace core {
namespace scoring {
namespace dna {

/// @details Auto-generated virtual destructor
DNA_BasePotential::~DNA_BasePotential() {}

typedef ObjexxFCL::FArray1D< Real > FArray1D_Real;

DNA_BasePotential::DNA_BasePotential():
	mean_     (    6, 2, 16 ),
	stddev_   (    6, 2, 16 ),
	stiffness_( 6, 6, 2, 16 )
{
	load_score_tables();
}

std::string
DNA_BasePotential::base_string( Residue const & rsd ) const
{
	using namespace chemical;
	switch ( rsd.aa() ) {
	case na_ade:
		return "A";
	case na_cyt:
		return "C";
	case na_gua:
		return "G";
	case na_thy:
		return "T";
	default:
		utility_exit_with_message("bad rsd type for DNA_BasePotential: "+rsd.name() );
	}
	return "X";
}

///////////////////////////////////////////////////////////////////////////////////
// i1 = 1,2
// i2 = 1,16
void
DNA_BasePotential::get_array_indices( InteractionType const & t, std::string const & bases, int & i1, int & i2 ) const
{
	char const b1( bases[0] ), b2( bases[1] );
	i1 = 0; i2 = 0;
	if ( t == BP_type ) { // base pair ////////////////////////////////
		i1 = 1;
		if ( b1 == 'A' ) {
			i2 = 1;
		} else if ( b1 == 'C' ) {
			i2 = 2;
		} else if ( b1 == 'G' ) {
			i2 = 3;
		} else if ( b1 == 'T' ) {
			i2 = 4;
		}
	} else { // base step /////////////////////////////////////////////
		i1 = 2;
		if ( b1 == 'A' ) {
			i2 = 0;
		} else if ( b1 == 'C' ) {
			i2 = 4;
		} else if ( b1 == 'G' ) {
			i2 = 8;
		} else if ( b1 == 'T' ) {
			i2 = 12;
		} else {
			utility_exit_with_message( "Unknown first DNA base label: " + std::string(1, b1) );
		}

		if ( b2 == 'A' ) {
			i2 += 1;
		} else if ( b2 == 'C' ) {
			i2 += 2;
		} else if ( b2 == 'G' ) {
			i2 += 3;
		} else if ( b2 == 'T' ) {
			i2 += 4;
		} else {
			utility_exit_with_message( "Unknown first DNA base label: " + std::string(1, b2) );
		}
	}
}


///////////////////////////////////////////////////////////////////////////////
void
DNA_BasePotential::set_mean_and_stddev(
	InteractionType const & type,
	std::string const & bases,
	int const index,
	Real mean,
	Real stddev
)
{
	int i1,i2;
	get_array_indices( type, bases, i1, i2 );

	mean_  ( index, i1, i2 ) = mean;
	stddev_( index, i1, i2 ) = stddev;
}


///////////////////////////////////////////////////////////////////////////////
void
DNA_BasePotential::set_stiffness(
	InteractionType const & type,
	std::string const & bases,
	int const index1,
	int const index2,
	Real const val
)
{
	int i1,i2;
	get_array_indices( type, bases, i1, i2 );

	stiffness_( index1, index2, i1, i2 ) = val;
}

///////////////////////////////////////////////////////////////////////////////
//
void
DNA_BasePotential::load_score_tables()
{

	// load base-pair and base-step parameters
	utility::io::izstream data;
	basic::database::open( data, "scoring/dna/dna_bs_bp.dat" );


	std::string line, type_name, bases, param_name;
	bool fail( false );
	while ( getline( data,line ) ) {
		std::istringstream l( line );
		std::string tag, step;

		l >> tag;
		if ( l.fail() ) continue;
		if ( tag == "SD:" ) {
			int param_index;
			Real mean, stddev;
			l >> type_name >> bases >> param_name >> param_index >> mean >> stddev;
		debug_assert( type_name == "BS" || type_name == "BP" );
			InteractionType const type( type_name == "BS" ? BS_type : BP_type );
			if ( l.fail() || param_index < 1 || param_index > 6 ) {
				std::cout << "failline:" << line << std::endl;
				fail = true;
				break;
			}
			set_mean_and_stddev( type /* BS or BP */, bases /* eg AT */, param_index /* 1-6 */, mean, stddev );


		} else if ( tag == "stiffness:" ) {
			int i;
			Real Fij;
			std::string tmp;
			l >> type_name >> bases >> tmp >> i;
		debug_assert( type_name == "BS" || type_name == "BP" );
			InteractionType const type( type_name == "BS" ? BS_type : BP_type );
			for ( int j=1; j<= 6; ++j ) {
				l >> Fij;
				set_stiffness(  type /* BS or BP */, bases /* eg AT */, i, j, Fij );
			}
			if ( l.fail() ) {
				std::cout << "failline:" << line << std::endl;
				fail = true;
				break;
			}
		}
	}

	if ( fail ) {
		utility_exit_with_message( "Unable to find/parse the DNA basepair and basestep params -- update your mini-DB?" );
	}
}

////////////////////////////////////////////////////////////////////////////////////
// STOLEN from Alex Morozov
// private
Real
DNA_BasePotential::base_score(
	InteractionType const & type,
	std::string const & bases,
	utility::vector1< Real > const & params
) const
{
debug_assert( params.size() == 6 );

	// max allowed deviation from average beyond which the score is capped
	// PB -- this seems very large! it's in standard deviations, right?
	Real const MAX_SCORE_DEV = 100.0;

	// accumulate score values
	Real score(0.0);

	for ( int i = 1; i <= 6; ++i ) {
		Real const max_delta_i = MAX_SCORE_DEV * stddev( type, bases, i );

		Real delta_i = params[i] - mean( type, bases, i );
		// reset outliers to max allowed values
		// (beyond which the potential is unreliable)
		if ( delta_i > max_delta_i ) delta_i = max_delta_i;
		if ( delta_i < -max_delta_i ) delta_i = -max_delta_i;

		for ( int j = 1; j <= 6; ++j ) {
			Real const max_delta_j = MAX_SCORE_DEV * stddev( type, bases, j );

			Real delta_j = params[j] - mean( type, bases, j );
			// reset outliers to max allowed values
			// (beyond which the potential is unreliable)
			if ( delta_j > max_delta_j ) delta_j = max_delta_j;
			if ( delta_j < -max_delta_j ) delta_j = -max_delta_j;

			score += stiffness( type, bases, i, j ) * delta_i * delta_j;
		}
	}

	return score;
}


Real
DNA_BasePotential::base_step_score(
	Residue const & rsd1,
	Residue const & rsd2
) const
{
	utility::vector1< Real > params(6);
debug_assert( rsd2.seqpos() == rsd1.seqpos() + 1 );

	get_base_step_params( rsd1, rsd2, params );

	return base_score( BS_type, base_string(rsd1)+base_string(rsd2), params );

}

void
DNA_BasePotential::eval_base_step_Z_scores(
	Residue const & rsd1,
	Residue const & rsd2,
	utility::vector1< Real > & z_scores
) const
{
	utility::vector1< Real > params(6);
	z_scores.resize(6);
debug_assert( rsd2.seqpos() == rsd1.seqpos() + 1 );

	get_base_step_params( rsd1, rsd2, params );

	std::string const bases( base_string(rsd1)+base_string(rsd2) );

	for ( Size i=1; i<= 6; ++i ) {
		z_scores[i] = ( params[i] - mean( BS_type, bases, i ) ) / stddev( BS_type, bases, i );
	}
}


void
DNA_BasePotential::eval_base_pair_Z_scores(
	Residue const & rsd1,
	Residue const & rsd2,
	utility::vector1< Real > & z_scores
) const
{
	utility::vector1< Real > params(6);
	z_scores.resize(6);
	get_base_pair_params( rsd1, rsd2, params );

	std::string const bases( base_string(rsd1)+base_string(rsd2) );

	for ( Size i=1; i<= 6; ++i ) {
		z_scores[i] = ( params[i] - mean( BP_type, bases, i ) ) / stddev( BP_type, bases, i );
	}
}


Real
DNA_BasePotential::base_pair_score(
	Residue const & rsd1,
	Residue const & rsd2
) const
{
	utility::vector1< Real > params(6);
//debug_assert( rsd1.seqpos() < rsd2.seqpos() );
	get_base_pair_params( rsd1, rsd2, params );
	std::string const bpstr( base_string(rsd1) + base_string(rsd2) );

	// ashworth scoring potential is made highly intolerant of non-Watson-Crick basepairs. This is justified because there are currently no basepair parameter scores for non-Watson-Crick basepairs in the default database parameter file (dna_bp_bs.dat)
	if ( bpstr == "AT" || bpstr == "CG" || bpstr == "GC" || bpstr == "TA" )
		return base_score( BP_type, bpstr, params );
	else return 1e9;
}


void
DNA_BasePotential::eval_base_step_derivative(
	Residue const & rsd1,
	Residue const & rsd2,
	Vector & F1,
	Vector & F2,
	Real const external_weight_factor // should probably be +1 or -1
) const
{
	using numeric::conversions::degrees;
	using numeric::conversions::radians;
	using numeric::arccos;
	using numeric::cross;

	bool const debug( true ), verbose( false );

debug_assert( rsd1.seqpos() + 1 == rsd2.seqpos() );

	kinematics::Stub const stub1( get_base_stub( rsd1, 1 ) ), stub2( get_base_stub( rsd2, 1 ) );

	// copy matrices
	Matrix M1( stub1.M ), M2( stub2.M );

debug_assert( is_orthonormal( M1, 1e-3 ) );
debug_assert( is_orthonormal( M2, 1e-3 ) );

	if ( dot( M1.col_z(), M2.col_z() ) < 0.0 ) {
		std::cout << "dna_bs_deriv: base flip!" << std::endl;
		//utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
	}

	// get angle between the z-axes
	Real const gamma( arccos( dot( M1.col_z(), M2.col_z() ) ) );

	Vector const rt( ( cross( M2.col_z(), M1.col_z() ) ).normalized() );

	Matrix R_gamma_2( rotation_matrix( rt, gamma/2.0f ) );

	M2 = R_gamma_2 * M2;
	M1 = R_gamma_2.transposed() * M1;

debug_assert( is_orthonormal( M1, 1e-3 ) );
debug_assert( is_orthonormal( M2, 1e-3 ) );

	// build mid-base-pair triad
debug_assert( M1.col_z().distance( M2.col_z() ) < 1e-3 );
debug_assert( std::abs( dot( rt, M1.col_z() ) ) < 1e-3 );

	Matrix MBT;
	MBT.col_z( M1.col_z() );

debug_assert( std::abs( dot( M1.col_x(), MBT.col_z() ) ) < 1e-3 );
debug_assert( std::abs( dot( M2.col_x(), MBT.col_z() ) ) < 1e-3 );
debug_assert( std::abs( dot( M1.col_y(), MBT.col_z() ) ) < 1e-3 );
debug_assert( std::abs( dot( M2.col_y(), MBT.col_z() ) ) < 1e-3 );

	// get
	MBT.col_y( ( 0.5f * ( M1.col_y() + M2.col_y() ) ).normalized() );
	MBT.col_x( ( 0.5f * ( M1.col_x() + M2.col_x() ) ).normalized() );

debug_assert( is_orthonormal( MBT, 1e-3 ) );

	// angular params

	// TWIST
	// x,y,z make rh coord system
	utility::vector1< Real > params( 6, 0.0 );

	params[1] = std::atan2( dot( M1.col_x(), M2.col_y() ),
													dot( M1.col_x(), M2.col_x() ) );
	// ROLL:
	params[2] = gamma * dot( rt, MBT.col_y() );

	// TILT
	params[3] = gamma * dot( rt, MBT.col_x() );

	// translational params
	Vector const displacement( stub1.v - stub2.v );

	params[4] = dot( displacement, MBT.col_x() ); // SHIFT
	params[5] = dot( displacement, MBT.col_y() ); // SLIDE
	params[6] = dot( displacement, MBT.col_z() ); // RISE

	// check sign conventions
debug_assert( params[1] * dot( MBT.col_z(), cross( M2.col_y(), M1.col_y() ) ) > 0);

	// convert to degrees
	Real const twist( params[1] );
	Real const roll ( params[2] );
	Real const tilt ( params[3] );

	params[1] = degrees( params[1] );
	params[2] = degrees( params[2] );
	params[3] = degrees( params[3] );


	if ( debug ) {
		// doublecheck our params calculation
		utility::vector1< Real > new_params;
		get_base_step_params( rsd1, rsd2, new_params );
		for ( Size i=1; i<= 6; ++i ) {
		debug_assert( std::abs( params[i] - new_params[i] ) < 1e-2 );
		}
	}

	// max allowed deviation from average beyond which the score is capped
	// PB -- this seems very large! it's in standard deviations, right?
	Real const MAX_SCORE_DEV = 100.0;

	FArray1D_Real dE_dp(6,0.0);

	bool out_of_bounds( false );
 	std::string const bases( base_string( rsd1 ) + base_string( rsd2 ) );

	for ( int i = 1; i <= 6; ++i ) {
		Real const max_delta_i = MAX_SCORE_DEV * stddev( BS_type, bases, i );

		Real delta_i = params[i] - mean( BS_type, bases, i );
		// reset outliers to max allowed values
		// (beyond which the potential is unreliable)
		if ( delta_i > max_delta_i ) {
			out_of_bounds = true;
			delta_i = max_delta_i;
		}
		if ( delta_i < -max_delta_i ) {
			out_of_bounds = true;
			delta_i = -max_delta_i;
		}

		for ( int j = 1; j <= 6; ++j ) {
			Real const max_delta_j = MAX_SCORE_DEV * stddev( BS_type, bases, j );

			Real delta_j = params[j] - mean( BS_type, bases, j );
			// reset outliers to max allowed values
			// (beyond which the potential is unreliable)
			if ( delta_j > max_delta_j ) {
				out_of_bounds = true;
				delta_j = max_delta_j;
			}
			if ( delta_j < -max_delta_j ) {
				delta_j = -max_delta_j;
				out_of_bounds = true;
			}

			Real const Fij( stiffness( BS_type, bases, i, j ) );
			dE_dp(i) += Fij * delta_j;
			dE_dp(j) += Fij * delta_i;
		}
	}


	if ( out_of_bounds ) {
		// this is not really a big deal
		std::cout << "WARNING:: out_of_bounds in dna base-pair derivative!" <<
			std::endl;
	}


	// adjust for the external weight factor:
	for ( int i=1; i<= 6; ++i ) dE_dp(i) *= external_weight_factor;


	//////////////////////////
	// twist deriv
	//////////////////////////
	//
	// let xi be the original x-axis vector in triad i
	// xi' the vector transformed so that the z-axes coincide
	//
	// cos( twist ) = dot( x1' , x2' )
	//              = dot( x1, R_gamma( x2 ) ) // apply R_gamma/2 to both sides
	//
	// ASSUMPTION: R_gamma doesn't depend on the torsion angle, which it
	//  does! but I think the primary contribution is still captured here.
	//
	// d/dphi( dot( x1, R_gamma(x2) )) = dot( cross(u_phi,x1), R_gamma(x2) )
	// = dot( -u_phi , cross( R_gamma(x2), x1 )

	Real const rad2deg( degrees(1.0));

	{
		Real const dE_dtwist( dE_dp(1) );

		// these are for preventing d_arccos from blowing up
		static Real const small_angle( radians( 0.1f ) );
		static Real const big_angle( radians( 179.9f ) );
		static Real const max_c( std::cos( small_angle ));
		static Real const min_c( std::cos( big_angle ));

		Vector x1( stub1.M.col_x() );
		Vector R_gamma_x2( rotation_matrix( rt, gamma) * Vector(stub2.M.col_x()) );

	debug_assert( Vector(stub1.M.col_z()).distance( rotation_matrix(rt,gamma) * Vector(stub2.M.col_z())) <1e-2);

		Real c( std::cos( twist ) );
	debug_assert( std::abs( c - dot( x1,R_gamma_x2 ) ) < 1e-2 );
		c = std::min( std::max( min_c, c ), max_c );

		Real const dtwist_dc = -1 / std::sqrt( 1 - c * c );

		// since the derivatives factor through cosine we are getting
		// the derivative of an unsigned angle

		int const sign_factor( twist < 0 ? -1 : 1 );

		if ( verbose ) std::cout << "dE_dtwist: " <<
										 F(9,3, dE_dtwist ) <<
										 F(9,3, dtwist_dc ) <<
										 F(9,3, cross(R_gamma_x2,x1).length() ) << std::endl;

		F1 += sign_factor * dE_dtwist * rad2deg * dtwist_dc*
			cross( R_gamma_x2, x1 );
	}

	///////////////////
	// roll + tilt
	///////////////////
	//
	// as confirmed in the base_pair_params routine,
	// if we let tilt' = sin(gamma) * tilt / gamma,
	// then tilt' = dot( cross( z2, z1 ), x_MBT )
	//
	// likewise roll' = dot( cross( z2, z1 ), y_MBT )
	//
	// but how do x_MBT and y_MBT vary as we rotate about the phi axis??
	//
	// since these are approximately averages between a moving vector and
	// a fixed vector, we can assume:
	//
	// d/dphi( x_MBT ) = cross( u_phi, x_MBT ) / 2.0
	//
	// have to see if this actually works
	//
	// in addition we can assume that sin(gamma)/gamma is constant
	// since gamma is small


	{
		Real const dE_droll( dE_dp(2) );
		Real const dE_dtilt( dE_dp(3) );

		Vector const z1( stub1.M.col_z() ), z2( stub2.M.col_z() );

		Real const gamma_sin_gamma( gamma / std::sin( gamma ) );

		runtime_assert( std::abs( roll - gamma_sin_gamma *
											dot( cross( z2, z1 ), MBT.col_y() ) ) < 1e-2 );

		runtime_assert( std::abs( tilt - gamma_sin_gamma *
											dot( cross( z2, z1 ), MBT.col_x() ) ) < 1e-2 );

		F1 += dE_dtilt * gamma_sin_gamma * rad2deg *
			( cross( z1, cross( z2, MBT.col_x() ) ) +        // 1st term: dz1/dphi
				0.5f * cross( MBT.col_x(), cross( z1, z2 ) ) );// 2nd term: dx_MBT/dphi

		F1 += dE_droll * gamma_sin_gamma * rad2deg *
			( cross( z1, cross( z2, MBT.col_y() ) ) +
				0.5f * cross( MBT.col_y(), cross( z1, z2 ) ) );

	}

	///////////////////////
	// translational derivs
	///////////////////////
	//
	// need d/dphi( dot( M-F, x_MBT ) ) where M is rotated, F is fixed,
	// and x_MBT is varying in some complicated way
	//
	// from cst_set.cc, comments to Angle_cst::helper(...)
	// the term w x_MBT constant contributes x_MBT to F2
	// and cross( x_MBT, M) to F1
	//
	// for the other term we use d/dphi( x_MBT ) = 0.5 * cross( u_phi, X_MBT)
	// and the standard vector rearrangement
	{
		Vector M( stub1.v ), F( stub2.v );

		F1 += dE_dp(4) * ( cross( MBT.col_x(), M - 0.5f * ( M - F ) ) );
		F2 += dE_dp(4) * ( MBT.col_x() );

		F1 += dE_dp(5) * ( cross( MBT.col_y(), M - 0.5f * ( M - F ) ) );
		F2 += dE_dp(5) * ( MBT.col_y() );

		F1 += dE_dp(6) * ( cross( MBT.col_z(), M - 0.5f * ( M - F ) ) );
		F2 += dE_dp(6) * ( MBT.col_z() );

	}

} // eval_base_step_derivative ////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
void
DNA_BasePotential::eval_base_pair_derivative(
	Residue const & rsd1,
	Residue const & rsd2,
	Vector & F1,
	Vector & F2,
	Real const external_weight_factor
) const
{
	using numeric::conversions::degrees;
	using numeric::conversions::radians;
	using numeric::arccos;
	using numeric::cross;

	bool const debug( true );

debug_assert( rsd1.seqpos() < rsd2.seqpos() );

	utility::vector1< Real > params( 6, 0.0 );
	Real const MAX_SCORE_DEV = 100.0; // should be same as in dna_bp_score


	kinematics::Stub const stub1( get_base_stub( rsd1, 1 ) ), stub2( get_base_stub( rsd2, 2 ) );

	// copy matrices
	Matrix M1( stub1.M ), M2( stub2.M );

debug_assert( is_orthonormal( M1, 1e-3 ) );
debug_assert( is_orthonormal( M2, 1e-3 ) );

	if ( dot( M1.col_z(), M2.col_z() ) < 0.0 ) {
		std::cout << "dna_bp_deriv: base flip!" << std::endl;
		//utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
	}

	// get angle between the y-axes
	Real const gamma( arccos( dot( M1.col_y(), M2.col_y() ) ) );

	Vector const bo( ( cross( M2.col_y(), M1.col_y() ) ).normalized() );

	Matrix R_gamma_2( rotation_matrix( bo, gamma/2.0f ) );

	M2 = R_gamma_2 * M2;
	M1 = R_gamma_2.transposed() * M1;

debug_assert( is_orthonormal( M1, 1e-3 ) );
debug_assert( is_orthonormal( M2, 1e-3 ) );
debug_assert( M1.col_y().distance( M2.col_y() ) < 1e-3 );
debug_assert( std::abs( dot( bo, M1.col_y() ) ) < 1e-3 );

	// build mid-base-pair triad
	Matrix MBT;
	MBT.col_y( M1.col_y() );
	MBT.col_x( ( 0.5f * ( M1.col_x() + M2.col_x() ) ).normalized() );
	MBT.col_z( ( 0.5f * ( M1.col_z() + M2.col_z() ) ).normalized() );

	/////////////////
	// angular params

	// propeller
	// z,x,y make rh coord system
	Real const omega = std::atan2( dot( M1.col_z(), M2.col_x() ),
																	dot( M1.col_z(), M2.col_z() ) );

debug_assert( ( std::abs( std::abs( omega ) -
											arccos( dot( M1.col_z(), M2.col_z() ) ) )<1e-2 ) &&
					( std::abs( std::abs( omega ) -
											arccos( dot( M1.col_x(), M2.col_x() ) ) )<1e-2 ) );

	// buckle:
	Real const kappa = gamma * dot( bo, MBT.col_x() );

	// opening
	Real const sigma = gamma * dot( bo, MBT.col_z() );

	///////////////////////
	// translational params
	Vector const displacement( stub1.v - stub2.v );

	params[4] = dot( displacement, MBT.col_x() );
	params[5] = dot( displacement, MBT.col_y() );
	params[6] = dot( displacement, MBT.col_z() );

	/////////////
	// check sign conventions
debug_assert( omega * dot( MBT.col_y(), cross( M2.col_x(), M1.col_x() ) ) > 0);

	// convert to degrees
	params[1] = degrees( omega );
	params[2] = degrees( kappa );
	params[3] = degrees( sigma );

	if ( debug ) {
		// doublecheck our params calculation
		utility::vector1< Real > new_params;
		get_base_pair_params( rsd1, rsd2, new_params );
		for ( Size i=1; i<= 6; ++i ) {
		debug_assert( std::abs( params[i] - new_params[i] ) < 1e-2 );
		}
	}

	// now do deriv stuff:
	FArray1D_Real dE_dp( 6, 0.0 );

	bool out_of_bounds( false );

	std::string const bases( base_string( rsd1 ) + base_string( rsd2 ) );
	for ( int i = 1; i <= 6; ++i ) {
		Real delta_i = params[i] - mean( BP_type, bases, i );

		// reset outliers to max allowed values
		// (beyond which the potential is unreliable)
		Real const max_delta_i( MAX_SCORE_DEV*stddev(BP_type,bases,i));
		if ( delta_i >  max_delta_i ) {
			delta_i =  max_delta_i;
			out_of_bounds = true;
		}
		if ( delta_i < -max_delta_i ) {
			delta_i = -max_delta_i;
			out_of_bounds = true;
		}

		for ( int j = 1; j <= 6; ++j ) {
			Real delta_j = params[j] - mean( BP_type, bases, j );

			// reset outliers to max allowed values
			// (beyond which the potential is unreliable)
			Real const max_delta_j( MAX_SCORE_DEV*stddev(BP_type,bases,j));
			if ( delta_j >  max_delta_j ) {
				delta_j =  max_delta_j;
				out_of_bounds = true;
			}
			if ( delta_j < -max_delta_j ) {
				delta_j = -max_delta_j;
				out_of_bounds = true;
			}

			// score += Gij(ind,i,j)*delta_i*delta_j;
			Real const Fij( stiffness( BP_type, bases, i, j ));
			dE_dp(i) += Fij * delta_j;
			dE_dp(j) += Fij * delta_i;
		}
	}

	if ( out_of_bounds ) {
		// this is not really a big deal
		std::cout << "WARNING:: out_of_bounds in dna base-pair derivative!" <<
			std::endl;
	}

	// apply external weights
	for ( int i=1; i<= 6; ++i ) {
		dE_dp(i) *= external_weight_factor;
	}


	//////////////////////////
	// propeller (omega) deriv
	//////////////////////////
	//
	// let xi be the original x-axis vector in triad i
	// xi' the vector transformed so that the y-axes coincide
	//
	// cos(omega) = dot( x1' , x2' )
	//            = dot( x1, R_gamma( x2 ) ) // apply R_gamma/2 to both sides
	//
	// ASSUMPTION: R_gamma doesn't depend on the torsion angle, which it
	//  does! but I think the primary contribution is still captured here.
	//
	// d/dphi( dot( x1, R_gamma(x2) )) = dot( cross(u_phi,x1), R_gamma(x2) )
	// = dot( -u_phi , cross( R_gamma(x2), x1 )

	Real const rad2deg( degrees(1.0));

	{
		Real const dE_domega( dE_dp(1) );

		// these are for preventing d_arccos from blowing up
		static Real const small_angle( radians( 0.1f ) );
		static Real const big_angle( radians( 179.9f ) );
		static Real const max_c( std::cos( small_angle ));
		static Real const min_c( std::cos( big_angle ));

		Vector x1( stub1.M.col_x() );
		Vector R_gamma_x2( rotation_matrix( bo, gamma) * Vector(stub2.M.col_x()) );

	debug_assert( Vector(stub1.M.col_y()).distance( rotation_matrix(bo,gamma) * Vector(stub2.M.col_y())) <1e-2);

		Real c( std::cos( omega ) );
	debug_assert( std::abs( c - dot( x1,R_gamma_x2 ) ) < 1e-2 );
		c = std::min( std::max( min_c, c ), max_c );

		Real const domega_dc = -1 / std::sqrt( 1 - c * c );

		// since the derivatives factor through cosine we are getting
		// the derivative of an unsigned angle

		int const sign_factor2( omega < 0 ? -1 : 1 );

		F1 += sign_factor2 * dE_domega * rad2deg * domega_dc *
			cross( R_gamma_x2, x1 );
	}

	///////////////////
	// buckle + opening
	///////////////////
	//
	// as confirmed in the base_pair_params routine,
	// if we let kappa' = sin(gamma) * kappa / gamma,
	// then kappa' = dot( cross( y2, y1 ), x_MBT )
	//
	// likewise sigma' = dot( cross( y2, y1 ), z_MBT )
	//
	// but how do x_MBT and z_MBT vary as we rotate about the phi axis??
	//
	// since these are approximately averages between a moving vector and
	// a fixed vector, we can assume:
	//
	// d/dphi( x_MBT ) = cross( u_phi, x_MBT ) / 2.0
	//
	// have to see if this actually works
	//
	// in addition we can assume that sin(gamma)/gamma is constant
	// since gamma is small


	{
		Real const dE_dkappa( dE_dp(2) );
		Real const dE_dsigma( dE_dp(3) );

		Vector const y1( stub1.M.col_y() ), y2( stub2.M.col_y() );

		Real const gamma_sin_gamma( gamma / std::sin( gamma ) );

		F1 += dE_dkappa * gamma_sin_gamma * rad2deg *
			( cross( y1, cross( y2, MBT.col_x() ) ) +        // 1st term: dy1/dphi
				0.5f * cross( MBT.col_x(), cross( y1, y2 ) ) );// 2nd term: dx_MBT/dphi

		F1 += dE_dsigma * gamma_sin_gamma * rad2deg *
			( cross( y1, cross( y2, MBT.col_z() ) ) +
				0.5f * cross( MBT.col_z(), cross( y1, y2 ) ) );

	}

	///////////////////////
	// translational derivs
	///////////////////////
	//
	// need d/dphi( dot( M-F, x_MBT ) ) where M is rotated, F is fixed,
	// and x_MBT is varying in some complicated way
	//
	// from cst_set.cc, comments to Angle_cst::helper(...)
	// the term w x_MBT constant contributes x_MBT to F2
	// and cross( x_MBT, M) to F1
	//
	// for the other term we use d/dphi( x_MBT ) = 0.5 * cross( u_phi, X_MBT)
	// and the standard vector rearrangement
	{
		Vector M( stub1.v ), F( stub2.v );

		F1 += dE_dp(4) * ( cross( MBT.col_x(), M - 0.5f * ( M - F ) ) );
		F2 += dE_dp(4) * ( MBT.col_x() );

		F1 += dE_dp(5) * ( cross( MBT.col_y(), M - 0.5f * ( M - F ) ) );
		F2 += dE_dp(5) * ( MBT.col_y() );

		F1 += dE_dp(6) * ( cross( MBT.col_z(), M - 0.5f * ( M - F ) ) );
		F2 += dE_dp(6) * ( MBT.col_z() );

	}

} // eval_base_pair_derivative //////////////////////////


} // namespace dna
}} // scoring core
