// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

//////////////////////////////////////////////////////////////////////
///
/// @brief
/// read the header file!
///
/// @references
/// Numerical Recipes in c++ 2nd edition
/// Ralf Mueller
///
///
/// @author Steven Combs, Ralf Mueller, Jens Meiler
/// ported to Rosetta by Andrew Leaver-Fay
///
/////////////////////////////////////////////////////////////////////////

// Unit headers
#include <numeric/interpolation/spline/TricubicSpline.hh>

// Package headers
#include <numeric/types.hh>
#include <numeric/interpolation/spline/CubicSpline.hh>
#include <numeric/interpolation/spline/BicubicSpline.hh>
#include <numeric/MathVector.hh>
#include <numeric/MathMatrix.hh>

namespace numeric {
namespace interpolation {
namespace spline {

Real sqr( Real x) { return x*x; }

//! train TricubicSpline
void TricubicSpline::train
(
	const BorderFlag BORDER[3],
	const double START[3],
	const double DELTA[3],
	const MathTensor< Real> & RESULTS,
	const bool LINCONT[3],
	const std::pair< Real, Real > FIRSTBE[3]
)
{
	//check, if the points are given in positive direction
	assert( DELTA[ 0]>0 );
	assert( DELTA[ 1]>0 );
	assert( DELTA[ 2]>0 );

	//determine values for all three dimensions
	const int dimz( RESULTS.ncols());
	const int dimy( RESULTS.nrows());
	const int dimx( RESULTS.nlayers());

	//assigning values
	border_[ 0 ]  = BORDER[ 0 ];
	border_[ 1 ]  = BORDER[ 1 ];
	border_[ 2 ]  = BORDER[ 2 ];
	start_[ 0 ]   = START[ 0 ];
	start_[ 1 ]   = START[ 1 ];
	start_[ 2 ]   = START[ 2 ];
	delta_[ 0 ]   = DELTA[ 0 ];
	delta_[ 1 ]   = DELTA[ 1 ];
	delta_[ 2 ]   = DELTA[ 2 ];

	values_   = RESULTS;
	dsecox_   = RESULTS;
	dsecoy_   = RESULTS;
	dsecoz_   = RESULTS;
	dsecoxy_  = RESULTS;
	dsecoxz_  = RESULTS;
	dsecoyz_  = RESULTS;
	dsecoxyz_ = RESULTS;

	LinCont_[ 0] = LINCONT[ 0];
	LinCont_[ 1] = LINCONT[ 1];
	LinCont_[ 2] = LINCONT[ 2];
	firstbe_[ 0] = FIRSTBE[ 0];
	firstbe_[ 1] = FIRSTBE[ 1];
	firstbe_[ 2] = FIRSTBE[ 2];

	//train seven times for fxx, fyy, fzz, fxxyy, fxxzz, fyyzz, fxxyyzz
	//reduction to Spline2D by training only 2D-layers at the same time
	for ( int layer( 0); layer < dimx; ++layer ) {
		MathMatrix< Real> values( dimy, dimz);
		for ( int row( 0); row < dimy; ++row ) {
			for ( int col( 0); col < dimz; ++col ) {
				values( row, col )= values_( layer, row, col);
			}
		}
		BicubicSpline bs;
		BorderFlag border[2]                  = { BORDER[ 1],   BORDER[ 2]};
		Real start[2]                         = { START[ 1],    START[ 2]};
		Real delta[2]                         = { DELTA[ 1],    DELTA[ 2]};
		bool lin_cont[2]                      = { LINCONT[ 1],  LINCONT[ 2]};
		std::pair< Real, Real > firstbe[2]    = { FIRSTBE[ 1],  FIRSTBE[ 2]};

		bs.train( border, start, delta, values, lin_cont, firstbe);
		dsecoy_.replace_layer(  layer, bs.get_dsecox());
		dsecoz_.replace_layer(  layer, bs.get_dsecoy());
		dsecoyz_.replace_layer( layer, bs.get_dsecoxy());
	}

	for ( int row( 0); row < dimy; ++row ) {
		for ( int col( 0); col < dimz; ++col ) {
			MathVector< Real > values( dimx), dsecoz( dimx), dsecoy( dimx), dsecoyz( dimx);
			for ( int layer( 0); layer < dimx; ++layer ) {
				values(   layer) = values_(  layer, row, col);
				dsecoy(   layer) = dsecoy_(  layer, row, col);
				dsecoz(   layer) = dsecoz_(  layer, row, col);
				dsecoyz(  layer) = dsecoyz_( layer, row, col);
			}
			CubicSpline cs, csz, csy, csyz;
			cs.train(   BORDER[ 0], START[ 0], DELTA[ 0], values , FIRSTBE[ 0]);
			csy.train(  BORDER[ 0], START[ 0], DELTA[ 0], dsecoy , FIRSTBE[ 0]);
			csz.train(  BORDER[ 0], START[ 0], DELTA[ 0], dsecoz , FIRSTBE[ 0]);
			csyz.train( BORDER[ 0], START[ 0], DELTA[ 0], dsecoyz, FIRSTBE[ 0]);
			for ( int layer( 0); layer < dimx; ++layer ) {
				dsecox_(   layer, row, col ) = cs.get_dsecox()(    layer);
				dsecoxy_(  layer, row, col ) = csy.get_dsecox()(   layer);
				dsecoxz_(  layer, row, col ) = csz.get_dsecox()(   layer);
				dsecoxyz_( layer, row, col ) = csyz.get_dsecox()(  layer);
			}
		}
	}
	return;
}


//! return value at certain (x, y, z)
double TricubicSpline::F( const double x, const double y, const double z) const
{
	const int dimx( values_.nlayers());
	const int dimy( values_.nrows());
	const int dimz( values_.ncols());

	//check if argument is in range for non-periodic splines
	if ( ( border_[ 0] != e_Periodic)
			&& ( ( x < start_[ 0]) || ( start_[ 0] + ( dimx-1) * delta_[ 0] < x)) ) {
		//BCL_Assert( LinCont_[ 0], "argument out of range for non-periodic spline!");
		//BCL_Message( util::Message::e_Debug, "argument out of range, using linear continuation");
		if ( x < start_[ 0] ) {
			return F( start_[ 0], y, z) + ( x-start_[ 0])*dFdx( start_[ 0], y, z);
		}
		if ( x > start_[ 0] + ( dimx-1) * delta_[ 0] ) {
			return F( start_[ 0] + ( dimx - 1) * delta_[ 0] , y, z)
				+ ( x - start_[ 0] - ( dimx-1) * delta_[ 0])
				* dFdx( start_[ 0] + ( dimx - 1) * delta_[ 0], y, z);
		}
	}

	//check if argument y is in range for non-periodic splines
	if ( ( border_[ 1] != e_Periodic)
			&& ( ( y  < start_[ 1]) || ( start_[ 1] + ( dimy - 1) * delta_[ 1] <  y)) ) {
		//BCL_Assert( LinCont_[ 1], "argument out of range for non-periodic spline!");
		//BCL_Message( util::Message::e_Debug, "argument out of range, using linear continuation");
		if ( y < start_[ 1] ) {
			return F( x, start_[ 1], z) + ( y - start_[ 1]) * dFdy( x, start_[ 1], z);
		}
		if ( y > start_[ 1] + ( dimy - 1) * delta_[ 1] ) {
			return F( x, start_[ 1] + ( dimy-1) * delta_[ 1], z)
				+ ( y - start_[ 1] - ( dimy-1) * delta_[ 1])
				* dFdy( x, start_[ 1] + ( dimy-1) * delta_[ 1], z);
		}
	}

	//check if argument z is in range for non-periodic splines
	if ( ( border_[ 2] != e_Periodic)
			&& ( ( z < start_[ 2] ) || ( start_[ 2] + ( dimz-1) * delta_[ 2] < z)) ) {
		//BCL_Assert( LinCont_[ 2], "argument out of range for non-periodic spline!");
		//BCL_Message( util::Message::e_Debug, "argument out of range, using linear continuation");
		if ( z < start_[ 2] ) {
			return F( x, y, start_[ 2]) + ( z-start_[ 2])*dFdz(x, y, start_[ 2]);
		}
		if ( z > start_[ 2] + ( dimz - 1) * delta_[ 2] ) {
			return F( x, y, start_[ 2] + ( dimz-1) * delta_[ 2])
				+ ( z-start_[ 2] - ( dimz-1) * delta_[ 2])
				* dFdz( x, y, start_[ 2] + ( dimz-1) * delta_[ 2]);
		}
	}

	//determine i with start_[ 0]+(i-1)*delta_[ 0] < x < start_[ 0]+i*delta_[ 0] for the correct supporting points
	int i( int( floor( ( x-start_[ 0]) / delta_[ 0])));
	while ( start_[ 0] + i * delta_[ 0] < x ) i++;

	//determine j with start_[ 1]+(i-1)*delta_[ 1] < y < start_[ 1]+i*delta_[ 1] for the correct supporting points
	int j( int( floor( ( y-start_[ 1]) / delta_[ 1])));
	while ( start_[ 1] + j * delta_[ 1] < y ) j++;

	// determine k with start_[ 2]+(k-1)*delta_z < z < start_[ 2]+k*delta_z for the correct supporting points
	int k( int( floor( ( z - start_[ 2]) / delta_[ 2])));
	while ( start_[ 2] + k * delta_[ 2] < z ) k++;

	// this formula was derived from the Spline2D formula
	// the idea is to 'combine every part of the 2D formula
	// with dzm, dzp, dz3m, dz3p'

	const double delta_aktx( x - start_[ 0] - (i-1) * delta_[ 0]);
	const double delta_akty( y - start_[ 1] - (j-1) * delta_[ 1]);
	const double delta_aktz( z - start_[ 2] - (k-1) * delta_[ 2]);

	const double dxp( delta_aktx/delta_[ 0]);
	const double dxm( 1 - dxp);
	const double dx3p( ( dxp*dxp*dxp - dxp) * sqr( delta_[ 0]) / 6);
	const double dx3m( ( dxm*dxm*dxm - dxm) * sqr( delta_[ 0]) / 6);

	const double dyp( delta_akty/delta_[ 1]);
	const double dym( 1 - dyp);
	const double dy3p( ( dyp*dyp*dyp - dyp) * sqr( delta_[ 1]) / 6);
	const double dy3m( ( dym*dym*dym - dym) * sqr( delta_[ 1]) / 6);

	const double dzp( delta_aktz/delta_[ 2]);
	const double dzm( 1 - dzp);
	const double dz3p( ( dzp*dzp*dzp - dzp) * sqr( delta_[ 2]) / 6);
	const double dz3m( ( dzm*dzm*dzm - dzm) * sqr( delta_[ 2]) / 6);

	//generate positive values to prevent some problems with the indices
	while ( i<1 ) i += dimx;
	while ( j<1 ) j += dimy;
	while ( k<1 ) k += dimz;

	return
		dzm
		*(dxm*(dym*values_( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*values_( (i-1)%dimx  , j%dimy, (k-1)%dimz))
		+ dxp*(dym*values_(i%dimx       , (j-1)%dimy, (k-1)%dimz)+dyp*values_(i%dimx       , j%dimy, (k-1)%dimz))
		+dx3m*(dym*dsecox_( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*dsecox_( (i-1)%dimx  , j%dimy, (k-1)%dimz))
		+dx3p*(dym*dsecox_(i%dimx       , (j-1)%dimy, (k-1)%dimz)+dyp*dsecox_(i%dimx       , j%dimy, (k-1)%dimz))
		+ dxm*(dy3m*dsecoy_( (i-1)%dimx , (j-1)%dimy, (k-1)%dimz)+dy3p*dsecoy_( (i-1)%dimx , j%dimy, (k-1)%dimz))
		+ dxp*(dy3m*dsecoy_(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dy3p*dsecoy_(i%dimx      , j%dimy, (k-1)%dimz))
		+dx3m*(dy3m*dsecoxy_( (i-1)%dimx, (j-1)%dimy, (k-1)%dimz)+dy3p*dsecoxy_( (i-1)%dimx, j%dimy, (k-1)%dimz))
		+dx3p*(dy3m*dsecoxy_(i%dimx     , (j-1)%dimy, (k-1)%dimz)+dy3p*dsecoxy_(i%dimx     , j%dimy, (k-1)%dimz)))

		+dzp
		*(dxm*(dym*values_( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*values_( (i-1)%dimx  , j%dimy, k%dimz))
		+ dxp*(dym*values_(i%dimx       , (j-1)%dimy, k%dimz)+dyp*values_(i%dimx       , j%dimy, k%dimz))
		+dx3m*(dym*dsecox_( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*dsecox_( (i-1)%dimx  , j%dimy, k%dimz))
		+dx3p*(dym*dsecox_(i%dimx       , (j-1)%dimy, k%dimz)+dyp*dsecox_(i%dimx       , j%dimy, k%dimz))
		+ dxm*(dy3m*dsecoy_( (i-1)%dimx , (j-1)%dimy, k%dimz)+dy3p*dsecoy_( (i-1)%dimx , j%dimy, k%dimz))
		+ dxp*(dy3m*dsecoy_(i%dimx      , (j-1)%dimy, k%dimz)+dy3p*dsecoy_(i%dimx      , j%dimy, k%dimz))
		+dx3m*(dy3m*dsecoxy_( (i-1)%dimx, (j-1)%dimy, k%dimz)+dy3p*dsecoxy_( (i-1)%dimx, j%dimy, k%dimz))
		+dx3p*(dy3m*dsecoxy_(i%dimx     , (j-1)%dimy, k%dimz)+dy3p*dsecoxy_(i%dimx     , j%dimy, k%dimz)))

		+dz3m
		*(dxm*(dym*dsecoz_( (i-1)%dimx   , (j-1)%dimy, (k-1)%dimz)+dyp*dsecoz_( (i-1)%dimx   , j%dimy, (k-1)%dimz))
		+ dxp*(dym*dsecoz_(i%dimx        , (j-1)%dimy, (k-1)%dimz)+dyp*dsecoz_(i%dimx        , j%dimy, (k-1)%dimz))
		+dx3m*(dym*dsecoxz_( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*dsecoxz_( (i-1)%dimx  , j%dimy, (k-1)%dimz))
		+dx3p*(dym*dsecoxz_(i%dimx       , (j-1)%dimy, (k-1)%dimz)+dyp*dsecoxz_(i%dimx       , j%dimy, (k-1)%dimz))
		+ dxm*(dy3m*dsecoyz_( (i-1)%dimx , (j-1)%dimy, (k-1)%dimz)+dy3p*dsecoyz_( (i-1)%dimx , j%dimy, (k-1)%dimz))
		+ dxp*(dy3m*dsecoyz_(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dy3p*dsecoyz_(i%dimx      , j%dimy, (k-1)%dimz))
		+dx3m*(dy3m*dsecoxyz_( (i-1)%dimx, (j-1)%dimy, (k-1)%dimz)+dy3p*dsecoxyz_( (i-1)%dimx, j%dimy, (k-1)%dimz))
		+dx3p*(dy3m*dsecoxyz_(i%dimx     , (j-1)%dimy, (k-1)%dimz)+dy3p*dsecoxyz_(i%dimx     , j%dimy, (k-1)%dimz)))

		+dz3p
		*(dxm*(dym*dsecoz_( (i-1)%dimx   , (j-1)%dimy, k%dimz)+dyp*dsecoz_( (i-1)%dimx   , j%dimy, k%dimz))
		+ dxp*(dym*dsecoz_(i%dimx        , (j-1)%dimy, k%dimz)+dyp*dsecoz_(i%dimx        , j%dimy, k%dimz))
		+dx3m*(dym*dsecoxz_( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*dsecoxz_( (i-1)%dimx  , j%dimy, k%dimz))
		+dx3p*(dym*dsecoxz_(i%dimx       , (j-1)%dimy, k%dimz)+dyp*dsecoxz_(i%dimx       , j%dimy, k%dimz))
		+ dxm*(dy3m*dsecoyz_( (i-1)%dimx , (j-1)%dimy, k%dimz)+dy3p*dsecoyz_( (i-1)%dimx , j%dimy, k%dimz))
		+ dxp*(dy3m*dsecoyz_(i%dimx      , (j-1)%dimy, k%dimz)+dy3p*dsecoyz_(i%dimx      , j%dimy, k%dimz))
		+dx3m*(dy3m*dsecoxyz_( (i-1)%dimx, (j-1)%dimy, k%dimz)+dy3p*dsecoxyz_( (i-1)%dimx, j%dimy, k%dimz))
		+dx3p*(dy3m*dsecoxyz_(i%dimx     , (j-1)%dimy, k%dimz)+dy3p*dsecoxyz_(i%dimx     , j%dimy, k%dimz))) ;
}

//! return partial derivative at certain (x, y, z) for x
double TricubicSpline::dFdx( const double x, const double y, const double z) const
{
	const int dimx( values_.nlayers());
	const int dimy( values_.nrows());
	const int dimz( values_.ncols());

	//check if argument x is in range for non-periodic splines
	if ( ( border_[ 0] != e_Periodic)
			&& ( ( x < start_[ 0]) || ( start_[ 0] + ( dimx-1) * delta_[ 0] < x)) ) {
		//BCL_Assert( LinCont_[ 0], "argument out of range for non-periodic spline!");
		//BCL_Message( util::Message::e_Debug, "argument out of range, using linear continuation");
		if ( x < start_[ 0] ) {
			return dFdx( start_[ 0], y, z);
		}
		if ( x > start_[ 0] + ( dimx-1) * delta_[ 0] ) {
			return dFdx( start_[ 0] + ( dimx-1) * delta_[ 0], y, z);
		}
	}

	//check if argument y is in range for non-periodic splines
	if ( ( border_[ 1] != e_Periodic)
			&& ( ( y < start_[ 1] || start_[ 1] + ( dimy-1) * delta_[ 1] < y)) ) {
		//BCL_Assert( LinCont_[ 1], "argument out of range for non-periodic spline!");
		//BCL_Message( util::Message::e_Debug, "argument out of range, using linear continuation");
		if ( y < start_[ 1] ) {
			return dFdx(x, start_[ 1], z);
		}
		if ( y > start_[ 1] + ( dimy - 1) * delta_[ 1] ) {
			return dFdx( x, start_[ 1] + ( dimy-1) * delta_[ 1], z);
		}
	}

	//check if argument z is in range for non-periodic splines
	if ( ( border_[ 2] != e_Periodic )
			&& ( ( z < start_[ 2] || start_[ 2] + ( dimz-1) * delta_[ 2] < z)) ) {
		//BCL_Assert( LinCont_[ 2], "argument out of range for non-periodic spline!");
		//BCL_Message( util::Message::e_Debug, "argument out of range, using linear continuation");
		if ( z < start_[ 2] ) {
			return dFdx( x, y, start_[ 2]);
		}
		if ( z > start_[ 2] + ( dimz-1) * delta_[ 2] ) {
			return dFdx( x, y, start_[ 2] + ( dimz-1) * delta_[ 2]);
		}
	}

	//determine i with start_+(i-1)*delta_ < x < start_+i*delta_ for the correct supporting points
	int i( int( floor( ( x - start_[ 0]) / delta_[ 0])));
	while ( start_[ 0] + i * delta_[ 0] < x ) i++;

	//the same for j
	int j( int( floor( ( y-start_[ 1]) / delta_[ 1])));
	while ( start_[ 1] + j * delta_[ 1] < y ) j++;

	//the same for k
	int k( int( floor( ( z - start_[ 2]) / delta_[ 2])));
	while ( start_[ 2] + k * delta_[ 2] < z ) k++;

	const double delta_aktx( x - start_[ 0] - ( i - 1) * delta_[ 0]);
	const double delta_akty( y - start_[ 1] - ( j - 1) * delta_[ 1]);
	const double delta_aktz( z - start_[ 2] - ( k - 1) * delta_[ 2]);

	const double dxp( delta_aktx / delta_[ 0]);
	const double dxm( 1 - dxp);

	const double dyp( delta_akty / delta_[ 1]);
	const double dym( 1 - dyp);
	const double dy3p( ( dyp * dyp * dyp - dyp) * sqr( delta_[ 1]) / 6);
	const double dy3m( ( dym * dym * dym - dym) * sqr( delta_[ 1]) / 6);

	const double dzp( delta_aktz / delta_[ 2]);
	const double dzm( 1 - dzp);
	const double dz3p( ( dzp * dzp * dzp - dzp) * sqr( delta_[ 2]) / 6);
	const double dz3m( ( dzm * dzm * dzm - dzm) * sqr( delta_[ 2]) / 6);

	//generate positive values to prevent some problems with the indizes
	while ( i<1 ) i += dimx;
	while ( j<1 ) j += dimy;
	while ( k<1 ) k += dimz;

	return
		dzm
		*(
		-(dym*values_( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*values_( (i-1)%dimx  , j%dimy, (k-1)%dimz))/delta_[ 0]
		+(dym*values_(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*values_(i%dimx      , j%dimy, (k-1)%dimz))/delta_[ 0]
		- (3 * dxm*dxm - 1) * delta_[ 0] / 6*(dym*dsecox_( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*dsecox_( (i-1)%dimx  , j%dimy, (k-1)%dimz))
		+ (3 * dxp*dxp - 1) * delta_[ 0] / 6*(dym*dsecox_(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*dsecox_(i%dimx      , j%dimy, (k-1)%dimz))
		-(dy3m*dsecoy_( (i-1)%dimx , (j-1)%dimy, (k-1)%dimz)+dy3p*dsecoy_( (i-1)%dimx , j%dimy, (k-1)%dimz))/delta_[ 0]
		+(dy3m*dsecoy_(i%dimx     , (j-1)%dimy, (k-1)%dimz)+dy3p*dsecoy_(i%dimx     , j%dimy, (k-1)%dimz))/delta_[ 0]
		- (3 * dxm*dxm - 1) * delta_[ 0] / 6*(dy3m*dsecoxy_( (i-1)%dimx, (j-1)%dimy, (k-1)%dimz)+dy3p*dsecoxy_( (i-1)%dimx, j%dimy, (k-1)%dimz))
		+ (3 * dxp*dxp - 1) * delta_[ 0] / 6*(dy3m*dsecoxy_(i%dimx    , (j-1)%dimy, (k-1)%dimz)+dy3p*dsecoxy_(i%dimx    , j%dimy, (k-1)%dimz))
		)

		+dzp
		*(
		-(dym*values_( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*values_( (i-1)%dimx  , j%dimy, k%dimz))/delta_[ 0]
		+(dym*values_(i%dimx      , (j-1)%dimy, k%dimz)+dyp*values_(i%dimx      , j%dimy, k%dimz))/delta_[ 0]
		- (3 * dxm*dxm - 1) * delta_[ 0] / 6*(dym*dsecox_( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*dsecox_( (i-1)%dimx  , j%dimy, k%dimz))
		+ (3 * dxp*dxp - 1) * delta_[ 0] / 6*(dym*dsecox_(i%dimx      , (j-1)%dimy, k%dimz)+dyp*dsecox_(i%dimx      , j%dimy, k%dimz))
		-(dy3m*dsecoy_( (i-1)%dimx , (j-1)%dimy, k%dimz)+dy3p*dsecoy_( (i-1)%dimx , j%dimy, k%dimz))/delta_[ 0]
		+(dy3m*dsecoy_(i%dimx     , (j-1)%dimy, k%dimz)+dy3p*dsecoy_(i%dimx     , j%dimy, k%dimz))/delta_[ 0]
		- (3 * dxm*dxm - 1) * delta_[ 0] / 6*(dy3m*dsecoxy_( (i-1)%dimx, (j-1)%dimy, k%dimz)+dy3p*dsecoxy_( (i-1)%dimx, j%dimy, k%dimz))
		+ (3 * dxp*dxp - 1) * delta_[ 0] / 6*(dy3m*dsecoxy_(i%dimx    , (j-1)%dimy, k%dimz)+dy3p*dsecoxy_(i%dimx    , j%dimy, k%dimz))
		)

		+dz3m
		*(
		-(dym*dsecoz_( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*dsecoz_( (i-1)%dimx  , j%dimy, (k-1)%dimz))/delta_[ 0]
		+(dym*dsecoz_(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*dsecoz_(i%dimx      , j%dimy, (k-1)%dimz))/delta_[ 0]
		- (3 * dxm*dxm - 1) * delta_[ 0] / 6*(dym*dsecoxz_( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*dsecoxz_( (i-1)%dimx  , j%dimy, (k-1)%dimz))
		+ (3 * dxp*dxp - 1) * delta_[ 0] / 6*(dym*dsecoxz_(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*dsecoxz_(i%dimx      , j%dimy, (k-1)%dimz))
		-(dy3m*dsecoyz_( (i-1)%dimx , (j-1)%dimy, (k-1)%dimz)+dy3p*dsecoyz_( (i-1)%dimx , j%dimy, (k-1)%dimz))/delta_[ 0]
		+(dy3m*dsecoyz_(i%dimx     , (j-1)%dimy, (k-1)%dimz)+dy3p*dsecoyz_(i%dimx     , j%dimy, (k-1)%dimz))/delta_[ 0]
		- (3 * dxm*dxm - 1) * delta_[ 0] / 6*(dy3m*dsecoxyz_( (i-1)%dimx, (j-1)%dimy, (k-1)%dimz)+dy3p*dsecoxyz_( (i-1)%dimx, j%dimy, (k-1)%dimz))
		+ (3 * dxp*dxp - 1) * delta_[ 0] / 6*(dy3m*dsecoxyz_(i%dimx    , (j-1)%dimy, (k-1)%dimz)+dy3p*dsecoxyz_(i%dimx    , j%dimy, (k-1)%dimz))
		)

		+dz3p
		*(
		-(dym*dsecoz_( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*dsecoz_( (i-1)%dimx  , j%dimy, k%dimz))/delta_[ 0]
		+(dym*dsecoz_(i%dimx      , (j-1)%dimy, k%dimz)+dyp*dsecoz_(i%dimx      , j%dimy, k%dimz))/delta_[ 0]
		- (3 * dxm*dxm - 1) * delta_[ 0] / 6*(dym*dsecoxz_( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*dsecoxz_( (i-1)%dimx  , j%dimy, k%dimz))
		+ (3 * dxp*dxp - 1) * delta_[ 0] / 6*(dym*dsecoxz_(i%dimx      , (j-1)%dimy, k%dimz)+dyp*dsecoxz_(i%dimx      , j%dimy, k%dimz))
		-(dy3m*dsecoyz_( (i-1)%dimx , (j-1)%dimy, k%dimz)+dy3p*dsecoyz_( (i-1)%dimx , j%dimy, k%dimz))/delta_[ 0]
		+(dy3m*dsecoyz_(i%dimx     , (j-1)%dimy, k%dimz)+dy3p*dsecoyz_(i%dimx     , j%dimy, k%dimz))/delta_[ 0]
		- (3 * dxm*dxm - 1) * delta_[ 0] / 6*(dy3m*dsecoxyz_( (i-1)%dimx, (j-1)%dimy, k%dimz)+dy3p*dsecoxyz_( (i-1)%dimx, j%dimy, k%dimz))
		+ (3 * dxp*dxp - 1) * delta_[ 0] / 6*(dy3m*dsecoxyz_(i%dimx    , (j-1)%dimy, k%dimz)+dy3p*dsecoxyz_(i%dimx    , j%dimy, k%dimz)));
}

//! return partial derivative at certain (x, y, z) for y
double TricubicSpline::dFdy( const double x, const double y, const double z) const
{
	const int dimx( values_.nlayers());
	const int dimy( values_.nrows());
	const int dimz( values_.ncols());

	//check if argument x is in range for non-periodic splines
	if ( ( border_[ 0] != e_Periodic )
			&& ( ( x < start_[ 0]) || (start_[ 0] + ( dimx-1) * delta_[ 0] < x)) ) {
		//BCL_Assert( LinCont_[ 0], "argument out of range for non-periodic spline!");
		//BCL_Message( util::Message::e_Debug, "argument out of range, using linear continuation");
		if ( x < start_[ 0] ) {
			return dFdy( start_[ 0], y, z);
		}
		if ( x > start_[ 0] + ( dimx-1) * delta_[ 0] ) {
			return dFdy( start_[ 0] + ( dimx-1) * delta_[ 0], y, z);
		}
	}

	//check if argument y is in range for non-periodic splines
	if ( ( border_[ 1] != e_Periodic )
			&& ( ( y < start_[ 1] || start_[ 1] + ( dimy-1) * delta_[ 1] < y)) ) {
		//BCL_Assert( LinCont_[ 1], "argument out of range for non-periodic spline!");
		//BCL_Message( util::Message::e_Debug, "argument out of range, using linear continuation");
		if ( y < start_[ 1] ) {
			return dFdy(x, start_[ 1], z);
		}
		if ( y > start_[ 1] + ( dimy-1) * delta_[ 1] ) {
			return dFdy(x, start_[ 1] + ( dimy - 1) * delta_[ 1], z);
		}
	}

	//check if argument z is in range for non-periodic splines
	if ( ( border_[ 2] != e_Periodic)
			&& ( ( z < start_[ 2] || start_[ 2] + ( dimz - 1) * delta_[ 2] < z )) ) {
		//BCL_Assert( LinCont_[ 2], "argument out of range for non-periodic spline!");
		//BCL_Message( util::Message::e_Debug, "argument out of range, using linear continuation");
		if ( z < start_[ 2] ) {
			return dFdy( x, y, start_[ 2]);
		}
		if ( z > start_[ 2] + ( dimz-1) * delta_[ 2] ) {
			return dFdy( x, y, start_[ 2] + ( dimz-1) * delta_[ 2]);
		}
	}

	//determine i with start_[ 0]+(i-1)*delta_[ 0] < x < start_[ 0]+i*delta_[ 0]
	//for the correct supporting points
	int i( int( floor( ( x - start_[ 0]) / delta_[ 0])));
	while ( start_[ 0] + i * delta_[ 0] < x ) i++;

	//determine j with start_[ 1]+(j-1)*delta_[ 1] < y < start_[ 1]+j*delta_[ 1]
	//for the correct supporting points
	int j( int( floor( ( y - start_[ 1]) / delta_[ 1])));
	while ( start_[ 1] + j * delta_[ 1] < y ) j++;

	//determine k with start_[ 2]+(k-1)*delta_[ 2] < z < start_[ 2]+k*delta_[ 2]
	//for the correct supporting points
	int k( int( floor( ( z - start_[ 2]) / delta_[ 2])));
	while ( start_[ 2] + k * delta_[ 2] < z ) k++;

	//generate some auxiliary variables
	const double delta_aktx( x - start_[ 0] - ( i - 1) * delta_[ 0]);
	const double delta_akty( y - start_[ 1] - ( j - 1) * delta_[ 1]);
	const double delta_aktz( z - start_[ 2] - ( k - 1) * delta_[ 2]);

	//see F(x, y, z) for a short explanation of the values
	const double dxp( delta_aktx / delta_[ 0]);
	const double dxm( 1 - dxp);
	const double dx3p( ( dxp * dxp * dxp - dxp) * sqr( delta_[ 0]) / 6);
	const double dx3m( ( dxm * dxm * dxm - dxm) * sqr( delta_[ 0]) / 6);

	const double dyp( delta_akty / delta_[ 1]);
	const double dym( 1 - dyp);

	const double dzp( delta_aktz / delta_[ 2]);
	const double dzm( 1 - dzp);
	const double dz3p( ( dzp * dzp * dzp - dzp) * sqr( delta_[ 2]) / 6);
	const double dz3m( ( dzm * dzm * dzm - dzm) * sqr( delta_[ 2]) / 6);

	//generate positive values to prevent some problems with the indizes
	while ( i <= 0 ) i += dimx;
	while ( j <= 0 ) j += dimy;
	while ( k <= 0 ) k += dimz;

	return
		dzm
		*(
		dxm*(-values_( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+values_( (i-1)%dimx  , j%dimy, (k-1)%dimz))/delta_[ 1]
		+ dxp*(-values_(i%dimx      , (j-1)%dimy, (k-1)%dimz)+values_(i%dimx      , j%dimy, (k-1)%dimz))/delta_[ 1]
		+dx3m*(-dsecox_( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dsecox_( (i-1)%dimx  , j%dimy, (k-1)%dimz))/delta_[ 1]
		+dx3p*(-dsecox_(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dsecox_(i%dimx      , j%dimy, (k-1)%dimz))/delta_[ 1]
		+ dxm*(-(3*dym*dym-1)*dsecoy_( (i-1)%dimx , (j-1)%dimy, (k-1)%dimz)+(3*dyp*dyp-1)*dsecoy_( (i-1)%dimx , j%dimy, (k-1)%dimz))* delta_[ 1]/ 6
		+ dxp*(-(3*dym*dym-1)*dsecoy_(i%dimx     , (j-1)%dimy, (k-1)%dimz)+(3*dyp*dyp-1)*dsecoy_(i%dimx     , j%dimy, (k-1)%dimz))* delta_[ 1]/ 6
		+dx3m*(-(3*dym*dym-1)*dsecoxy_( (i-1)%dimx, (j-1)%dimy, (k-1)%dimz)+(3*dyp*dyp-1)*dsecoxy_( (i-1)%dimx, j%dimy, (k-1)%dimz))* delta_[ 1]/ 6
		+dx3p*(-(3*dym*dym-1)*dsecoxy_(i%dimx    , (j-1)%dimy, (k-1)%dimz)+(3*dyp*dyp-1)*dsecoxy_(i%dimx    , j%dimy, (k-1)%dimz))* delta_[ 1]/ 6
		)

		+dzp
		*(
		dxm*(-values_( (i-1)%dimx  , (j-1)%dimy, k%dimz)+values_( (i-1)%dimx  , j%dimy, k%dimz))/delta_[ 1]
		+ dxp*(-values_(i%dimx      , (j-1)%dimy, k%dimz)+values_(i%dimx      , j%dimy, k%dimz))/delta_[ 1]
		+dx3m*(-dsecox_( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dsecox_( (i-1)%dimx  , j%dimy, k%dimz))/delta_[ 1]
		+dx3p*(-dsecox_(i%dimx      , (j-1)%dimy, k%dimz)+dsecox_(i%dimx      , j%dimy, k%dimz))/delta_[ 1]
		+ dxm*(-(3*dym*dym-1)*dsecoy_( (i-1)%dimx , (j-1)%dimy, k%dimz)+(3*dyp*dyp-1)*dsecoy_( (i-1)%dimx , j%dimy, k%dimz))* delta_[ 1]/ 6
		+ dxp*(-(3*dym*dym-1)*dsecoy_(i%dimx     , (j-1)%dimy, k%dimz)+(3*dyp*dyp-1)*dsecoy_(i%dimx     , j%dimy, k%dimz))* delta_[ 1]/ 6
		+dx3m*(-(3*dym*dym-1)*dsecoxy_( (i-1)%dimx, (j-1)%dimy, k%dimz)+(3*dyp*dyp-1)*dsecoxy_( (i-1)%dimx, j%dimy, k%dimz))* delta_[ 1]/ 6
		+dx3p*(-(3*dym*dym-1)*dsecoxy_(i%dimx    , (j-1)%dimy, k%dimz)+(3*dyp*dyp-1)*dsecoxy_(i%dimx    , j%dimy, k%dimz))* delta_[ 1]/ 6
		)

		+dz3m
		*(
		dxm*(-dsecoz_( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dsecoz_( (i-1)%dimx  , j%dimy, (k-1)%dimz))/delta_[ 1]
		+ dxp*(-dsecoz_(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dsecoz_(i%dimx      , j%dimy, (k-1)%dimz))/delta_[ 1]
		+dx3m*(-dsecoxz_( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dsecoxz_( (i-1)%dimx  , j%dimy, (k-1)%dimz))/delta_[ 1]
		+dx3p*(-dsecoxz_(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dsecoxz_(i%dimx      , j%dimy, (k-1)%dimz))/delta_[ 1]
		+ dxm*(-(3*dym*dym-1)*dsecoyz_( (i-1)%dimx , (j-1)%dimy, (k-1)%dimz)+(3*dyp*dyp-1)*dsecoyz_( (i-1)%dimx , j%dimy, (k-1)%dimz))* delta_[ 1]/ 6
		+ dxp*(-(3*dym*dym-1)*dsecoyz_(i%dimx     , (j-1)%dimy, (k-1)%dimz)+(3*dyp*dyp-1)*dsecoyz_(i%dimx     , j%dimy, (k-1)%dimz))* delta_[ 1]/ 6
		+dx3m*(-(3*dym*dym-1)*dsecoxyz_( (i-1)%dimx, (j-1)%dimy, (k-1)%dimz)+(3*dyp*dyp-1)*dsecoxyz_( (i-1)%dimx, j%dimy, (k-1)%dimz))* delta_[ 1]/ 6
		+dx3p*(-(3*dym*dym-1)*dsecoxyz_(i%dimx    , (j-1)%dimy, (k-1)%dimz)+(3*dyp*dyp-1)*dsecoxyz_(i%dimx    , j%dimy, (k-1)%dimz))* delta_[ 1]/ 6
		)

		+dz3p
		*(
		dxm*(-dsecoz_( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dsecoz_( (i-1)%dimx  , j%dimy, k%dimz))/delta_[ 1]
		+ dxp*(-dsecoz_(i%dimx      , (j-1)%dimy, k%dimz)+dsecoz_(i%dimx      , j%dimy, k%dimz))/delta_[ 1]
		+dx3m*(-dsecoxz_( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dsecoxz_( (i-1)%dimx  , j%dimy, k%dimz))/delta_[ 1]
		+dx3p*(-dsecoxz_(i%dimx      , (j-1)%dimy, k%dimz)+dsecoxz_(i%dimx      , j%dimy, k%dimz))/delta_[ 1]
		+ dxm*(-(3*dym*dym-1)*dsecoyz_( (i-1)%dimx , (j-1)%dimy, k%dimz)+(3*dyp*dyp-1)*dsecoyz_( (i-1)%dimx , j%dimy, k%dimz))* delta_[ 1]/ 6
		+ dxp*(-(3*dym*dym-1)*dsecoyz_(i%dimx     , (j-1)%dimy, k%dimz)+(3*dyp*dyp-1)*dsecoyz_(i%dimx     , j%dimy, k%dimz))* delta_[ 1]/ 6
		+dx3m*(-(3*dym*dym-1)*dsecoxyz_( (i-1)%dimx, (j-1)%dimy, k%dimz)+(3*dyp*dyp-1)*dsecoxyz_( (i-1)%dimx, j%dimy, k%dimz))* delta_[ 1]/ 6
		+dx3p*(-(3*dym*dym-1)*dsecoxyz_(i%dimx    , (j-1)%dimy, k%dimz)+(3*dyp*dyp-1)*dsecoxyz_(i%dimx    , j%dimy, k%dimz))* delta_[ 1]/ 6);
}

//! return partial derivative at certain (x, y, z) for z
double TricubicSpline::dFdz( const double x, const double y, const double z) const
{
	const int dimx( values_.nlayers());
	const int dimy( values_.nrows());
	const int dimz( values_.ncols());

	//check if argument x is in range for non-periodic splines
	if ( ( border_[ 0] != e_Periodic)
			&& ( ( x < start_[ 0]) || ( start_[ 0] + ( dimx - 1) * delta_[ 0] < x)) ) {
		//BCL_Assert( LinCont_[ 0], "argument out of range for non-periodic spline!");
		//BCL_Message( util::Message::e_Debug, "argument out of range, using linear continuation");
		if ( x < start_[ 0] ) {
			return dFdz( start_[ 0], y, z);
		}
		if ( x > start_[ 0] + ( dimx - 1) * delta_[ 0] ) {
			return dFdz( start_[ 0] + ( dimx - 1) * delta_[ 0], y, z);
		}
	}

	//check if argument y is in range for non-periodic splines
	if ( ( border_[ 1] != e_Periodic)
			&& ( ( y < start_[ 1] || start_[ 1] + ( dimy - 1) * delta_[ 1] < y )) ) {
		//BCL_Assert( LinCont_[ 1], "argument out of range for non-periodic spline!");
		//BCL_Message( util::Message::e_Debug, "argument out of range, using linear continuation");
		if ( y < start_[ 1] ) {
			return dFdz( x, start_[ 1], z);
		}
		if ( y > start_[ 1] + ( dimy - 1) * delta_[ 1] ) {
			return dFdz( x, start_[ 1] + ( dimy - 1) * delta_[ 1], z);
		}
	}

	//check if argument z is in range for non-periodic splines
	if ( ( border_[ 2] != e_Periodic)
			&& ( ( z < start_[ 2]) || ( start_[ 2] + ( dimz - 1) * delta_[ 2] < z)) ) {
		//BCL_Assert( LinCont_[ 2], "argument out of range for non-periodic spline!");
		//BCL_Message( util::Message::e_Debug, "argument out of range, using linear continuation");
		if ( z < start_[ 2] ) {
			return dFdz(x, y, start_[ 2]);
		}
		if ( z > start_[ 2] + ( dimz - 1) * delta_[ 2] ) {
			return dFdz( x, y, start_[ 2] + ( dimz - 1) * delta_[ 2]);
		}
	}

	//determine i with start_+(i-1)*delta_ < x < start_+i*delta_ for the correct supporting points
	int i( int( floor( ( x - start_[ 0]) / delta_[ 0])));
	while ( start_[ 0] + i * delta_[ 0] < x ) i++;

	//the same for j and y
	int j( int( floor( ( y - start_[ 1]) / delta_[ 1])));
	while ( start_[ 1] + j * delta_[ 1] < y ) j++;

	//the same for k and z
	int k( int( floor( ( z - start_[ 2]) / delta_[ 2])));
	while ( start_[ 2] + k * delta_[ 2] < z ) k++;

	//generate some auxiliary variables
	const double delta_aktx( x-start_[ 0]-(i-1)*delta_[ 0]);
	const double delta_akty( y-start_[ 1]-(j-1)*delta_[ 1]);
	const double delta_aktz( z-start_[ 2]-(k-1)*delta_[ 2]);

	const double dxp( delta_aktx / delta_[ 0]);
	const double dxm( 1 - dxp);
	const double dx3p( ( dxp * dxp * dxp - dxp) * sqr( delta_[ 0]) / 6);
	const double dx3m( ( dxm * dxm * dxm - dxm) * sqr( delta_[ 0]) / 6);

	const double dyp( delta_akty / delta_[ 1]);
	const double dym( 1 - dyp);
	const double dy3p( ( dyp * dyp * dyp - dyp) * sqr( delta_[ 1]) / 6);
	const double dy3m( ( dym * dym * dym - dym) * sqr( delta_[ 1]) / 6);

	const double dzp( delta_aktz / delta_[ 2]);
	const double dzm( 1 - dzp);

	//generate positive values to prevent some problems with the indizes
	while ( i <= 0 ) i += dimx;
	while ( j <= 0 ) j += dimy;
	while ( k <= 0 ) k += dimz;

	return
		-(dxm*(dym*values_( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*values_( (i-1)%dimx  , j%dimy, (k-1)%dimz))
		+ dxp*(dym*values_(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*values_(i%dimx      , j%dimy, (k-1)%dimz))
		+dx3m*(dym*dsecox_( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*dsecox_( (i-1)%dimx  , j%dimy, (k-1)%dimz))
		+dx3p*(dym*dsecox_(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*dsecox_(i%dimx      , j%dimy, (k-1)%dimz))
		+ dxm*(dy3m*dsecoy_( (i-1)%dimx , (j-1)%dimy, (k-1)%dimz)+dy3p*dsecoy_( (i-1)%dimx , j%dimy, (k-1)%dimz))
		+ dxp*(dy3m*dsecoy_(i%dimx     , (j-1)%dimy, (k-1)%dimz)+dy3p*dsecoy_(i%dimx     , j%dimy, (k-1)%dimz))
		+dx3m*(dy3m*dsecoxy_( (i-1)%dimx, (j-1)%dimy, (k-1)%dimz)+dy3p*dsecoxy_( (i-1)%dimx, j%dimy, (k-1)%dimz))
		+dx3p*(dy3m*dsecoxy_(i%dimx    , (j-1)%dimy, (k-1)%dimz)+dy3p*dsecoxy_(i%dimx    , j%dimy, (k-1)%dimz)))
		/delta_[ 2]

		+(dxm*(dym*values_( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*values_( (i-1)%dimx  , j%dimy, k%dimz))
		+ dxp*(dym*values_(i%dimx      , (j-1)%dimy, k%dimz)+dyp*values_(i%dimx      , j%dimy, k%dimz))
		+dx3m*(dym*dsecox_( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*dsecox_( (i-1)%dimx  , j%dimy, k%dimz))
		+dx3p*(dym*dsecox_(i%dimx      , (j-1)%dimy, k%dimz)+dyp*dsecox_(i%dimx      , j%dimy, k%dimz))
		+ dxm*(dy3m*dsecoy_( (i-1)%dimx , (j-1)%dimy, k%dimz)+dy3p*dsecoy_( (i-1)%dimx , j%dimy, k%dimz))
		+ dxp*(dy3m*dsecoy_(i%dimx     , (j-1)%dimy, k%dimz)+dy3p*dsecoy_(i%dimx     , j%dimy, k%dimz))
		+dx3m*(dy3m*dsecoxy_( (i-1)%dimx, (j-1)%dimy, k%dimz)+dy3p*dsecoxy_( (i-1)%dimx, j%dimy, k%dimz))
		+dx3p*(dy3m*dsecoxy_(i%dimx    , (j-1)%dimy, k%dimz)+dy3p*dsecoxy_(i%dimx    , j%dimy, k%dimz)))
		/delta_[ 2]

		-(3*dzm*dzm-1)*delta_[ 2]/ 6
		*(dxm*(dym*dsecoz_( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*dsecoz_( (i-1)%dimx  , j%dimy, (k-1)%dimz))
		+ dxp*(dym*dsecoz_(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*dsecoz_(i%dimx      , j%dimy, (k-1)%dimz))
		+dx3m*(dym*dsecoxz_( (i-1)%dimx  , (j-1)%dimy, (k-1)%dimz)+dyp*dsecoxz_( (i-1)%dimx  , j%dimy, (k-1)%dimz))
		+dx3p*(dym*dsecoxz_(i%dimx      , (j-1)%dimy, (k-1)%dimz)+dyp*dsecoxz_(i%dimx      , j%dimy, (k-1)%dimz))
		+ dxm*(dy3m*dsecoyz_( (i-1)%dimx , (j-1)%dimy, (k-1)%dimz)+dy3p*dsecoyz_( (i-1)%dimx , j%dimy, (k-1)%dimz))
		+ dxp*(dy3m*dsecoyz_(i%dimx     , (j-1)%dimy, (k-1)%dimz)+dy3p*dsecoyz_(i%dimx     , j%dimy, (k-1)%dimz))
		+dx3m*(dy3m*dsecoxyz_( (i-1)%dimx, (j-1)%dimy, (k-1)%dimz)+dy3p*dsecoxyz_( (i-1)%dimx, j%dimy, (k-1)%dimz))
		+dx3p*(dy3m*dsecoxyz_(i%dimx    , (j-1)%dimy, (k-1)%dimz)+dy3p*dsecoxyz_(i%dimx    , j%dimy, (k-1)%dimz)))

		+(3*dzp*dzp-1)*delta_[ 2]/ 6
		*(dxm*(dym*dsecoz_( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*dsecoz_( (i-1)%dimx  , j%dimy, k%dimz))
		+ dxp*(dym*dsecoz_(i%dimx      , (j-1)%dimy, k%dimz)+dyp*dsecoz_(i%dimx      , j%dimy, k%dimz))
		+dx3m*(dym*dsecoxz_( (i-1)%dimx  , (j-1)%dimy, k%dimz)+dyp*dsecoxz_( (i-1)%dimx  , j%dimy, k%dimz))
		+dx3p*(dym*dsecoxz_(i%dimx      , (j-1)%dimy, k%dimz)+dyp*dsecoxz_(i%dimx      , j%dimy, k%dimz))
		+ dxm*(dy3m*dsecoyz_( (i-1)%dimx , (j-1)%dimy, k%dimz)+dy3p*dsecoyz_( (i-1)%dimx , j%dimy, k%dimz))
		+ dxp*(dy3m*dsecoyz_(i%dimx     , (j-1)%dimy, k%dimz)+dy3p*dsecoyz_(i%dimx     , j%dimy, k%dimz))
		+dx3m*(dy3m*dsecoxyz_( (i-1)%dimx, (j-1)%dimy, k%dimz)+dy3p*dsecoxyz_( (i-1)%dimx, j%dimy, k%dimz))
		+dx3p*(dy3m*dsecoxyz_(i%dimx    , (j-1)%dimy, k%dimz)+dy3p*dsecoxyz_(i%dimx    , j%dimy, k%dimz)));
}


}//end namespace spline
}//end namespace interpolation
}//end namespace numeric

