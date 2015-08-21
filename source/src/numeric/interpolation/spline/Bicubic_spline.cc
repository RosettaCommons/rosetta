// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
///
/////////////////////////////////////////////////////////////////////////

// Unit headers
#include <numeric/interpolation/spline/Bicubic_spline.hh>

// Package headers
#include <numeric/types.hh>
#include <numeric/interpolation/spline/Cubic_spline.hh>
#include <numeric/MathVector_operations.hh>
#include <numeric/MathMatrix.hh>

// C++ headers
#include <iostream>

namespace numeric {
namespace interpolation {
namespace spline {


////////////////
// operations //
////////////////


// BORDER determines the behavior of the spline at the borders (natural, first derivative, periodic)
// START the start of the interval the spline is defined on
// DELTA the distance between two support points of the spline
// RESULTS the function values at the support points of the spline
// FIRSTBE values for the first order derivative at begin and end of spline (only FIRSTDER)
// train BicubicSpline
void BicubicSpline::train(
	const BorderFlag BORDER[2],
	const Real START[2],
	const Real DELTA[2],
	const MathMatrix< Real> &RESULTS,
	const bool LINCONT[2],
	const std::pair< Real, Real> FIRSTBE[2]
)
{
	//check, if the points are given in positive direction


	//determine values for all dimensions
	const Size dimx( RESULTS.get_number_rows());
	const Size dimy( RESULTS.get_number_cols());

	//assigning values
	border_[ 0] = BORDER[ 0];
	border_[ 1] = BORDER[ 1];
	start_[ 0]  = START[ 0];
	start_[ 1]  = START[ 1];
	delta_[ 0]  = DELTA[ 0];
	delta_[ 1]  = DELTA[ 1];

	values_ = RESULTS;
	dsecox_ = RESULTS;
	dsecoy_ = RESULTS;
	dsecoxy_ = RESULTS;
	LinCont_[ 0] = LINCONT[ 0];
	LinCont_[ 1] = LINCONT[ 1];
	firstbe_[ 0] = FIRSTBE[ 0];
	firstbe_[ 1] = FIRSTBE[ 1];

	/*for ( Size jj = 0; jj < 36; ++jj ) {
	for ( Size kk = 0; kk < 36; ++kk ) {
	std::cout << values_[ jj ][ kk ] << " ";
	}
	std::cout << std::endl;
	}*/

	//train three times for fxx, fyy, fxxyy
	//reduction to Spline by training only rows/columns at the same time
	for ( Size row( 0); row < dimx; ++row ) {
		CubicSpline cs;
		//std::cout << "Doing row " << row << " for dsecoy" << std::endl << "Old row is ";
		//for ( Size i = 0; i < dimx; ++i) std::cout <<RESULTS.get_row( row)(i) << " ";
		cs.train( BORDER[ 1], START[ 1], DELTA[ 1], RESULTS.get_row( row), FIRSTBE[ 1]);
		dsecoy_.replace_row( row, cs.get_dsecox());
		//std::cout << std::endl << " and new row is ";
		//for ( Size i = 0; i < dimx; ++i) std::cout <<dsecoy_.get_row( row)(i) << " ";
		//std::cout << std::endl;
	}

	for ( Size col( 0); col < dimy; ++col ) {
		CubicSpline cs;
		//std::cout << "Doing col " << col << " for dsecox" << std::endl << "Old col is ";
		//for ( Size i = 0; i < dimy; ++i) std::cout <<RESULTS.get_col( col)(i) << " ";
		cs.train( BORDER[ 0], START[ 0], DELTA[ 0], RESULTS.get_col( col), FIRSTBE[ 0]);
		dsecox_.replace_col( col, cs.get_dsecox());
		//std::cout << std::endl << " and new col is ";
		//for ( Size i = 0; i < dimy; ++i) std::cout <<dsecox_.get_col( col)(i) << " ";
		//std::cout << std::endl;
	}

	for ( Size row( 0); row < dimx; ++row ) {
		CubicSpline cs;
		//std::cout << "Doing row " << row << " for dsecoxy" << std::endl << "Old row is ";
		//for ( Size i = 0; i < dimx; ++i) std::cout <<dsecox_.get_row( row)(i) << " ";
		cs.train( BORDER[ 1], START[ 1], DELTA[ 1], dsecox_.get_row( row), FIRSTBE[ 1]);
		dsecoxy_.replace_row( row, cs.get_dsecox());
		//std::cout << std::endl << " and new row is ";
		//for ( Size i = 0; i < dimx; ++i) std::cout <<dsecoxy_.get_row( row)(i) << " ";
		//std::cout << std::endl;
	}
	return;
}


/// @return value at certain (x, y)
Real BicubicSpline::F( const MathVector< Real> &ARGUMENTS) const
{
	const Real x( ARGUMENTS( 0));
	const Real y( ARGUMENTS( 1));

	// check that there are two argument values given
	return F( x, y );
}

Real BicubicSpline::F( Real x, Real y ) const
{
	const int dimx( values_.get_number_rows());
	const int dimy( values_.get_number_cols());

	// check if argument x is in range for non-periodic splines
	if ( x < start_[ 0] ) {
		switch( border_[ 0])
				{
				case e_FirstDer : //BCL_Assert( LinCont_[ 0], "argument out of range for non-periodic spline!");
					//BCL_Message( basic::Message::e_Debug, "argument out of range, using linear continuation");
					return F( MakeVector( start_[ 0], y))+( x-start_[ 0] )*firstbe_[ 0].first;

				case e_Natural : //BCL_Assert( LinCont_[ 0], "argument out of range for non-periodic spline!");
					//BCL_Message( basic::Message::e_Debug, "argument out of range, using linear continuation");
					return F( MakeVector( start_[ 0], y))+( x-start_[ 0] )*dFdx( MakeVector( start_[ 0], y));

				case e_Periodic : break;
				}
	}

	if ( start_[ 0] + ( dimx - 1) * delta_[ 0] < x ) {
		switch( border_[ 0])
				{
				case e_FirstDer : //BCL_Assert( LinCont_[ 0], "argument out of range for non-periodic spline!");
					//BCL_Message( basic::Message::e_Debug, "argument out of range, using linear continuation");
					return F( MakeVector( start_[ 0] + ( dimx-1 ) * delta_[ 0] , y ))+( x - start_[ 0] - ( dimx - 1) * delta_[ 0]) * firstbe_[ 0].second;

				case e_Natural : //BCL_Assert( LinCont_[ 0], "argument out of range for non-periodic spline!");
					//BCL_Message( basic::Message::e_Debug, "argument out of range, using linear continuation");
					return F( MakeVector( start_[ 0] + ( dimx-1 ) * delta_[ 0] , y ))+( x - start_[ 0] - ( dimx - 1) * delta_[ 0]) * dFdx( MakeVector( start_[ 0] + ( dimx-1 ) * delta_[ 0], y));

				case e_Periodic : break;
				}
	}

	//check if argument y is in range for non-periodic splines
	if ( y  < start_[ 1] ) {
		switch( border_[ 1])
				{
				case e_FirstDer : //BCL_Assert( LinCont_[ 1], "argument out of range for non-periodic spline!");
					//BCL_Message( basic::Message::e_Debug, "argument out of range, using linear continuation");
					return F( MakeVector( x, start_[ 1]))+(y-start_[ 1]) * firstbe_[ 1].first;

				case e_Natural : //BCL_Assert( LinCont_[ 1], "argument out of range for non-periodic spline!");
					// BCL_Message( basic::Message::e_Debug, "argument out of range, using linear continuation");
					return F( MakeVector( x, start_[ 1]))+(y-start_[ 1]) * dFdy( MakeVector( x, start_[ 1]));

				case e_Periodic : break;
				}
	}

	if ( start_[ 1] + ( dimy-1 ) * delta_[ 1] <  y ) {
		switch( border_[ 1])
				{
				case e_FirstDer : //BCL_Assert( LinCont_[ 1], "argument out of range for non-periodic spline!");
					//BCL_Message( basic::Message::e_Debug, "argument out of range, using linear continuation");
					return F( MakeVector( x, start_[ 1] + ( dimy-1 ) * delta_[ 1])) + ( y - start_[ 1] - ( dimy - 1) * delta_[ 1])*firstbe_[ 1].second;

				case e_Natural : //BCL_Assert( LinCont_[ 1], "argument out of range for non-periodic spline!");
					//BCL_Message( basic::Message::e_Debug, "argument out of range, using linear continuation");
					return F( MakeVector( x, start_[ 1] + ( dimy-1 ) * delta_[ 1])) + ( y - start_[ 1] - ( dimy - 1) * delta_[ 1])*dFdy( MakeVector( x, start_[ 1] + ( dimy-1 ) * delta_[ 1]));

				case e_Periodic : break;
				}
	}

	//determine i with start_[ 0]+(i-1)*delta_[ 0] < x < start_[ 0]+i*delta_[ 0] for the correct supporting points
	int i( int (floor( (x-start_[ 0])/delta_[ 0])+1));

	//determine j with start_[ 1]+(j-1)*delta_[ 1] < y < start_[ 1]+j*delta_[ 1] for the correct supporting points
	int j( int (floor( (y-start_[ 1])/delta_[ 1])+1));

	const Real dxp( ( x-start_[ 0])/delta_[ 0] - floor( ( x-start_[ 0]) / delta_[ 0]));
	const Real dxm( 1 - dxp);
	const Real dx3p( ( dxp*dxp*dxp - dxp) * sqr( delta_[ 0]) / 6); // =0 at the grid points, adds cubic part of the spline
	const Real dx3m( ( dxm*dxm*dxm - dxm) * sqr( delta_[ 0]) / 6); // =0 at the grid points, adds cubic part of the spline

	const Real dyp( ( y-start_[ 1])/delta_[ 1] - floor( ( y-start_[ 1]) / delta_[ 1]));
	const Real dym( ( 1 - dyp));
	const Real dy3p( ( dyp*dyp*dyp - dyp) * sqr( delta_[ 1]) / 6); // =0 at the grid points, adds cubic part of the spline
	const Real dy3m( ( dym*dym*dym - dym) * sqr( delta_[ 1]) / 6); // =0 at the grid points, adds cubic part of the spline

	//generate positive values to prevent some problems with the indices
	while ( i < 1 ) i += dimx;
	while ( j < 1 ) j += dimy;

	return
		dxm * ( dym * values_( (   i - 1) % dimx, ( j - 1) % dimy) + dyp  * values_( (  i - 1) % dimx, j % dimy))
		+ dxp * ( dym * values_(    i      % dimx, ( j - 1) % dimy) + dyp  * values_(   i      % dimx, j % dimy))
		+dx3m * ( dym * dsecox_( (   i - 1) % dimx, ( j - 1) % dimy) + dyp  * dsecox_( (  i - 1) % dimx, j % dimy))
		+dx3p * ( dym * dsecox_(    i      % dimx, ( j - 1) % dimy) + dyp  * dsecox_(   i      % dimx, j % dimy))
		+ dxm * ( dy3m * dsecoy_( (  i - 1) % dimx, ( j - 1) % dimy) + dy3p * dsecoy_( (  i - 1) % dimx, j % dimy))
		+ dxp * ( dy3m * dsecoy_(   i      % dimx, ( j - 1) % dimy) + dy3p * dsecoy_(   i      % dimx, j % dimy))
		+dx3m * ( dy3m * dsecoxy_( ( i - 1) % dimx, ( j - 1) % dimy) + dy3p * dsecoxy_( ( i - 1) % dimx, j % dimy))
		+dx3p * ( dy3m * dsecoxy_(  i      % dimx, ( j - 1) % dimy) + dy3p * dsecoxy_(  i      % dimx, j % dimy));
}


/// @return partial derivative at certain (x, y) for x
Real BicubicSpline::dFdx( const MathVector< Real> &ARGUMENTS) const
{
	//BCL_Assert( ARGUMENTS.GetSize() == 2, "number of arguments doesn't match");

	const Real x( ARGUMENTS(0));
	const Real y( ARGUMENTS(1));
	return dFdx( x, y );
}

/// @return partial derivative at certain (x, y) for x
Real BicubicSpline::dFdx( Real x, Real y ) const
{

	const int dimx( values_.get_number_rows());
	const int dimy( values_.get_number_cols());

	//check if argument x is in range for non-periodic splines
	if ( x < start_[ 0] ) {
		switch( border_[ 0] )
				{
				case e_FirstDer : //BCL_Assert( LinCont_[ 0], "argument out of range for non-periodic spline!");BCL_Message( basic::Message::e_Verbose, "argument out of range, using linear continuation");return firstbe_[ 0].first;
				case e_Natural : //BCL_Assert( LinCont_[ 0], "argument out of range for non-periodic spline!");BCL_Message( basic::Message::e_Verbose, "argument out of range, using linear continuation");return dFdx( MakeVector( start_[ 0], y));
				case e_Periodic : break;
				}
	}

	if ( start_[ 0] + ( dimx-1 ) * delta_[ 0] < x ) {
		switch( border_[ 0] )
				{
				case e_FirstDer : //BCL_Assert( LinCont_[ 0], "argument out of range for non-periodic spline!");BCL_Message( basic::Message::e_Verbose, "argument out of range, using linear continuation");return firstbe_[ 0].second;
				case e_Natural : //BCL_Assert( LinCont_[ 0], "argument out of range for non-periodic spline!");BCL_Message( basic::Message::e_Verbose, "argument out of range, using linear continuation");return dFdx( MakeVector( start_[ 0] + ( dimx-1 ) * delta_[ 0], y));
				case e_Periodic : break;
				}
	}

	//check if argument y is in range for non-periodic splines
	if ( y < start_[ 1] ) {
		switch( border_[ 1] )
				{
				case e_FirstDer : //BCL_Assert( LinCont_[ 1], "argument out of range for non-periodic spline!");BCL_Message( basic::Message::e_Verbose, "argument out of range, using linear continuation");return dFdx( MakeVector( x, start_[ 1]));
				case e_Natural : //BCL_Assert( LinCont_[ 1], "argument out of range for non-periodic spline!");BCL_Message( basic::Message::e_Verbose, "argument out of range, using linear continuation");return dFdx( MakeVector( x, start_[ 1]));
				case e_Periodic : break;
				}
	}

	if ( start_[ 1] + ( dimy-1 ) * delta_[ 1] < y ) {
		switch( border_[ 1] )
				{
				case e_FirstDer : //BCL_Assert( LinCont_[ 0], "argument out of range for non-periodic spline!");BCL_Message( basic::Message::e_Verbose, "argument out of range, using linear continuation");return dFdx( MakeVector( x, start_[ 1] + ( dimy-1 ) * delta_[ 1]));
				case e_Natural : //BCL_Assert( LinCont_[ 0], "argument out of range for non-periodic spline!");BCL_Message( basic::Message::e_Verbose, "argument out of range, using linear continuation");return dFdx( MakeVector( x, start_[ 1] + ( dimy-1 ) * delta_[ 1]));
				case e_Periodic : break;
				}
	}

	// determine i with start_[ 0]+(i-1)*delta_[ 0] < x < start_[ 0]+i*delta_[ 0] for the correct supporting points
	int i( int( floor( ( x - start_[ 0]) / delta_[ 0])));
	while ( start_[ 0] + i * delta_[ 0] < x ) { i++;}

	// determine j with start_[ 1]+(j-1)*delta_[ 1] < y < start_[ 1]+j*delta_[ 1] for the correct supporting points
	int j( int( floor( ( y - start_[ 1]) / delta_[ 1])));
	while ( start_[ 1] + j * delta_[ 1] < y ) { j++;}

	//see F(x, y) for a short explanation of the values
	const Real delta_aktx( x-start_[ 0] - ( i - 1) * delta_[ 0]);
	const Real delta_akty( y-start_[ 1] - ( j - 1) * delta_[ 1]);

	const Real dxp( delta_aktx / delta_[ 0]);
	const Real dxm( 1 - dxp);

	const Real dyp( delta_akty / delta_[ 1]);
	const Real dym( 1 - dyp);
	const Real dy3p( ( dyp * dyp * dyp - dyp) * sqr( delta_[ 1]) / 6);
	const Real dy3m( ( dym * dym * dym - dym) * sqr( delta_[ 1]) / 6);

	//generate positive values to prevent some problems with the indices
	while ( i < 1 ) { i += dimx;}
	while ( j < 1 ) { j += dimy;}
	//  BCL_Message( basic::Message::e_Critical, "Final i orig: " + basic::Format()( i));

	return
		-( dym * values_( ( i - 1) % dimx, ( j - 1) % dimy) + dyp * values_( ( i - 1) % dimx  , j % dimy)) / delta_[ 0]
		+( dym * values_( i % dimx      , ( j - 1) % dimy) + dyp * values_( i % dimx        , j % dimy)) / delta_[ 0]
		- ( 3 * dxm * dxm - 1) * delta_[ 0] / 6 *( dym*dsecox_( ( i - 1) % dimx, ( j - 1) % dimy) + dyp * dsecox_( ( i - 1) % dimx, j%dimy))
		+ ( 3 * dxp * dxp - 1) * delta_[ 0] / 6 *( dym*dsecox_( i % dimx      , ( j - 1) % dimy) + dyp * dsecox_( i % dimx      , j%dimy))
		-( dy3m * dsecoy_( ( i - 1) % dimx , ( j - 1) % dimy) + dy3p * dsecoy_( ( i-1) % dimx , j % dimy)) / delta_[ 0]
		+( dy3m * dsecoy_( i % dimx       , ( j - 1) % dimy) + dy3p * dsecoy_( i % dimx     , j % dimy)) / delta_[ 0]
		- ( 3 * dxm * dxm - 1) * delta_[ 0] / 6 *( dy3m*dsecoxy_( ( i - 1)%dimx, ( j - 1) % dimy) + dy3p * dsecoxy_( ( i - 1) % dimx, j % dimy))
		+ ( 3 * dxp * dxp - 1) * delta_[ 0] / 6 *( dy3m*dsecoxy_( i % dimx    , ( j - 1) % dimy) + dy3p * dsecoxy_( i % dimx      , j % dimy));
}


/// @return partial derivative at certain (x, y) for y
Real BicubicSpline::dFdy( const MathVector< Real> &ARGUMENTS) const
{
	// BCL_Assert( ARGUMENTS.GetSize() == 2, "number of arguments doesn't match");

	const Real x( ARGUMENTS( 0));
	const Real y( ARGUMENTS( 1));
	return dFdy( x, y );
}

/// @return partial derivative at certain (x, y) for y
Real BicubicSpline::dFdy( Real const x, Real const y) const
{

	const int dimx( values_.get_number_rows());
	const int dimy( values_.get_number_cols());

	//check if argument x is in range for non-periodic splines
	if ( x < start_[ 0] ) {
		switch( border_[ 0])
				{
				case e_FirstDer : //BCL_Assert( LinCont_[ 0], "argument out of range for non-periodic spline!");BCL_Message( basic::Message::e_Verbose, "argument out of range, using linear continuation");return dFdy( MakeVector( start_[ 0], y));
				case e_Natural : //BCL_Assert( LinCont_[ 0], "argument out of range for non-periodic spline!");BCL_Message( basic::Message::e_Verbose, "argument out of range, using linear continuation");return dFdy( MakeVector( start_[ 0], y));
				case e_Periodic : break;
				}
	}

	if ( start_[ 0] + ( dimx - 1) * delta_[ 0] < x ) {
		switch( border_[ 0])
				{
				case e_FirstDer : //BCL_Assert( LinCont_[ 0], "argument out of range for non-periodic spline!");BCL_Message( basic::Message::e_Verbose, "argument out of range, using linear continuation");return dFdy( MakeVector( start_[ 0] + ( dimx-1 ) * delta_[ 0], y));
				case e_Natural : //BCL_Assert( LinCont_[ 0], "argument out of range for non-periodic spline!");BCL_Message( basic::Message::e_Verbose, "argument out of range, using linear continuation");return dFdy( MakeVector( start_[ 0] + ( dimx-1 ) * delta_[ 0], y));
				case e_Periodic : break;
				}
	}

	//check if argument y is in range for non-periodic splines
	if ( y < start_[ 1] ) {
		switch( border_[ 1])
				{
				case e_FirstDer : //BCL_Assert( LinCont_[ 1], "argument out of range for non-periodic spline!");BCL_Message( basic::Message::e_Verbose, "argument out of range, using linear continuation");return firstbe_[ 1].first;
				case e_Natural : //BCL_Assert( LinCont_[ 1], "argument out of range for non-periodic spline!");BCL_Message( basic::Message::e_Verbose, "argument out of range, using linear continuation");return dFdy( MakeVector( x, start_[ 1]));
				case e_Periodic : break;
				}
	}

	if ( start_[ 1] + ( dimy - 1) * delta_[ 1] < y ) {
		switch( border_[ 1])
				{
				case e_FirstDer : //BCL_Assert( LinCont_[ 0], "argument out of range for non-periodic spline!");BCL_Message( basic::Message::e_Verbose, "argument out of range, using linear continuation");return firstbe_[ 1].second;
				case e_Natural : //BCL_Assert( LinCont_[ 0], "argument out of range for non-periodic spline!");BCL_Message( basic::Message::e_Verbose, "argument out of range, using linear continuation");return dFdy( MakeVector( x, start_[ 1] + ( dimy-1 ) * delta_[ 1]));
				case e_Periodic : break;
				}
	}

	//determine i with start_[ 0]+(i-1)*delta_[ 0] < x < start_[ 0]+i*delta_[ 0] for the correct supporting points
	int i( int( floor( ( x - start_[ 0])/delta_[ 0])));
	while ( start_[ 0] + i * delta_[ 0] < x ) { i++;}
	if ( !i ) {
		while ( start_[ 0] + i * delta_[ 0] > x ) { i--;}
		i++;
	}

	//determine j with start_[ 1]+(j-1)*delta_[ 1] < y < start_[ 1]+j*delta_[ 1] for the correct supporting points
	int j( int( floor( ( y - start_[ 1])/delta_[ 1])));
	while ( start_[ 1] + j * delta_[ 1] < y ) {j++;}
	if ( !j ) {
		while  ( start_[ 1]+j*delta_[ 1]>y ) {j--;}
		j++;
	}

	//see F(x, y) for a short explanation of the values
	const Real delta_aktx( x - start_[ 0] - ( i - 1) * delta_[ 0]);
	const Real delta_akty( y - start_[ 1] - ( j - 1) * delta_[ 1]);

	const Real dxp( delta_aktx / delta_[ 0]);
	const Real dxm( 1 - dxp);
	const Real dx3p( ( dxp * dxp * dxp - dxp) * sqr( delta_[ 0]) / 6);
	const Real dx3m( ( dxm * dxm * dxm - dxm) * sqr( delta_[ 0]) / 6);

	const Real dyp( delta_akty / delta_[ 1]);
	const Real dym( 1 - dyp);

	//generate positive values to prevent some problems with the indices
	while ( i < 1 ) { i += dimx;}
	while ( j < 1 ) { j += dimy;}

	return
		dxm *( -values_( ( i-1)%dimx  , (j-1)%dimy)+values_( (i-1)%dimx  , j%dimy))/delta_[ 1]
		+ dxp *( -values_( i%dimx      , (j-1)%dimy)+values_(i%dimx      , j%dimy))/delta_[ 1]
		+dx3m *( -dsecox_( ( i-1)%dimx  , (j-1)%dimy)+dsecox_( (i-1)%dimx  , j%dimy))/delta_[ 1]
		+dx3p *( -dsecox_( i%dimx      , (j-1)%dimy)+dsecox_(i%dimx      , j%dimy))/delta_[ 1]
		+ dxm *( -( 3 * dym * dym - 1) * dsecoy_( ( i-1)%dimx , ( j - 1)% dimy) +( 3 * dyp * dyp - 1) * dsecoy_( ( i-1)%dimx , j % dimy)) * delta_[ 1]/ 6
		+ dxp *( -( 3 * dym * dym - 1) * dsecoy_( i%dimx     , ( j - 1)% dimy) +( 3 * dyp * dyp - 1) * dsecoy_( i%dimx     , j % dimy)) * delta_[ 1]/ 6
		+dx3m *( -( 3 * dym * dym - 1) * dsecoxy_( ( i-1)%dimx, ( j - 1)% dimy) +( 3 * dyp * dyp - 1) * dsecoxy_( ( i-1)%dimx, j % dimy)) * delta_[ 1]/ 6
		+dx3p *( -( 3 * dym * dym - 1) * dsecoxy_( i%dimx    , ( j - 1)% dimy) +( 3 * dyp * dyp - 1) * dsecoxy_( i%dimx    , j % dimy)) * delta_[ 1]/ 6;
}


/// @return value and derivative at certain (x, y)
std::pair<Real, MathVector<Real> > BicubicSpline::FdF( const MathVector< Real> &ARGUMENTS) const
{
	//BCL_Assert( ARGUMENTS.GetSize() == 2, "number of arguments doesn't match");

	Real x = ARGUMENTS(0);
	Real y = ARGUMENTS(1);

	int dimx = values_.get_number_rows();
	int dimy = values_.get_number_cols();

	//auxiliary variables for the function value and the derivatives
	Real fvalue( 0), dfdxvalue( 0), dfdyvalue( 0);

	//check if argument is in range for non-periodic splines
	if ( ( ( border_[ 0] != e_Periodic ) && ( x < start_[ 0] || start_[ 0] + ( dimx-1 ) * delta_[ 0] < x))
			|| ( ( border_[ 1] != e_Periodic ) && (y < start_[ 1] || start_[ 1] + ( dimy-1 ) * delta_[ 1] < y)) ) {
		if ( x < start_[ 0] ) {
			//BCL_Assert( LinCont_[ 0], "argument out of range for non-periodic spline!");
			//BCL_Message( basic::Message::e_Debug, "argument out of range, using linear continuation");
			fvalue    = F(MakeVector( start_[ 0], y))+(x-start_[ 0])*dFdx( MakeVector(start_[ 0], y));
			dfdxvalue = dFdx( MakeVector( start_[ 0], y));
			dfdyvalue = dFdy( MakeVector( start_[ 0], y));
		}
		if ( x > start_[ 0] + ( dimx-1 ) * delta_[ 0] ) {
			//BCL_Assert( LinCont_[ 0], "argument out of range for non-periodic spline!");
			//BCL_Message( basic::Message::e_Debug, "argument out of range, using linear continuation");
			fvalue    = F(    MakeVector( start_[ 0] + ( dimx-1 ) * delta_[ 0] , y))+(x-start_[ 0] - ( dimx-1 ) * delta_[ 0])*dFdx( MakeVector( start_[ 0] + ( dimx-1 ) * delta_[ 0], y));
			dfdxvalue = dFdx( MakeVector( start_[ 0] + ( dimx-1 ) * delta_[ 0] , y));
			dfdyvalue = dFdy( MakeVector( start_[ 0] + ( dimx-1 ) * delta_[ 0] , y));
		}
		if ( y < start_[ 1] ) {
			//BCL_Assert( LinCont_[ 1], "argument out of range for non-periodic spline!");
			//BCL_Message( basic::Message::e_Debug, "argument out of range, using linear continuation");
			fvalue    = F(    MakeVector( x, start_[ 1]))+(y-start_[ 1])*dFdy( MakeVector( x, start_[ 1]));
			dfdxvalue = dFdx( MakeVector( x, start_[ 1]));
			dfdyvalue = dFdy( MakeVector( x, start_[ 1]));
		}
		if ( y > start_[ 1] + ( dimy-1 ) * delta_[ 1] ) {
			//BCL_Assert( LinCont_[ 1], "argument out of range for non-periodic spline!");
			//BCL_Message( basic::Message::e_Debug, "argument out of range, using linear continuation");
			fvalue    = F(    MakeVector( x, start_[ 1] + ( dimy-1 ) * delta_[ 1]))+(y-start_[ 1] - ( dimy-1 ) * delta_[ 1])*dFdy( MakeVector( x, start_[ 1] + ( dimy-1 ) * delta_[ 1]));
			dfdxvalue = dFdx( MakeVector( x, start_[ 1] + ( dimy-1 ) * delta_[ 1]));
			dfdyvalue = dFdy( MakeVector( x, start_[ 1] + ( dimy-1 ) * delta_[ 1]));
		}
	} else {
		//determine i with start_[ 0]+(i-1)*delta_[ 0] < x < start_[ 0]+i*delta_[ 0] for the correct supporting points
		int    i(int (floor( (x-start_[ 0])/delta_[ 0])));
		while  ( start_[ 0]+i*delta_[ 0]<x ) i++;
		if ( !i ) {
			while  ( start_[ 0]+i*delta_[ 0]>x ) i--;
			i++;
		}

		//determine j with start_[ 1]+(j-1)*delta_[ 1] < y < start_[ 1]+j*delta_[ 1] for the correct supporting points
		int    j(int (floor( (y-start_[ 1])/delta_[ 1])));
		while  ( start_[ 1]+j*delta_[ 1]<y ) j++;
		if ( !j ) {
			while  ( start_[ 1]+j*delta_[ 1]>y ) j--;
			j++;
		}

		//see method F(x,y) for detailed formula

		Real delta_aktx = x-start_[ 0]-(i-1)*delta_[ 0];
		Real delta_akty = y-start_[ 1]-(j-1)*delta_[ 1];

		Real dxp(delta_aktx/delta_[ 0]);
		Real dxm( 1 - dxp);
		Real dx3p( ( dxp*dxp*dxp - dxp) * sqr( delta_[ 0]) / 6);
		Real dx3m( ( dxm*dxm*dxm - dxm) * sqr( delta_[ 0]) / 6);

		Real dyp(delta_akty/delta_[ 1]);
		Real dym( 1 - dyp);
		Real dy3p( ( dyp*dyp*dyp - dyp) * sqr( delta_[ 1]) / 6);
		Real dy3m( ( dym*dym*dym - dym) * sqr( delta_[ 1]) / 6);

		fvalue =
			dxm*(dym*values_( (i-1)%dimx  , (j-1)%dimy)+dyp*values_( (i-1)%dimx  , j%dimy))
			+ dxp*(dym*values_(i%dimx      , (j-1)%dimy)+dyp*values_(i%dimx      , j%dimy))
			+dx3m*(dym*dsecox_( (i-1)%dimx  , (j-1)%dimy)+dyp*dsecox_( (i-1)%dimx  , j%dimy))
			+dx3p*(dym*dsecox_(i%dimx      , (j-1)%dimy)+dyp*dsecox_(i%dimx      , j%dimy))
			+ dxm*(dy3m*dsecoy_( (i-1)%dimx , (j-1)%dimy)+dy3p*dsecoy_( (i-1)%dimx , j%dimy))
			+ dxp*(dy3m*dsecoy_(i%dimx     , (j-1)%dimy)+dy3p*dsecoy_(i%dimx     , j%dimy))
			+dx3m*(dy3m*dsecoxy_( (i-1)%dimx, (j-1)%dimy)+dy3p*dsecoxy_( (i-1)%dimx, j%dimy))
			+dx3p*(dy3m*dsecoxy_(i%dimx    , (j-1)%dimy)+dy3p*dsecoxy_(i%dimx    , j%dimy))
			;

		dfdxvalue =
			-(dym*values_( (i-1)%dimx  , (j-1)%dimy)+dyp*values_( (i-1)%dimx  , j%dimy))/delta_[ 0]
			+(dym*values_(i%dimx      , (j-1)%dimy)+dyp*values_(i%dimx      , j%dimy))/delta_[ 0]
			- (3 * dxm*dxm - 1) * delta_[ 0] / 6*(dym*dsecox_( (i-1)%dimx  , (j-1)%dimy)+dyp*dsecox_( (i-1)%dimx  , j%dimy))
			+ (3 * dxp*dxp - 1) * delta_[ 0] / 6*(dym*dsecox_(i%dimx      , (j-1)%dimy)+dyp*dsecox_(i%dimx      , j%dimy))
			-(dy3m*dsecoy_( (i-1)%dimx , (j-1)%dimy)+dy3p*dsecoy_( (i-1)%dimx , j%dimy))/delta_[ 0]
			+(dy3m*dsecoy_(i%dimx     , (j-1)%dimy)+dy3p*dsecoy_(i%dimx     , j%dimy))/delta_[ 0]
			- (3 * dxm*dxm - 1) * delta_[ 0] / 6*(dy3m*dsecoxy_( (i-1)%dimx, (j-1)%dimy)+dy3p*dsecoxy_( (i-1)%dimx, j%dimy))
			+ (3 * dxp*dxp - 1) * delta_[ 0] / 6*(dy3m*dsecoxy_(i%dimx    , (j-1)%dimy)+dy3p*dsecoxy_(i%dimx    , j%dimy))
			;

		dfdyvalue =
			dxm*(-values_( (i-1)%dimx  , (j-1)%dimy)+values_( (i-1)%dimx  , j%dimy))/delta_[ 1]
			+ dxp*(-values_(i%dimx      , (j-1)%dimy)+values_(i%dimx      , j%dimy))/delta_[ 1]
			+dx3m*(-dsecox_( (i-1)%dimx  , (j-1)%dimy)+dsecox_( (i-1)%dimx  , j%dimy))/delta_[ 1]
			+dx3p*(-dsecox_(i%dimx      , (j-1)%dimy)+dsecox_(i%dimx      , j%dimy))/delta_[ 1]
			+ dxm*(-(3*dym*dym-1)*dsecoy_( (i-1)%dimx , (j-1)%dimy)+(3*dyp*dyp-1)*dsecoy_( (i-1)%dimx , j%dimy))* delta_[ 1]/ 6
			+ dxp*(-(3*dym*dym-1)*dsecoy_(i%dimx     , (j-1)%dimy)+(3*dyp*dyp-1)*dsecoy_(i%dimx     , j%dimy))* delta_[ 1]/ 6
			+dx3m*(-(3*dym*dym-1)*dsecoxy_( (i-1)%dimx, (j-1)%dimy)+(3*dyp*dyp-1)*dsecoxy_( (i-1)%dimx, j%dimy))* delta_[ 1]/ 6
			+dx3p*(-(3*dym*dym-1)*dsecoxy_(i%dimx    , (j-1)%dimy)+(3*dyp*dyp-1)*dsecoxy_(i%dimx    , j%dimy))* delta_[ 1]/ 6
			;
	}

	Real dfvalues[] = { dfdxvalue, dfdyvalue};
	MathVector<Real> dfvector( 2, dfvalues);

	return std::pair< Real, MathVector< Real> >( fvalue, dfvector);
}


}//end namespace spline
}//end namespace interpolation
}//end namespace numeric

