// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file check_cubic_conv
/// @brief checking new functions for Keys-style/Catmull-Rom cubic interpolation

// libRosetta headers
#include <protocols/viewer/viewers.hh>
#include <devel/init.hh>

#include <numeric/MathNTensor.hh>
#include <numeric/interpolation/polycubic_catmull_rom.hh>

#include <basic/Tracer.hh>


#include <utility/tools/make_vector1.hh>
#include <utility/excn/Exceptions.hh>
#include <numeric/constants.hh>

static basic::Tracer TR( "check_cubic_conv" );

using namespace utility;
using namespace numeric;
using namespace numeric::interpolation;
using numeric::constants::d::pi;

/////////////////////////////////////////////////////////////////////////////////
// Eventually move this into a unit test
/////////////////////////////////////////////////////////////////////////////////
Real
periodic_test_function( Real const & X, Real const & Y)
{
	// F = 3 + sin( 2*pi * ( X + Y ) ) + sin( 6*pi*X ) + cos( 4* pi*Y);
	return 3 + sin( 2 * pi * ( X + Y ) ) + sin( 6 * pi * X ) + cos( 4* pi * Y);
}

Real
cubic_test_function( Real const & X, Real const & Y,
	Real const & Z = 0.0, Real const & W = 0.0)
{
	return 3*X + Y + 3*X*X - X*X*X - 3*Y*Y*Y - 2*X*Y + Y*Z + W*W*W - W*X + 0.5*W;
}

template< typename T, numeric::Size N >
void
test_numerical_deriv(
	MathNTensor< T, N > const & F,
	utility::fixedsizearray1< Real, N > const & x,
	Real const & minval,
	Real const & binwidth,
	CatmullRomSplineBoundaryType const & boundary,
	Real & val,
	utility::fixedsizearray1< T, N > & deriv_analytic )
{
	utility::fixedsizearray1< Real, N > minval_array( minval );
	utility::fixedsizearray1< Real, N > binwidth_array( binwidth );
	utility::fixedsizearray1< CatmullRomSplineBoundaryType, N > boundary_array( boundary );
	TR << std::endl;
	val = polycubic_interpolate_catmull_rom( F, minval_array, binwidth_array, x, boundary_array, deriv_analytic);
	TR << "Value at ";
	for ( Size n = 1; n <= N; n++ ) TR << ' ' << x[n];
	TR << ": "<< val << std::endl;
	TR << "Analytical deriv: ";
	for ( Size n = 1; n <= N; n++ ) TR << ' ' << deriv_analytic[ n ];
	TR << " (boundary: " << to_string( boundary ) << ")" << std::endl;
	TR << std::endl;

	Real delta( 1.0e-8 );
	utility::fixedsizearray1< Real, N > deriv_numeric;
	for ( Size n = 1; n <= N; n++ ) {
		utility::fixedsizearray1< Real, N > x_perturb( x );
		x_perturb[ n ] += delta;
		Real const val_perturb = polycubic_interpolate_catmull_rom( F, minval_array, binwidth_array, x_perturb, boundary_array );
		deriv_numeric[ n ] = ( val_perturb - val ) / delta;
	}
	for ( Size n = 1; n <= N; n++ ) {
		runtime_assert( std::abs( deriv_numeric[ n ] - deriv_analytic[ n ] ) < 1.0e-3 );
	}
}

// @brief tests, based on vals & derivs computed in a prototype MATLAB implementation
void
check_cubic_conv_test()
{
	using utility::tools::make_vector1;

	Real minval( 0.0 );
	Real binwidth( 0.1 );
	Size n_bins( 10 );
	utility::vector1< Real> ctrs( 10, 0 ); // 0, 0.1, ... 0.9
	for ( Size i = 1; i <= n_bins; i++ ) ctrs[ i ] = (i - 1) * binwidth + minval;
	utility::fixedsizearray1< Size, 1 > n_bins_1D( n_bins );
	utility::fixedsizearray1< Size, 2 > n_bins_2D( n_bins );
	utility::fixedsizearray1< Size, 3 > n_bins_3D( n_bins );
	utility::fixedsizearray1< Size, 4 > n_bins_4D( n_bins );
	MathNTensor< Real, 1 > F_1D( n_bins_1D );
	MathNTensor< Real, 2 > F_2D( n_bins_2D );
	MathNTensor< Real, 3 > F_3D( n_bins_3D );
	MathNTensor< Real, 4 > F_4D( n_bins_4D );
	for ( Size i = 0; i < n_bins; i++ ) {
		F_1D( i ) = periodic_test_function( 0.0, ctrs[ i+1 ] );
		for ( Size j = 0; j < n_bins; j++ ) {
			F_2D( i, j ) = periodic_test_function( ctrs[ i+1 ], ctrs[ j+1 ] );
			for ( Size k = 0; k < n_bins; k++ ) {
				F_3D( i, j, k ) = cubic_test_function( ctrs[ i+1 ], ctrs[ j+1 ], ctrs[ k+1 ] );
				for ( Size l = 0; l < n_bins; l++ ) {
					F_4D( i, j, k, l ) = cubic_test_function( ctrs[ i+1 ], ctrs[ j+1 ], ctrs[ k+1 ], ctrs[ l+1 ] );
				}
			}
		}
	}

	Real val( 0.0 );
	{
		TR << "--- Checking 1D periodic ---" << std::endl;
		utility::fixedsizearray1< Real, 1 > x, deriv;
		x[1] = 0.1211;
		test_numerical_deriv( F_1D, x, minval, binwidth, PERIODIC, val, deriv );
		runtime_assert( std::abs( val - 3.767083 ) < 1.0e-3 );
		runtime_assert( std::abs( deriv[1] - -7.692789 ) < 1.0e-3 );
		// val = 3.767083; deriv =  -7.692789

		x[1] = -1.23;
		test_numerical_deriv( F_1D, x, minval, binwidth, PERIODIC, val, deriv );
		runtime_assert( std::abs( val -  1.084389 ) < 1.0e-3 );
		runtime_assert( std::abs( deriv[1] - 2.962611 ) < 1.0e-3 );
		// val = 1.084389; deriv = 2.962611

		test_numerical_deriv( F_1D, x, minval, binwidth, FLAT, val, deriv );
		runtime_assert( std::abs( val - 4.000000 ) < 1.0e-3 );
		runtime_assert( std::abs( deriv[1] -  0.000000 ) < 1.0e-3 );
		// val = 4.000000; deriv = 0.000000

		test_numerical_deriv( F_1D, x, minval, binwidth, LINEAR, val, deriv );
		runtime_assert( std::abs( val - 5.269332 ) < 1.0e-3 );
		runtime_assert( std::abs( deriv[1] -  -1.031978 ) < 1.0e-3 );
		// val = 5.269332; deriv = -1.031978

		x = 1.005;
		test_numerical_deriv( F_1D, x, minval, binwidth, LINEAR, val, deriv );
		runtime_assert( std::abs( val - 4.276602 ) < 1.0e-3 );
		runtime_assert( std::abs( deriv[1] -  14.813053 ) < 1.0e-3 );
		// val = 4.276602; deriv = 14.813053
	}

	{
		TR << "--- Checking 2D periodic ---" << std::endl;
		runtime_assert( std::abs( periodic_test_function( 0.5, 0.5 ) - 4.0 ) < 1.0e-3 );
		runtime_assert( std::abs( F_2D( 5, 5 ) - 4.0 ) < 1.0e-3 );
		utility::fixedsizearray1< Real, 2 > x, deriv;
		x[1] = 0.1211; x[2] = 0.5555;
		test_numerical_deriv( F_2D, x, minval, binwidth, PERIODIC, val, deriv );
		runtime_assert( std::abs( val - 3.594802 ) < 1.0e-3 );
		runtime_assert( std::abs( deriv[1] - -16.701366 ) < 1.0e-3 );
		runtime_assert( std::abs( deriv[2] - -11.475952 ) < 1.0e-3 );
		// -16.701366 -11.475952
	}

	{
		TR << "--- Checking 3D cubic ---" << std::endl;
		utility::fixedsizearray1< Real, 3 > x, deriv;
		x[1] = 0.1211; x[2] = 0.5555; x[3] = 0.4461;
		test_numerical_deriv( F_3D, x, minval, binwidth, PERIODIC, val, deriv );
		runtime_assert( std::abs( val - 0.560022 ) < 1.0e-3 );
		runtime_assert( std::abs( deriv[1] -  2.571593 ) < 1.0e-3 );
		runtime_assert( std::abs( deriv[2] - -1.558867  ) < 1.0e-3 );
		runtime_assert( std::abs( deriv[3] - 0.555500  ) < 1.0e-3 );
		// 2.571593 -1.558867 0.555500

		x[1] = -1.230000; x[2] =  1.500000; x[3] = -0.233333;
		// Value at  -1.230000 1.500000 -0.233333: 3.370585
		// Compare  4.843900 -2.053333 0.500000 (analytic) to
		//          4.843900 -2.053333 0.500000 (numeric) with boundary periodic
		test_numerical_deriv( F_3D, x, minval, binwidth, PERIODIC, val, deriv );
		runtime_assert( std::abs( val - 3.370585 ) < 1.0e-3 );
		runtime_assert( std::abs( deriv[1] -  4.843900 ) < 1.0e-3 );
		runtime_assert( std::abs( deriv[2] - -2.053333  ) < 1.0e-3 );
		runtime_assert( std::abs( deriv[3] -  0.500000  ) < 1.0e-3 );

		// Value at  -1.230000 1.500000 -0.233333: -1.287000
		// Compare  0.000000 0.000000 -0.000000 (analytic) to
		//          0.000000 0.000000 -0.000000 (numeric) with boundary flat
		test_numerical_deriv( F_3D, x, minval, binwidth, FLAT, val, deriv );
		runtime_assert( std::abs( val - -1.287000 ) < 1.0e-3 );
		runtime_assert( std::abs( deriv[1] -  0.000000 ) < 1.0e-3 );
		runtime_assert( std::abs( deriv[2] -  0.000000  ) < 1.0e-3 );
		runtime_assert( std::abs( deriv[3] -  0.000000  ) < 1.0e-3 );

		// Value at  -1.230000 1.500000 -0.233333: -5.299700
		// Compare  0.290000 -3.283333 1.500000 (analytic) to
		//          0.290000 -3.283333 1.500000 (numeric) with boundary linear
		test_numerical_deriv( F_3D, x, minval, binwidth, LINEAR, val, deriv );
		runtime_assert( std::abs( val - -5.299700 ) < 1.0e-3 );
		runtime_assert( std::abs( deriv[1] -  0.290000  ) < 1.0e-3 );
		runtime_assert( std::abs( deriv[2] - -3.283333  ) < 1.0e-3 );
		runtime_assert( std::abs( deriv[3] -  1.500000  ) < 1.0e-3 );

	}


	{
		TR << "--- Checking 4D cubic ---" << std::endl;
		utility::fixedsizearray1< Real, 4 > x, deriv;
		x[1] = 0.1211; x[2] = 0.5555; x[3] = 0.4461; x[4] = 0.0022;
		test_numerical_deriv( F_4D, x, minval, binwidth, PERIODIC, val, deriv );
		runtime_assert( std::abs( val -  0.549209 ) < 1.0e-3 );
		runtime_assert( std::abs( deriv[1] - 2.579914 ) < 1.0e-3 );
		runtime_assert( std::abs( deriv[2] - -1.558867  ) < 1.0e-3 );
		runtime_assert( std::abs( deriv[3] - 0.555500  ) < 1.0e-3 );
		runtime_assert( std::abs( deriv[4] - -4.677011   ) < 1.0e-3 );
		// 2.579914 -1.558867 0.555500 -4.677011
	}
}
///////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	check_cubic_conv_test();
	protocols::viewer::clear_conformation_viewers();
	exit( 0 );
}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
		////////////////////////////////////////////////////////////////////////////
		// setup
		////////////////////////////////////////////////////////////////////////////
		devel::init(argc, argv);

		protocols::viewer::viewer_main( my_main );
	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
	return 0;

}
