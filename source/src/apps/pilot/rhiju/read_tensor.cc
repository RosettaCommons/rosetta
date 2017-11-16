// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file read_tensor.cc
/// @brief checking if we can read in a tensor from a .bin.gz and associated .json

// libRosetta headers
#include <protocols/viewer/viewers.hh>
#include <devel/init.hh>

#include <numeric/MathNTensor_io.hh>
#include <numeric/interpolation/spline/PolycubicSpline.hh>
#include <numeric/interpolation/spline/PolycubicSpline.tmpl.hh>
#include <numeric/interpolation/interpolation.hh>

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/database/open.hh>

#include <core/types.hh>
#include <core/scoring/loop_graph/evaluator/util.hh>
#include <core/scoring/loop_graph/evaluator/SixDTransRotPotential.hh>

#include <utility/io/izstream.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/json_spirit/json_spirit.h>
#include <utility/json_spirit/json_spirit_tools.hh>
#include <utility/excn/Exceptions.hh>
#include <numeric/xyz.io.hh>

static basic::Tracer TR( "read_tensor" );

using namespace basic::options::OptionKeys;
using namespace basic::options;
using namespace utility;
using namespace numeric;
using namespace numeric::interpolation;
using namespace numeric::interpolation::spline;

OPT_KEY( String, tensor_file )


/////////////////////////////////////////////////////////////////////////////////
// Eventually move this into a unit test
/////////////////////////////////////////////////////////////////////////////////
void
check_tensor( numeric::MathNTensor< double, 6 > const & T,
	utility::vector1< numeric::Size > const & checkbins );

void
check_multilinear_interpolation( numeric::MathNTensor< double, 6 > const & T,
	utility::fixedsizearray1< numeric::Real, 6 > const & minval,
	utility::fixedsizearray1< numeric::Real, 6 > const & binwidth,
	utility::fixedsizearray1< numeric::Real, 6 > const & xs );

void
check_spline( numeric::MathNTensor< double, 6 > const & T,
	utility::fixedsizearray1< numeric::Real, 6 > const & minval,
	utility::fixedsizearray1< numeric::Real, 6 > const & binwidth,
	utility::fixedsizearray1< numeric::Real, 6 > const & xs );

/////////////////////////////////
void
check_6D_potential( numeric::MathNTensor< double, 6 > const & T,
	utility::json_spirit::mObject const & json,
	utility::fixedsizearray1< numeric::Real, 6 > const & xs )
{
	using core::Vector;
	using namespace numeric;
	core::scoring::loop_graph::evaluator::SixDTransRotPotential potential( T, json );
	Vector const t( xs[1], xs[2], xs[3] );
	Vector v( xs[4], xs[5], xs[6] );
	Vector v_norm = v.normalized();
	std::pair< Vector, Vector > deriv;
	utility::vector1< Real > xvals;
	for ( Real x =  120.0; x <= 180.0; x += 10.0 ) xvals.push_back( x );
	for ( Real x = -180.0; x <= -120.0; x += 10.0 ) xvals.push_back( x );

	TR << "Rotation axis: " << v_norm << std::endl;;

	TR << "Not enforcing continuity at pi" << std::endl;
	potential.set_enforce_continuity_at_pi( false );
	for ( Real const & x : xvals ) {
		v = x * v_norm;
		TR << "evaluating 6D potential at rotation angle: " << x << " " << potential.evaluate( t, v, true /*compute_deriv*/, deriv );
		TR << "  rot_deriv " << deriv.second << std::endl;
	}

	TR << "Enforcing continuity at pi" << std::endl;
	potential.set_enforce_continuity_at_pi( true );
	for ( Real const & x : xvals ) {
		v = x * v_norm;
		TR << "evaluating 6D potential at rotation angle: " << x << " " << potential.evaluate( t, v, true /*compute_deriv*/, deriv );
		TR << "  rot_deriv " << deriv.second << std::endl;
	}

}

/////////////////////////////////
void
read_tensor_test()
{
	using utility::tools::make_vector1;

	std::string const filename = option[ tensor_file ]();
	if ( !file::file_exists( filename ) ) utility_exit_with_message( "Supply a valid file with -tensor_file <*.bin.gz> or <*.txt.gz>" );

	TR <<  basic::database::full_name( "scoring/loop_close/6D_potentials/rna/loop_01/potential.txt.gz" ) << std::endl;

	// loop_01/potential.bin.gz --> value at this bin should be 2.4974.
	utility::vector1< numeric::Size > checkbins( make_vector1( 8, 11, 11, 6, 6, 6 )  );

	{
		numeric::MathNTensor< double, 6 > T;
		read_tensor_from_file( filename,  T );
		check_tensor( T, checkbins );

		utility::json_spirit::mObject json;
		read_tensor_from_file( filename, T, json );
		TR <<  write(json) << std::endl;
		check_tensor( T, checkbins );

		utility::fixedsizearray1< numeric::Real, 6 > minval, binwidth;
		core::scoring::loop_graph::evaluator::get_minval_binwidth( T, json, minval, binwidth );

		utility::fixedsizearray1< numeric::Real, 6 > xs;
		for ( numeric::Size n = 1; n <= 6; n++ ) xs[ n ] = binwidth[ n ] * ( checkbins[ n ] - 1 ) + minval[ n ];

		// check_spline( T, minval, binwidth, xs );

		check_multilinear_interpolation( T, minval, binwidth, xs );

		check_6D_potential( T, json, xs );
	}

}

////////////////////////////////
// testing accessor functions.
////////////////////////////////
void
check_tensor( numeric::MathNTensor< double, 6 > const & T,
	utility::vector1< numeric::Size > const & checkbins )
{
	using numeric::Size;
	using utility::tools::make_vector1;
	for ( numeric::Size n = 1; n <= 6; n++ ) TR << ' ' << T.n_bins( n );
	TR << std::endl;

	TR << T( make_vector1( checkbins[1]-1, checkbins[2]-1, checkbins[3]-1,
		checkbins[4]-1, checkbins[5]-1, checkbins[6]-1 ) ) << std::endl;
	TR << T(checkbins[1]-1, checkbins[2]-1, checkbins[3]-1,
		checkbins[4]-1, checkbins[5]-1, checkbins[6]-1  ) << std::endl;

}

void
check_multilinear_interpolation( numeric::MathNTensor< double, 6 > const & T,
	utility::fixedsizearray1< numeric::Real, 6 > const & minval,
	utility::fixedsizearray1< numeric::Real, 6 > const & binwidth,
	utility::fixedsizearray1< numeric::Real, 6 > const & xs )
{
	utility::fixedsizearray1< numeric::Real, 6 > xs_perturb, analytical_deriv, numerical_deriv;

	TR << "Checking ";
	for ( numeric::Size n = 1; n <= 6; n++ ) TR << ' ' << xs[ n ];
	TR << std::endl;

	numeric::Real val = multilinear_interpolation( T, minval, binwidth, xs, analytical_deriv );
	TR << "Value: " << val << std::endl;

	TR << "Analytical deriv:";
	for ( numeric::Size n = 1; n <= 6; n++ ) TR << ' ' << analytical_deriv[ n ];
	TR << std::endl;

	// perturb by a bit.
	Real delta( 1.0e-2 );
	for ( numeric::Size n = 1; n <= 6; n++ ) {
		xs_perturb = xs;
		xs_perturb[n] += delta;
		numerical_deriv[ n ] = ( multilinear_interpolation( T, minval, binwidth, xs_perturb ) - val ) / delta;
	}

	TR << "Numerical  deriv:";
	for ( numeric::Size n = 1; n <= 6; n++ ) TR << ' ' << numerical_deriv[ n ];
	TR << std::endl;

}

// try polycubic spline too
// This actually took 26 mins to run, so let's not make it default.
// And even saving the derivs to disk would cost 2^N * size( Tensor ), with N = 6, which would
// come out to ~1 Gb even after gzipping, I think.
void
check_spline( numeric::MathNTensor< double, 6 > const & T,
	utility::fixedsizearray1< numeric::Real, 6 > const & minval,
	utility::fixedsizearray1< numeric::Real, 6 > const & binwidth,
	utility::fixedsizearray1< numeric::Real, 6 > const & xs )
{
	using namespace utility::json_spirit;
	using namespace utility::tools;

	// adapted from core/scoring/MainchainScoreTable
	utility::fixedsizearray1< BorderFlag, 6 > borderflags; //Periodic boundaries.
	utility::fixedsizearray1< bool, 6 > lincont; //Meaningless argument for a polycubic spline with periodic boundary conditions.
	utility::fixedsizearray1< std::pair< numeric::Real, numeric::Real >, 6 > unused( std::pair<numeric::Real,numeric::Real>(0.0, 0.0) ); //Another meaningless argument for a polycubic spline with periodic boundary conditions.
	for ( numeric::Size i=1; i<=6; ++i ) {
		borderflags[i]   = e_Periodic;
		lincont[i]       = false;
	}
	PolycubicSplineOP<6> spline( new PolycubicSpline<6> );
	TR << "Training spline... " << std::endl;
	spline->train( borderflags, minval, binwidth, T, lincont, unused );

	TR << "Evaluate spline at ";
	for ( numeric::Size n = 1; n <= 6; n++ ) TR << " " << xs[ n ];
	TR << ": " << spline->F( xs ) << std::endl;
}

///////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	read_tensor_test();
	protocols::viewer::clear_conformation_viewers();
	exit( 0 );
}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
		NEW_OPT( tensor_file, "tensor_file", "default.bin.gz" );

		////////////////////////////////////////////////////////////////////////////
		// setup
		////////////////////////////////////////////////////////////////////////////
		devel::init(argc, argv);

		protocols::viewer::viewer_main( my_main );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}
