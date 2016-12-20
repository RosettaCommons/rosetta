// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// pilot app for testing Gaussian Chain Func - rhiju, Dec. 2016
// Move to unit test when done.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// libRosetta headers
#include <core/scoring/func/GaussianChainFunc.hh>
#include <core/scoring/func/GaussianChainGeneralFunc.hh>
#include <utility/tools/make_vector1.hh>
#include <core/types.hh>
#include <protocols/viewer/viewers.hh>
#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.rhiju.generalize_gaussian_chain_func_test" );


using namespace utility;
using namespace utility::tools;
using namespace core;

Real
get_gaussian_chain_func( vector1< Distance > const & dists,
												 Real const & gaussian_variance,
												 Real const & loop_fixed_cost,
												 Size const idx,
												 bool const force_use_general )
{
	using namespace core::scoring::func;
	Distance dist( dists[ idx ] );
	vector1< Distance > other_dists;
	for ( Size k = 1; k <= dists.size(); k++ ) {
		if ( k != idx ) other_dists.push_back( dists[ k ] );
	}
	FuncOP func;
	if ( force_use_general ) {
		func = FuncOP( new GaussianChainGeneralFunc( gaussian_variance, loop_fixed_cost, other_dists ) );
	} else {
		func = FuncOP( new GaussianChainFunc( gaussian_variance, loop_fixed_cost, other_dists ) );
	}

	// test derivs.
	Distance delta = 1.0e-6;
	TR << "Checking deriv " << ( func->func(dist + delta) - func->func(dist ) )/delta << " vs " << func->dfunc( dist ) << std::endl;

	return ( func->func( dist ) );
}

// option key includes
///////////////////////////////////////////////////////////////////////////////
void
generalize_gaussian_chain_func_test()
{
	vector1< Real > all_dists = make_vector1(  13.7201, 9.0501, 10.7559 , 16.4339 , 9.93362 );
	Real const gaussian_variance = 20.0 * 20.0;
	Real loop_fixed_cost = 0.0;
	for ( Size n = 1; n <= all_dists.size(); n++ ) {

		vector1< Real > dists;
		for ( Size k = 1; k <= n; k++ ) dists.push_back( all_dists[ k ] );

		Real val1 = get_gaussian_chain_func( dists, gaussian_variance, loop_fixed_cost, 1, false );
		Real val2 = get_gaussian_chain_func( dists, gaussian_variance, loop_fixed_cost, std::min( Size(2), n ), false );
		Real val3 = get_gaussian_chain_func( dists, gaussian_variance, loop_fixed_cost, 1, true );
		Real val4 = get_gaussian_chain_func( dists, gaussian_variance, loop_fixed_cost, std::min( Size(2), n ), true );
		TR << val1 << " vs " << val2 << " vs " << val3 << " vs " << val4 << " better be equal!" << std::endl;

	}

}


///////////////////////////////////////////////////////////////
void*
my_main( void* )
{

	using namespace basic::options;

	generalize_gaussian_chain_func_test();
	protocols::viewer::clear_conformation_viewers();
	exit( 0 );

}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	try {

		using namespace basic::options;

		////////////////////////////////////////////////////////////////////////////
		// setup
		////////////////////////////////////////////////////////////////////////////
		devel::init(argc, argv);


		////////////////////////////////////////////////////////////////////////////
		// end of setup
		////////////////////////////////////////////////////////////////////////////

		protocols::viewer::viewer_main( my_main );


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
