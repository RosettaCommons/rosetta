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


// libRosetta headers
#include <core/scoring/func/GaussianChainFunc.hh>
#include <utility/tools/make_vector1.hh>
#include <iostream>
#include <ObjexxFCL/format.hh>

using namespace core;
using namespace core::scoring::func;
using ObjexxFCL::format::F;
using utility::tools::make_vector1;

///////////////////////////////////////////////////////////////
// This could become a unit test.
///////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////
void
double_func_test( Distance const distance2  ){

	utility::vector1< Distance > other_distances;
	other_distances.push_back( distance2 );

	Real const gaussian_variance( 5.0 * 5.0 );
	FuncOP func1 = new GaussianChainFunc( gaussian_variance ); // single func
	FuncOP func2 = new GaussianChainFunc( gaussian_variance, other_distances );
	GaussianChainFuncOP func2_approx = new GaussianChainFunc( gaussian_variance, other_distances );
	func2_approx->set_force_combined_gaussian_approximation( true );

	// sweep through distance1.
	std::cout << "-- GAUSSIAN_CHAIN_DOUBLE_FUNC --" << std::endl;
	std::cout << "Fixing Distance2 at " << distance2 << std::endl;
	for ( Size i = 1; i <= 20; i++ ){
		Distance const d = Real( i ) * 1.0;
		std::cout << "Distance1 " << F(6,2,d) << "   FUNC_EXACT " << F(10,5,func2->func( d )) << " FUNC_APPROX " << F(10,5,func2_approx->func( d )) << " SINGLE_FUNC " << F(10,5,func1->func( d )) << std::endl;
	}
	std::cout << std::endl;
}

void
triple_func_test( Distance const distance2, Distance const distance3  ){
	utility::vector1< Distance > other_distances;
	other_distances.push_back( distance2 );
	other_distances.push_back( distance3 );

	Real const gaussian_variance( 5.0 * 5.0 );
	FuncOP func1 = new GaussianChainFunc( gaussian_variance ); // single func
	FuncOP func2 = new GaussianChainFunc( gaussian_variance, make_vector1( distance2 ) );
	FuncOP func3 = new GaussianChainFunc( gaussian_variance, other_distances );
	GaussianChainFuncOP func3_approx = new GaussianChainFunc( gaussian_variance, other_distances );
	func3_approx->set_force_combined_gaussian_approximation( true );

	FuncOP func0 = new GaussianChainFunc( gaussian_variance ); // single func

	// sweep through distance1.
	std::cout << "-- GAUSSIAN_CHAIN_TRIPLE_FUNC --" << std::endl;
	std::cout << "Fixing Distance2 at " << distance2 << ", Distance3 at " << distance3 << std::endl;
	for ( Size i = 1; i <= 20; i++ ){
		Distance const d = Real( i ) * 1.0;
		std::cout << "Distance1 " << F(6,2,d) << "   FUNC_EXACT " << F(10,5,func3->func( d )) << " FUNC_APPROX " << F(10,5,func3_approx->func( d )) << " DOUBLE_FUNC " << F(10,5,func2->func( d )) << " SINGLE_FUNC " << F(10,5,func1->func( d ) ) << std::endl;
	}
	std::cout << std::endl;
}


void
quadruple_func_test( Distance const distance2, Distance const distance3, Distance const distance4  ){
	utility::vector1< Distance > other_distances;
	other_distances.push_back( distance2 );
	other_distances.push_back( distance3 );
	other_distances.push_back( distance4 );

	Real const gaussian_variance( 5.0 * 5.0 );
	//	FuncOP func1 = new GaussianChainFunc( gaussian_variance ); // single func
	//	FuncOP func2 = new GaussianChainFunc( gaussian_variance, make_vector1( distance2 ) );
	FuncOP func4 = new GaussianChainFunc( gaussian_variance, other_distances );
	GaussianChainFuncOP func4_approx = new GaussianChainFunc( gaussian_variance, other_distances );
	func4_approx->set_force_combined_gaussian_approximation( true );

	FuncOP func0 = new GaussianChainFunc( gaussian_variance ); // single func

	// sweep through distance1.
	std::cout << "-- GAUSSIAN_CHAIN_QUADRUPLE_FUNC --" << std::endl;
	std::cout << "Fixing Distance2 at " << distance2 << ", Distance3 at " << distance3  << ", Distance4 at " << distance4 << std::endl;
	for ( Size i = 1; i <= 20; i++ ){
		Distance const d = Real( i ) * 1.0;
		std::cout << "Distance1 " << F(6,2,d) << "   FUNC_EXACT " << F(10,5,func4->func( d )) << " FUNC_APPROX " << F(10,5,func4_approx->func( d )) << std::endl; //<< " DOUBLE_FUNC " << F(10,5,func2->func( d )) << " SINGLE_FUNC " << F(10,5,func1->func( d ) ) << std::endl;
	}
	std::cout << std::endl;
}

void
func_test(){

	double_func_test( 0.5 );
	double_func_test( 15.0 );

	triple_func_test( 0.5, 0.5 );
	triple_func_test( 15.0, 0.5 );
	triple_func_test( 15.0, 15.0 );

	quadruple_func_test( 0.01, 0.01, 0.01 );
	quadruple_func_test( 15.0, 15.0, 15.0 );

}


///////////////////////////////////////////////////////////////////////////////
int
main()
{
	func_test();
	exit( 0 );
}
