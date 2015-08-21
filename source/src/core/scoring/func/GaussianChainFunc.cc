// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/func/GaussianChainFunc.hh
/// @brief Definition for functions used in loop closure terms.
/// @author Rhiju Das


#include <core/scoring/func/GaussianChainFunc.hh>
#include <core/scoring/func/GaussianChainSingleFunc.hh>
#include <core/scoring/func/GaussianChainDoubleFunc.hh>
#include <core/scoring/func/GaussianChainTripleFunc.hh>
#include <core/scoring/func/GaussianChainQuadrupleFunc.hh>
#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <utility/exit.hh>

// C++ Headers
#include <iostream>
#include <cmath>

//////////////////////////////////////////////////////////////////////////
// ---------------------------------
// Closure energies for loop cycles.
// ---------------------------------
//
//          D1
//       ------~
//      ~       | D4
//     ~        |
//    ~         ~
//    |        ~
// D2 |       /
//    |      /  D3
//     ~~~~~/
//
//  D1, D2, D3, D4 mark rigid joints, with lengths D1, ... D4.
//  Squiggles (~~~) mark gaussian chains,
//   with total variance of gaussian_variance.
//
//  In simplest case (one joint D, gaussian variance sigma^2 ):
//
//  P( closure ) = (capture volume) x ( 1 / 2 pi sigma^2 )^(3/2)
//                                  x exp( - D / 2 sigma^2 )
//
//  [For folks used to thinking about mean end-to-end distance L,
//        sigma^2 = L^2 / 3 ].
//
//  Following includes explicit formulae for one, two, three, and four joints.
//  There is also a general formula which would be straightforward to implement,
//  which I have worked out (and is basically implemented in the four-joint
//  function GaussianChainQuadrupleFunc.cc)
//
// See core/scoring/loop_graph/LoopGraph.cc for use case.
//
// More notes at:
//   https://docs.google.com/a/stanford.edu/file/d/0B6gpwdY_Bgd0OHdzVWJGTHBvTzg/edit
//
// Rhiju, Sep. 2013
//
//////////////////////////////////////////////////////////////////////////


namespace core {
namespace scoring {
namespace func {

GaussianChainFunc::GaussianChainFunc( Real const gaussian_variance,
	Real const loop_fixed_cost,
	utility::vector1< Real > const & other_distances ):
	gaussian_variance_( gaussian_variance ),
	loop_fixed_cost_( loop_fixed_cost ),
	other_distances_( other_distances )
{
	initialize_parameters();
	initialize_func();
}

GaussianChainFunc::GaussianChainFunc( Real const gaussian_variance,
	Real const loop_fixed_cost ):
	gaussian_variance_( gaussian_variance ),
	loop_fixed_cost_( loop_fixed_cost )
{
	initialize_parameters();
	initialize_func();
}

/////////////////////////////////////////////////////
void
GaussianChainFunc::initialize_parameters(){
	kB_T_ = 1.0; // choice of energy units.
	force_combined_gaussian_approximation_ = false;
}

/////////////////////////////////////////////////////
FuncOP
GaussianChainFunc::clone() const
{
	return FuncOP( new GaussianChainFunc( gaussian_variance_, loop_fixed_cost_, other_distances_ ) );
}

/////////////////////////////////////////////////////
void
GaussianChainFunc::initialize_func(){

	func_.reset();

	if ( !force_combined_gaussian_approximation_ ) {
		if ( other_distances_.size() == 0 ) {
			func_ = FuncOP( new GaussianChainSingleFunc( gaussian_variance_, loop_fixed_cost_ ) );
		} else if ( other_distances_.size() == 1 ) {
			func_ = FuncOP( new GaussianChainDoubleFunc( gaussian_variance_, loop_fixed_cost_, other_distances_[1] ) );
		} else if ( other_distances_.size() == 2 ) {
			func_ = FuncOP( new GaussianChainTripleFunc( gaussian_variance_, loop_fixed_cost_, other_distances_[1], other_distances_[2] ) );
		} else if ( other_distances_.size() == 3 ) {
			func_ = FuncOP( new GaussianChainQuadrupleFunc( gaussian_variance_, loop_fixed_cost_, other_distances_[1], other_distances_[2], other_distances_[3] ) );
		}
	}

	if ( func_ == 0 ) {
		// for more than four joints, treat as one effective long gaussian chain.
		// not a bad approximation actually.
		Real gaussian_variance_total  = gaussian_variance_;
		for ( Size n = 1; n <= other_distances_.size(); n++ ) {
			// Note that the radius of gyration is sqrt( 3 ) * gaussian_variance.
			gaussian_variance_total += (other_distances_[ n ] * other_distances_[ n ] / 3.0);
		}
		func_ = FuncOP( new GaussianChainSingleFunc( gaussian_variance_total,  loop_fixed_cost_ ) );
	}

	runtime_assert( func_.get() != NULL );
}

/////////////////////////////////////////////////////
Real
GaussianChainFunc::func( Real const z ) const
{
	return ( func_->func( z ) );
}

/////////////////////////////////////////////////////
Real
GaussianChainFunc::dfunc( Real const z ) const
{
	return ( func_->dfunc( z ) );
}

/////////////////////////////////////////////////////
void
GaussianChainFunc::read_data( std::istream & in ) {
	in >> loop_fixed_cost_ >> gaussian_variance_;
	Real dist;
	if ( in.good() ) {
		in >> dist;
		other_distances_.push_back( dist );
	}
	initialize_func();
}

void
GaussianChainFunc::show_definition( std::ostream &out ) const {
	func_->show_definition( out );
}


} // namespace constraints
} // namespace scoring
} // namespace core
