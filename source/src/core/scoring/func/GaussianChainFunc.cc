// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/func/GaussianChainFunc.cc
/// @brief Definition for functions used in loop closure terms.
/// @author Rhiju Das


#include <core/scoring/func/GaussianChainFunc.hh>
#include <core/scoring/func/GaussianChainSingleFunc.hh>
#include <core/scoring/func/GaussianChainDoubleFunc.hh>
#include <core/scoring/func/GaussianChainTripleFunc.hh>
#include <core/scoring/func/GaussianChainQuadrupleFunc.hh>
#include <core/scoring/func/GaussianChainGeneralFunc.hh>
#include <core/types.hh>
#include <utility>
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


#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#include <utility/vector1.srlz.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


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

bool GaussianChainFunc::operator == ( Func const & other ) const
{
	if ( ! same_type_as_me( other ) ) return false;
	if ( ! other.same_type_as_me( *this ) ) return false;

	auto const & other_downcast( static_cast< GaussianChainFunc const & > (other) );
	if ( gaussian_variance_     != other_downcast.gaussian_variance_     ) return false;
	if ( loop_fixed_cost_       != other_downcast.loop_fixed_cost_       ) return false;
	if ( other_distances_       != other_downcast.other_distances_       ) return false;
	if ( kB_T_                  != other_downcast.kB_T_                  ) return false;
	if ( force_combined_gaussian_approximation_ != other_downcast.force_combined_gaussian_approximation_ ) return false;

	return func_ == other_downcast.func_ || ( func_ && other_downcast.func_ && *func_ == *other_downcast.func_ );
}

bool GaussianChainFunc::same_type_as_me( Func const & other ) const
{
	return dynamic_cast< GaussianChainFunc const * > ( &other );
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
		} else {
			func_ = FuncOP( new GaussianChainGeneralFunc( gaussian_variance_, loop_fixed_cost_, other_distances_ ) );
		}
	}

	if ( func_ == nullptr ) {
		// at an early stage, was treating as one effective long gaussian chain,
		// but note that GaussianChainGeneralFunc (see above) replaces this.
		Real gaussian_variance_total  = gaussian_variance_;
		for ( Size n = 1; n <= other_distances_.size(); n++ ) {
			// Note that the radius of gyration is sqrt( 3 ) * gaussian_variance.
			gaussian_variance_total += (other_distances_[ n ] * other_distances_[ n ] / 3.0);
		}
		func_ = FuncOP( new GaussianChainSingleFunc( gaussian_variance_total,  loop_fixed_cost_ ) );
	}

	runtime_assert( func_.get() != nullptr );
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

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::scoring::func::GaussianChainFunc::GaussianChainFunc() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::func::GaussianChainFunc::save( Archive & arc ) const {
	arc( cereal::base_class< Func >( this ) );
	arc( CEREAL_NVP( gaussian_variance_ ) ); // Real
	arc( CEREAL_NVP( loop_fixed_cost_ ) ); // Real
	arc( CEREAL_NVP( other_distances_ ) ); // utility::vector1<Real>
	arc( CEREAL_NVP( kB_T_ ) ); // Real
	arc( CEREAL_NVP( force_combined_gaussian_approximation_ ) ); // _Bool
	arc( CEREAL_NVP( func_ ) ); // FuncOP
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::func::GaussianChainFunc::load( Archive & arc ) {
	arc( cereal::base_class< Func >( this ) );
	arc( gaussian_variance_ ); // Real
	arc( loop_fixed_cost_ ); // Real
	arc( other_distances_ ); // utility::vector1<Real>
	arc( kB_T_ ); // Real
	arc( force_combined_gaussian_approximation_ ); // _Bool
	arc( func_ ); // FuncOP
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::func::GaussianChainFunc );
CEREAL_REGISTER_TYPE( core::scoring::func::GaussianChainFunc )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_func_GaussianChainFunc )
#endif // SERIALIZATION
