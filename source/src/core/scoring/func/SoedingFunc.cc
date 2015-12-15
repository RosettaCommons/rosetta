// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/func/SoedingFunc.hh
/// @brief Definition for functions used in definition of constraints.
/// @author James Thompson


#include <core/scoring/func/SoedingFunc.hh>
#include <core/scoring/constraints/util.hh>
#include <core/types.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/james.OptionKeys.gen.hh>

// C++ Headers
#include <iostream>

#include <utility/vector1.hh>

//Auto Headers
#include <cmath>


#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace func {

using namespace core::scoring::constraints;

bool SoedingFunc::operator == ( Func const & other ) const
{
	if ( ! same_type_as_me( other ) ) return false;
	if ( ! other.same_type_as_me( *this ) ) return false;

	SoedingFunc const & other_downcast( static_cast< SoedingFunc const & > (other) );
	if ( w1_    != other_downcast.w1_    ) return false;
	if ( mean1_ != other_downcast.mean1_ ) return false;
	if ( sdev1_ != other_downcast.sdev1_ ) return false;
	if ( w2_    != other_downcast.w2_    ) return false;
	if ( mean2_ != other_downcast.mean2_ ) return false;
	if ( sdev2_ != other_downcast.sdev2_ ) return false;
	return true;
}

bool SoedingFunc::same_type_as_me( Func const & other ) const
{
	return dynamic_cast< SoedingFunc const * > ( &other );
}

void
SoedingFunc::read_data( std::istream & in ) {
	in >> w1_ >> mean1_ >> sdev1_ >> w2_ >> mean2_ >> sdev2_;
}

Real
SoedingFunc::compute_func( Real const x ) const {
	using std::exp;
	using std::log;

	Real const numerator(
		dgaussian(x,mean1_,sdev1_,w1_) + dgaussian(x,mean2_,sdev2_,w2_)
	);
	Real const denominator(
		dgaussian(x, mean2_, sdev2_,w2_)
	);

	//Real const score( log(numerator) - log(denominator) );
	Real const score( log(denominator) - log(numerator) );
	return score;
}

Real
SoedingFunc::func( Real const x ) const {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	Real const score( compute_func(x) );

	if ( option[ james::debug ]() ) {
		//return std::min( score, compute_func(10) );
		return std::min( score, 0.0 );
	}
	return score;
} // func

Real
SoedingFunc::dfunc( Real const x ) const {
	return estimate_dfunc(x);
} // dfunc

void SoedingFunc::show_definition( std::ostream & out ) const {
	out << "SOEDINGFUNC " << w1_ << " " << mean1_ << " " << sdev1_ << " "
		<< w2_ << " " << mean2_ << " " << sdev2_;
} // show_definition

} // namespace constraints
} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::func::SoedingFunc::save( Archive & arc ) const {
	arc( cereal::base_class< Func >( this ) );
	arc( CEREAL_NVP( w1_ ) ); // Real
	arc( CEREAL_NVP( mean1_ ) ); // Real
	arc( CEREAL_NVP( sdev1_ ) ); // Real
	arc( CEREAL_NVP( w2_ ) ); // Real
	arc( CEREAL_NVP( mean2_ ) ); // Real
	arc( CEREAL_NVP( sdev2_ ) ); // Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::func::SoedingFunc::load( Archive & arc ) {
	arc( cereal::base_class< Func >( this ) );
	arc( w1_ ); // Real
	arc( mean1_ ); // Real
	arc( sdev1_ ); // Real
	arc( w2_ ); // Real
	arc( mean2_ ); // Real
	arc( sdev2_ ); // Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::func::SoedingFunc );
CEREAL_REGISTER_TYPE( core::scoring::func::SoedingFunc )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_func_SoedingFunc )
#endif // SERIALIZATION
