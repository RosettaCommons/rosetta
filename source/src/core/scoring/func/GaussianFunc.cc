// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/func/GaussianFunc.hh
/// @brief Definition for functions used in definition of constraints.
/// @author James Thompson


#include <core/scoring/func/GaussianFunc.hh>
#include <core/scoring/constraints/util.hh>
#include <core/types.hh>

#include <utility/pointer/ReferenceCount.hh>

// C++ Headers

#include <iostream>

#include <utility/vector1.hh>


#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace func {

using namespace core::scoring::constraints;

bool GaussianFunc::operator == ( Func const & other ) const
{
	if ( ! same_type_as_me( other ) ) return false;
	if ( ! other.same_type_as_me( *this ) ) return false;

	auto const & other_downcast( static_cast< GaussianFunc const & > (other) );
	if ( mean_          != other_downcast.mean_          ) return false;
	if ( sd_            != other_downcast.sd_            ) return false;
	if ( use_log_score_ != other_downcast.use_log_score_ ) return false;
	if ( weight_        != other_downcast.weight_        ) return false;

	return true;
}

bool GaussianFunc::same_type_as_me( Func const & other ) const
{
	return dynamic_cast< GaussianFunc const * > ( &other );
}

void
GaussianFunc::read_data( std::istream& in ) {
	in  >> mean_ >> sd_;

	if ( !in.eof() ) {
		std::string tag;
		in >> tag;

		if ( tag == "NOLOG" ) {
			use_log_score_ = false;
		} else if ( tag == "WEIGHT" ) {
			in >> weight_;
			runtime_assert_string_msg( !in.fail(), "Error parsing Gaussian function definition.  A \"WEIGHT\" statement must be followed by a floating-point value!" );
		}

		if ( !in.eof() ) {
			in >> tag;
			if ( tag == "WEIGHT" ) {
				in >> weight_;
				runtime_assert_string_msg( !in.fail(), "Error parsing Gaussian function definition.  A \"WEIGHT\" statement must be followed by a floating-point value!" );
			}
		}
	}

}

Real
GaussianFunc::func( Real const x ) const {
	if ( use_log_score_ ) return - weight_*logdgaussian( x, mean_, sd_, 1 );
	else                  return weight_*dgaussian( x, mean_, sd_, 1 );
} // func

Real
GaussianFunc::dfunc( Real const x ) const {
	if ( use_log_score_ ) return - weight_*logdgaussian_deriv( x, mean_, sd_, 1 );
	else                  return weight_*gaussian_deriv( x, mean_, sd_, 1 );
} // dfunc

void GaussianFunc::show_definition( std::ostream& out ) const {
	out << "GAUSSIANFUNC " << mean_ << ' ' << sd_ << ' ' << weight_ << "\n";
}

} // namespace constraints
} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::scoring::func::GaussianFunc::GaussianFunc() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::func::GaussianFunc::save( Archive & arc ) const {
	arc( cereal::base_class< Func >( this ) );
	arc( CEREAL_NVP( mean_ ) ); // Real
	arc( CEREAL_NVP( sd_ ) ); // Real
	arc( CEREAL_NVP( use_log_score_ ) ); // _Bool
	arc( CEREAL_NVP( weight_ ) ); // Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::func::GaussianFunc::load( Archive & arc ) {
	arc( cereal::base_class< Func >( this ) );
	arc( mean_ ); // Real
	arc( sd_ ); // Real
	arc( use_log_score_ ); // _Bool
	arc( weight_ ); // Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::func::GaussianFunc );
CEREAL_REGISTER_TYPE( core::scoring::func::GaussianFunc )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_func_GaussianFunc )
#endif // SERIALIZATION
