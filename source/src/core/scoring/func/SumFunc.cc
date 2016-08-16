// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/func/SumFunc.cc
/// @brief Weighted constraint function that reweights other constraints
/// by a constant scalar value.
/// @author James Thompson, Greg Taylor


#include <core/scoring/func/SumFunc.hh>
#include <core/scoring/func/FuncFactory.hh>

#include <core/types.hh>

#include <utility/pointer/ReferenceCount.hh>


// C++ Headers

#include <sstream>


#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#include <utility/vector1.srlz.hh>

// Cereal headers
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace func {

bool SumFunc::operator == ( Func const & other ) const
{
	if ( ! same_type_as_me( other ) ) return false;
	if ( ! other.same_type_as_me( *this ) ) return false;

	SumFunc const & other_downcast( static_cast< SumFunc const & > (other) );
	if ( funcs_.size() != other_downcast.funcs_.size() ) return false;
	for ( Size ii = 1; ii <= funcs_.size(); ++ii ) {
		if ( funcs_[ ii ] == other_downcast.funcs_[ ii ] ) continue;
		if ( funcs_[ii] && other_downcast.funcs_[ ii ] && ! (*funcs_[ii] == *other_downcast.funcs_[ii] ) ) return false;
	}
	return true;
}

bool SumFunc::same_type_as_me( Func const & other ) const
{
	return dynamic_cast< SumFunc const * > ( &other );
}


Real
SumFunc::func( Real const x ) const {
	Real f( 0.0 );
	for ( utility::vector1< FuncOP >::const_iterator it = funcs_.begin(), end = funcs_.end();
			it != end; ++it
			) {
		f += (*it)->func( x );
	}

	return f;
}

Real
SumFunc::dfunc( Real const x ) const {
	Real df( 0.0 );
	for ( utility::vector1< FuncOP >::const_iterator it = funcs_.begin(), end = funcs_.end();
			it != end; ++it
			) {
		df += (*it)->dfunc( x );
	}

	return df;
}

void
SumFunc::read_data( std::istream& in ) {
	FuncFactory func_factory;
	std::string func_type;

	Size n_funcs( 0 );
	in >> n_funcs;
	if ( in.fail() ) return;

	debug_assert( n_funcs >= 1 );

	for ( Size i = 1; i <= n_funcs; ++i ) {
		in >> func_type;
		FuncOP current_func;

		current_func = func_factory.new_func( func_type );
		current_func->read_data( in );
		add_func( current_func );
	}
} // read_data

} // namespace constraints
} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::func::SumFunc::save( Archive & arc ) const {
	arc( cereal::base_class< Func >( this ) );
	arc( CEREAL_NVP( funcs_ ) ); // utility::vector1<FuncOP>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::func::SumFunc::load( Archive & arc ) {
	arc( cereal::base_class< Func >( this ) );
	arc( funcs_ ); // utility::vector1<FuncOP>
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::func::SumFunc );
CEREAL_REGISTER_TYPE( core::scoring::func::SumFunc )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_func_SumFunc )
#endif // SERIALIZATION
