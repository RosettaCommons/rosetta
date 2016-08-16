// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/keys/KeyLess.hh
/// @brief  Key comparison functor template
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
///
/// @note
///  @li Types must have key() member functions and comparable keys


#ifndef INCLUDED_utility_keys_KeyLess_hh
#define INCLUDED_utility_keys_KeyLess_hh


namespace utility {
namespace keys {


/// @brief Key member comparison functor template
template< typename T, typename U >
class KeyLess
{


public: // Methods


	/// @brief Functor operator
	inline
	bool
	operator ()( T const & t, U const & u ) const
	{
		return ( t.key() < u.key() );
	}


}; // KeyLess


/// @brief Key member comparison functor template for pointers
template< typename T, typename U >
class PointerKeyLess
{


public: // Methods


	/// @brief Functor operator
	inline
	bool
	operator ()( T const & t, U const & u ) const
	{
		return ( t->key() < u->key() );
	}


}; // PointerKeyLess


} // namespace keys
} // namespace utility


#endif // INCLUDED_utility_keys_KeyLess_HH
