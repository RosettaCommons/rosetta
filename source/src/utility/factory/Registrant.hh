// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/factory/Registrant.hh
/// @brief  Factory registrant
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
///
/// @remarks
///  @li Registers a concrete product class with its factory when constructed
///  @li Use Key pointers when keys are globals that may not be constructed yet


#ifndef INCLUDED_utility_factory_Registrant_hh
#define INCLUDED_utility_factory_Registrant_hh


// Unit headers
#include <utility/factory/Registrant.fwd.hh>

// Package headers
#include <utility/factory/Factory.hh>


namespace utility {
namespace factory {


/// @brief Factory registrant
template< typename P >
class Registrant
{


public: // Types


	typedef  P  Product; // Product
	typedef  Factory< Product >  ProductFactory; // Factory
	typedef  typename ProductFactory::Key  Key; // Product lookup key
	typedef  typename ProductFactory::KeyP  KeyP; // Product lookup key pointer
	typedef  typename ProductFactory::KeyFxn  KeyFxn; // Product lookup key function
	typedef  typename ProductFactory::Create  Create; // Product creation function


private: // Creation


	/// @brief Copy constructor
	Registrant( Registrant const & ); // Undefined


public: // Creation


	/// @brief 1 Key constructor
	inline
	Registrant(
		Key const & key,
		Create create
	)
	{
		ProductFactory::add( key, create );
	}


	/// @brief 1 Key pointer constructor
	inline
	Registrant(
		KeyP const & key_p,
		Create create
	)
	{
		ProductFactory::add( key_p, create );
	}


	/// @brief 1 Key function constructor
	inline
	Registrant(
		KeyFxn key_fxn,
		Create create
	)
	{
		ProductFactory::add( key_fxn, create );
	}


	/// @brief 2 Key constructor
	inline
	Registrant(
		Key const & key1,
		Key const & key2,
		Create create
	)
	{
		ProductFactory::add( key1, create );
		ProductFactory::add( key2, create );
	}


	/// @brief 2 Key pointer constructor
	inline
	Registrant(
		KeyP const & key1_p,
		KeyP const & key2_p,
		Create create
	)
	{
		ProductFactory::add( key1_p, create );
		ProductFactory::add( key2_p, create );
	}


	/// @brief 2 Key function constructor
	inline
	Registrant(
		KeyFxn key1_fxn,
		KeyFxn key2_fxn,
		Create create
	)
	{
		ProductFactory::add( key1_fxn, create );
		ProductFactory::add( key2_fxn, create );
	}


	/// @brief 3 Key constructor
	inline
	Registrant(
		Key const & key1,
		Key const & key2,
		Key const & key3,
		Create create
	)
	{
		ProductFactory::add( key1, create );
		ProductFactory::add( key2, create );
		ProductFactory::add( key3, create );
	}


	/// @brief 3 Key pointer constructor
	inline
	Registrant(
		KeyP const & key1_p,
		KeyP const & key2_p,
		KeyP const & key3_p,
		Create create
	)
	{
		ProductFactory::add( key1_p, create );
		ProductFactory::add( key2_p, create );
		ProductFactory::add( key3_p, create );
	}


	/// @brief 3 Key function constructor
	inline
	Registrant(
		KeyFxn key1_fxn,
		KeyFxn key2_fxn,
		KeyFxn key3_fxn,
		Create create
	)
	{
		ProductFactory::add( key1_fxn, create );
		ProductFactory::add( key2_fxn, create );
		ProductFactory::add( key3_fxn, create );
	}


private: // Assignment


	/// @brief Copy assignment
	Registrant &
	operator =( Registrant const & ); // Undefined


}; // Registrant


} // namespace factory
} // namespace utility


#endif // INCLUDED_utility_factory_Registrant_HH
