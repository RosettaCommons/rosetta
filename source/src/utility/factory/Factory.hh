// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/factory/Factory.hh
/// @brief  Pluggable factory
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
///
/// @remarks
///  @li Creates objects of registered concrete products in a common hierarchy
///  @li Product base class must typedef some types (see Types section)
///  @li Use Key pointers when keys are globals that may not be constructed yet
///  @li Use Key functions when keys are globals that may not be constructed yet or for generic registrants
///  @li Supports creation functions with 0 and 1 argument(s)
///  @li Place static Registrant members in concrete product classes to register with Factory


#ifndef INCLUDED_utility_factory_Factory_hh
#define INCLUDED_utility_factory_Factory_hh


// Unit headers
#include <utility/factory/Factory.fwd.hh>

// C++ headers
#include <utility/assert.hh>
#include <map>
#include <utility>


namespace utility {
namespace factory {


/// @brief Pluggable factory
template< typename P >
class Factory
{


public: // Types


	typedef  P  Product; // Product
	typedef  typename Product::FactoryKey  Key; // Product lookup key
	typedef  typename Product::FactoryKeyP  KeyP; // Product lookup key pointer
	typedef  typename Product::FactoryKeyFxn  KeyFxn; // Product key function
	typedef  typename Product::FactoryProductP  ProductP; // Product creation return pointer
	typedef  typename Product::FactoryArg  Arg; // Product creation argument
	typedef  typename Product::FactoryCreate  Create; // Product creation function


private: // Types


	typedef  std::map< Key, Create >  Registry; // Map from keys to creation functions
	typedef  std::map< KeyP, Create >  PtrRegistry; // Map from key pointers to creation functions
	typedef  std::map< KeyFxn, Create >  FxnRegistry; // Map from key functions to creation functions


private: // Creation


	/// @brief Default constructor
	inline
	Factory()
	= default;


public: // Static methods


	/// @brief Add a Product to the registry
	inline
	static
	void
	add( Key const & key, Create create )
	{
		registry().insert( std::make_pair( key, create ) );
	}


	/// @brief Add a Product to the pointer registry
	inline
	static
	void
	add( KeyP const & key_p, Create create )
	{
		ptr_registry().insert( std::make_pair( key_p, create ) );
	}


	/// @brief Add a Product to the function registry
	inline
	static
	void
	add( KeyFxn key_fxn, Create create )
	{
		fxn_registry().insert( std::make_pair( key_fxn, create ) );
	}


	/// @brief Remove a Product from the registry
	inline
	static
	void
	remove( Key const & key )
	{
		registry().erase( key );
	}


	/// @brief Remove a Product from the pointer registry
	inline
	static
	void
	remove( KeyP const & key_p )
	{
		ptr_registry().erase( key_p );
	}


	/// @brief Remove a Product from the function registry
	inline
	static
	void
	remove( KeyFxn key_fxn )
	{
		fxn_registry().erase( key_fxn );
	}


	/// @brief Has a Product with key?
	inline
	static
	bool
	has( Key const & key )
	{
		Registry & reg( registry() );
		return ( reg.find( key ) != reg.end() );
	}


	/// @brief Create a Product from key
	inline
	static
	ProductP
	create( Key const & key )
	{
		Registry & reg( registry() );
	debug_assert( reg.find( key ) != reg.end() );
		return reg.find( key )->second();
	}


	/// @brief Create a Product from key and 1 argument
	inline
	static
	ProductP
	create( Key const & key, Arg & arg )
	{
		Registry & reg( registry() );
	debug_assert( reg.find( key ) != reg.end() );
		return reg.find( key )->second( arg );
	}


private: // Static methods


	/// @brief Key registry
	/// @note  Precondition: Keys referenced by pointer and function registries must have been constructed
	inline
	static
	Registry &
	registry()
	{
		static Registry registry_; // Function local to avoid initialization order issues
		ptr_registry_transfer( registry_ );
		fxn_registry_transfer( registry_ );
		return registry_;
	}


	/// @brief Key pointer registry
	inline
	static
	PtrRegistry &
	ptr_registry()
	{
		static PtrRegistry ptr_registry_; // Function local to avoid initialization order issues
		return ptr_registry_;
	}


	/// @brief Key function registry
	inline
	static
	FxnRegistry &
	fxn_registry()
	{
		static FxnRegistry fxn_registry_; // Function local to avoid initialization order issues
		return fxn_registry_;
	}


	/// @brief Transfer key pointer registry to key registry
	inline
	static
	void
	ptr_registry_transfer( Registry & reg )
	{
		PtrRegistry & ptr_reg( ptr_registry() );
		if ( ! ptr_reg.empty() ) {
			for ( typename PtrRegistry::const_iterator i = ptr_reg.begin(), e = ptr_reg.end(); i != e; ++i ) {
				reg[ *(i->first) ] = i->second;
			}
			ptr_reg.clear();
		}
	}


	/// @brief Transfer key function registry to key registry
	inline
	static
	void
	fxn_registry_transfer( Registry & reg )
	{
		FxnRegistry & fxn_reg( fxn_registry() );
		if ( ! fxn_reg.empty() ) {
			for ( typename FxnRegistry::const_iterator i = fxn_reg.begin(), e = fxn_reg.end(); i != e; ++i ) {
				reg[ i->first() ] = i->second;
			}
			fxn_reg.clear();
		}
	}


}; // Factory


} // namespace factory
} // namespace utility


#endif // INCLUDED_utility_factory_Factory_HH
