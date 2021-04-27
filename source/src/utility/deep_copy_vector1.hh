// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/deep_copy_vector1.hh
/// @brief  A vector1 which does a deep copy of its contents
///
/// See utility/pointer/deep_copy.hh for more information about the deep copy framwork.
///
/// One limitation of the DeepCopyOP framework is that it doesn't work well in containers.
/// Which means if you have a vector1< MyClassOP > that you want to deep copy, well, things don't work so well.
/// This class attempts to fix that issue, by providing a vector1 class which is deep copied by default when
/// the containing class is default copied.
///
/// Like the DeepCopyOP, this class should only be used for member variables.
/// Don't pass it around. It should freely convert into a vector1 reference if you need to pass around a reference to it.
///
/// The other thing to note is that this class only really makes sense for OP contents.
/// If you're containing by value (instead of by OP), the default semantics of vector1 is what you want.
///
/// @author Rocco Moretti (rmorettiase@gmail.com)


#ifndef INCLUDED_utility_deep_copy_vector1_HH
#define INCLUDED_utility_deep_copy_vector1_HH


// Unit headers
#include <utility/vector1.hh>

// Package headers
#include <utility/vectorL.hh>

namespace utility {

template< typename T >
class deep_copy_vector1 : public vector1< T >
{

public:

	typedef  vector1< T >  super;

	// We can re-use most of the contructors
	// We only need to special-case the copy constructors.
	using vector1<T>::vector1;

	// @details It's a little concerning that the constructor isn't being inherited
	deep_copy_vector1() :
		vector1<T>()
	{}

	/// @brief Copy constructor
	deep_copy_vector1( deep_copy_vector1 const & v ) :
		vector1<T>()
	{
		this->reserve( v.size() );
		for ( auto const & e: v ) {
			this->push_back( deep_copy( *e ) ); // Clone each entry.
		}
	}

	/// @brief Move constructor
	/// We can steal the values
	deep_copy_vector1( deep_copy_vector1 && v ):
		vector1<T>( std::move(v) )
	{}

	/// @brief Construct from vector1
	/// Right now the semantics are to clone, but I'm not sure if that's the desired effect.
	deep_copy_vector1( vector1< T > const & v ) :
		vector1<T>()
	{
		this->reserve( v.size() );
		for ( auto const & e: v ) {
			this->push_back( deep_copy( *e ) ); // Clone each entry.
		}
	}

	/// @brief Construct from other vectorL
	/// Right now the semantics are to clone, but I'm not sure if that's the desired effect.
	template< typename super::ssize_type L_, typename T_, typename A_ >
	deep_copy_vector1( vectorL< L_, T_, A_ > const & v ) :
		vector1<T>()
	{
		this->reserve( v.size() );
		for ( auto const & e: v ) {
			this->push_back( deep_copy( *e ) ); // Clone each entry.
		}
	}

	/// @brief std::vector constructor
	/// Right now the semantics are to clone, but I'm not sure if that's the desired effect.
	template< typename T_, typename A_ >
	deep_copy_vector1( std::vector< T_, A_ > const & v ) :
		vector1<T>()
	{
		this->reserve( v.size() );
		for ( auto const & e: v ) {
			this->push_back( deep_copy( *e ) ); // Clone each entry.
		}
	}

	/// @brief Destructor
	~deep_copy_vector1() override = default;

public: // Assignment

	/// @brief Copy assignment
	deep_copy_vector1 &
	operator =( deep_copy_vector1 const & v )
	{
		if ( this == &v ) { return *this; }
		this->clear();
		this->reserve( v.size() );
		for ( auto const & e: v ) {
			this->push_back( deep_copy( *e ) ); // Clone each entry.
		}
		return *this;
	}

	/// @brief Move assignment
	/// (We can steal most of the values)
	deep_copy_vector1 &
	operator =( deep_copy_vector1 && v )
	{
		super::operator=( std::move(v) );
		return *this;
	}

	/// @brief General VectorL assignment
	/// Right now the semantics are to clone, but I'm not sure if that's the desired effect.
	template< typename super::ssize_type L_, typename T_, typename A_ >
	deep_copy_vector1 &
	operator =( vectorL< L_, T_, A_ > const & v )
	{
		this->clear();
		this->reserve( v.size() );
		for ( auto const & e: v ) {
			this->push_back( deep_copy( *e ) ); // Clone each entry.
		}
		return *this;
	}

	/// @brief std::vector assignment
	/// Right now the semantics are to clone, but I'm not sure if that's the desired effect.
	template< typename T_, typename A_ >
	deep_copy_vector1 &
	operator =( std::vector< T_, A_ > const & v )
	{
		this->clear();
		this->reserve( v.size() );
		for ( auto const & e: v ) {
			this->push_back( deep_copy( *e ) ); // Clone each entry.
		}
		return *this;
	}

};


} // namespace utility


#endif // INCLUDED_utility_vector1_HH
