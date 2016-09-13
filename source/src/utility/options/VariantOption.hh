// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/options/VariantOption.hh
/// @brief  Variant option class
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_utility_options_VariantOption_hh
#define INCLUDED_utility_options_VariantOption_hh


// Unit headers
#include <utility/options/VariantOption.fwd.hh>
#include <utility/exit.hh>

// C++ headers
#include <utility/assert.hh>
#include <cstddef>
#include <string>


namespace utility {
namespace options {


/// @brief Variant option class
template< typename O >
class VariantOption
{


public: // Types


	// STL/boost style
	typedef  O  option_type;
	typedef  typename O::key_type  key_type;
	typedef  std::size_t  size_type;

	// Project style
	typedef  O  Option;
	typedef  typename O::Key  Key;
	typedef  std::size_t  Size;


public: // Creation


	/// @brief Default constructor
	inline
	VariantOption() :
		option_p_( 0 )
	{}


	/// @brief Copy constructor
	inline
	VariantOption( VariantOption const & var ) :
		option_p_( var.option_p_ ? var.option_p_->clone() : 0 )
	{}


	/// @brief Option constructor
	inline
	VariantOption( Option const & option_a ) :
		option_p_( option_a.clone() )
	{}


	/// @brief Destructor
	inline
	~VariantOption() throw() // throw() is needed for ICC
	{
		delete option_p_;
	}


public: // Assignment


	/// @brief Copy assignment
	inline
	VariantOption &
	operator =( VariantOption const & var )
	{
		if ( this != &var ) {
			delete option_p_; option_p_ = ( var.option_p_ ? var.option_p_->clone() : 0 );
		}
		return *this;
	}


public: // Conversion


	/// @brief Option conversion
	inline
	operator Option const &() const
	{
		runtime_assert( option_p_ );
		return *option_p_;
	}


	/// @brief Option conversion
	inline
	operator Option &()
	{
		runtime_assert( option_p_ );
		return *option_p_;
	}

	/// @brief Option conversion
	inline
	operator Option const *() const
	{
		runtime_assert( option_p_ );
		return option_p_;
	}


	/// @brief Option conversion
	inline
	operator Option *()
	{
		runtime_assert( option_p_ );
		return option_p_;
	}


public: // Methods


	/// @brief Activate
	inline
	VariantOption &
	activate()
	{
		runtime_assert( option_p_ );
		option_p_->activate();
		return *this;
	}


	/// @brief Deactivate
	inline
	VariantOption &
	deactivate()
	{
		runtime_assert( option_p_ );
		option_p_->deactivate();
		return *this;
	}


	/// @brief Clear
	inline
	VariantOption &
	clear()
	{
		runtime_assert( option_p_ );
		option_p_->clear();
		return *this;
	}


	/// @brief Value assignment from a command line string
	inline
	VariantOption &
	cl_value( std::string const & value_str )
	{
		runtime_assert( option_p_ );
		option_p_->cl_value( value_str );
		return *this;
	}


public: // Properties


	/// @brief Key
	inline
	Key const &
	key() const
	{
		runtime_assert( option_p_ );
		return option_p_->key();
	}


	/// @brief ID
	inline
	std::string const &
	id() const
	{
		runtime_assert( option_p_ );
		return option_p_->id();
	}


	/// @brief Identifier
	inline
	std::string const &
	identifier() const
	{
		runtime_assert( option_p_ );
		return option_p_->identifier();
	}


	/// @brief Code
	inline
	std::string const &
	code() const
	{
		runtime_assert( option_p_ );
		return option_p_->code();
	}


	/// @brief Name
	inline
	std::string const &
	name() const
	{
		runtime_assert( option_p_ );
		return option_p_->name();
	}


	/// @brief Description
	inline
	std::string const &
	description() const
	{
		runtime_assert( option_p_ );
		return option_p_->description();
	}


	/// @brief Active?
	inline
	bool
	active() const
	{
		runtime_assert( option_p_ );
		return option_p_->active();
	}


	/// @brief User-specified?
	inline
	bool
	user() const
	{
		runtime_assert( option_p_ );
		return option_p_->user();
	}


	/// @brief Default size (number of default values)
	inline
	Size
	default_size() const
	{
		runtime_assert( option_p_ );
		return option_p_->default_size();
	}


	/// @brief Number of default values (default size)
	inline
	Size
	n_default_value() const
	{
		runtime_assert( option_p_ );
		return option_p_->n_default_value();
	}


	/// @brief Size (number of values)
	inline
	Size
	size() const
	{
		runtime_assert( option_p_ );
		return option_p_->size();
	}


	/// @brief Number of values (size)
	inline
	Size
	n_value() const
	{
		runtime_assert( option_p_ );
		return option_p_->n_value();
	}


	/// @brief Option
	inline
	Option const &
	operator ()() const
	{
		runtime_assert( option_p_ );
		return *option_p_;
	}


	/// @brief Option
	inline
	Option &
	operator ()()
	{
		runtime_assert( option_p_ );
		return *option_p_;
	}


public: // Comparison


	/// @brief VariantOption < VariantOption
	friend
	inline
	bool
	operator <( VariantOption const & a, VariantOption const & b )
	{
		runtime_assert( a.option_p_ );
		runtime_assert( b.option_p_ );
		return ( *a.option_p_ < *b.option_p_ );
	}


private: // Fields


	/// @brief Pointer to option
	Option * option_p_;


}; // VariantOption


// Friend function namespace declarations


/// @brief VariantOption < VariantOption
template< typename O >
bool
operator <( VariantOption< O > const & a, VariantOption< O > const & b );


} // namespace options
} // namespace utility


#endif // INCLUDED_utility_options_VariantOption_HH
