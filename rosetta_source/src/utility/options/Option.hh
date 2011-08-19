// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/options/Option.hh
/// @brief  Program option interface class
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Modified by Sergey Lyskov (Sergey.Lyskov@jhu.edu)


#ifndef INCLUDED_utility_options_Option_hh
#define INCLUDED_utility_options_Option_hh


// Unit headers
#include <utility/options/Option.fwd.hh>

// Package headers
#include <utility/options/keys/OptionKey.hh>

// C++ headers
#include <cstddef>
#include <string>


namespace utility {
namespace options {


/// @brief Program option interface class
class Option
{


protected: // Types


	/// @brief Option state enumeration
	enum State {
		INACTIVE, // No default or user-specified value
		DEFAULT,  // Default value active
		USER      // User-specified value active
	};


public: // Types


	// STL/boost style
	typedef  OptionKey  key_type;
	typedef  std::size_t  size_type;

	// Project style
	typedef  OptionKey  Key;
	typedef  std::size_t  Size;


protected: // Creation


	/// @brief Default constructor
	inline
	Option() :  been_accessed_(false)
	{}


	/// @brief Copy constructor
	inline
	Option( Option const & ) :  been_accessed_(false)
	{}


public: // Creation


	/// @brief Clone this
	virtual
	Option *
	clone() const = 0;


	/// @brief Destructor
	inline
	virtual
	~Option()
	{}


protected: // Assignment


	/// @brief Copy assignment
	inline
	Option &
	operator =( Option const & )
	{
		return *this;
	}


public: // Methods


	/// @brief Activate
	virtual
	Option &
	activate() = 0;


	/// @brief Deactivate
	virtual
	Option &
	deactivate() = 0;


	/// @brief Set to default value, if any
	virtual
	Option &
	to_default() = 0;


	/// @brief Clear
	virtual
	Option &
	clear() = 0;


	/// @brief Value assignment from a command line string
	virtual
	Option &
	cl_value( std::string const & value_str ) = 0;


	/// @brief Legal specifications check: Report and return error state
	virtual
	bool
	legal_specs_report() const = 0;


	/// @brief Legal value limits check: Report and return error state
	virtual
	bool
	legal_limits_report() const = 0;


	/// @brief Legal size limits check: Report and return error state
	virtual
	bool
	legal_size_report() const = 0;


	/// @brief Legal default value check: Report and return error state
	virtual
	bool
	legal_default_report() const = 0;


	/// @brief Legal default value check
	virtual
	void
	legal_default_check() const = 0;


	/// @brief Legal value check: Report and return error state
	virtual
	bool
	legal_report() const = 0;


	/// @brief Legal value check
	virtual
	void
	legal_check() const = 0;


	/// @brief Required specified option check: Report and return error state
	virtual
	bool
	specified_report() const = 0;


	/// @brief Required specified option check
	virtual
	void
	specified_check() const = 0;


public: // Properties


	/// @brief Key
	virtual
	Key const &
	key() const = 0;


	/// @brief ID
	virtual
	std::string const &
	id() const = 0;


	/// @brief Identifier
	virtual
	std::string const &
	identifier() const = 0;


	/// @brief Code
	virtual
	std::string const &
	code() const = 0;


	/// @brief Name
	virtual
	std::string const &
	name() const = 0;


	/// @brief Description
	virtual
	std::string const &
	description() const = 0;


	/// @brief short_Description
	virtual
	std::string const &
	short_description() const = 0;


	/// @brief Legal or inactive default value?
	virtual
	bool
	legal_default() const = 0;


	/// @brief Legal value?
	virtual
	bool
	legal() const = 0;


	/// @brief Has a default?
	virtual
	bool
	has_default() const = 0;


	/// @brief Default active?
	virtual
	bool
	default_active() const = 0;


	/// @brief Default inactive?
	virtual
	bool
	default_inactive() const = 0;


	/// @brief Active?  That is, the option has some value, either the default one or specified on the command line.
	virtual
	bool
	active() const = 0;


	/// @brief User-specified?  That is, the option value was specified on the command line.
	/// You should probably use active() instead in almost all cases!
	virtual
	bool
	user() const = 0;


	/// @brief Is a string readable as this option's value type?
	virtual
	bool
	is_value( std::string const & value_str ) const = 0;


	/// @brief Is a string readable as this option's value type and a legal command line value?
	virtual
	bool
	is_cl_value( std::string const & value_str ) const = 0;


	/// @brief Can another value be added and stay within any size constraints?
	virtual
	bool
	can_hold_another() const = 0;


	/// @brief Default size (number of default values)
	virtual
	Size
	default_size() const = 0;


	/// @brief Number of default values (default size)
	virtual
	Size
	n_default_value() const = 0;


	/// @brief Size (number of values)
	virtual
	Size
	size() const = 0;


	/// @brief Number of values (size)
	virtual
	Size
	n_value() const = 0;


	/// @brief Option type code string representation
	virtual
	std::string
	type_string() const = 0;


	/// @brief Legal value string representation
	virtual
	std::string
	legal_string() const = 0;


	/// @brief Size constraint string representation
	virtual
	std::string
	size_constraint_string() const = 0;


	/// @brief Default value string representation
	virtual
	std::string
	default_string() const = 0;


	/// @brief Value string representation
	virtual
	std::string
	value_string() const = 0;


	/// @brief =Value string representation
	virtual
	std::string
	equals_string() const = 0;



	/// @brief Set access property to true.
	void been_accessed() const { been_accessed_ = true; }
	void set_accessed( bool setting ) const { been_accessed_ = setting; }

	/// @brief Return true if option value was anyhow accessed.
	bool is_been_accessed() const { return been_accessed_; }

public: // Comparison


	/// @brief Option < Option
	/// @note  Key-based ordering
	/// @note  Needed for use as option in associative containers
	friend
	inline
	bool
	operator <( Option const & a, Option const & b )
	{
		return ( a.key() < b.key() );
	}



private: // Private data members

	/// @brief flag, will be true if application was trying to anyhow access/check option value.
	///        Used to create option usage reports.
	///        False by default, any access functions ie: user(), active(), value(), operator()() will set it to true.
	mutable bool  been_accessed_;

}; // Option


// Friend function namespace declarations


/// @brief Option < Option
bool
operator <( Option const & a, Option const & b );


} // namespace options
} // namespace utility


#endif // INCLUDED_utility_options_Option_HH
