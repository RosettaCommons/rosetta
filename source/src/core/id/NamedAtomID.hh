// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/id/NamedAtomID.hh
/// @author Oliver Lange


#ifndef INCLUDED_core_id_NamedAtomID_hh
#define INCLUDED_core_id_NamedAtomID_hh

// Unit headers
#include <core/id/NamedAtomID.fwd.hh>

// C++ headers
#ifdef WIN32
#include <string>
#endif

#include <core/types.hh>
#include <ostream>
#include <utility>


namespace core {
namespace id {


/////////////////////////////////////////////////////////////////////////////
// NamedAtomID
/////////////////////////////////////////////////////////////////////////////

// data to uniquely specify an atom
// could change -- pointer?

/// @brief  Atom identifier class
class NamedAtomID
{

public: // Creation

	/// @brief Default constructor
	inline
	NamedAtomID() :
		atom_( "" ),
		rsd_( 0 )
	{}

	/// @brief Copy constructor
	inline
	NamedAtomID( NamedAtomID const & ) = default;

	/// @brief Property constructor
	inline
	NamedAtomID(
		std::string  atom_in,
		Size const rsd_in
	) :
		atom_(std::move( atom_in )),
		rsd_( rsd_in )
	{}

public: // Properties

	inline
	Size
	rsd() const { return rsd_; }

	inline
	Size &
	rsd() { return rsd_; }

	inline
	std::string const&
	atom() const { return atom_; }

	inline
	std::string &
	atom() { return atom_; }

	/// @brief Is this id valid?
	/// @note Must return false for BOGUS_ATOM_ID
	inline
	bool
	valid() const { return ( atom_.size() ) && ( rsd() > 0 ); }

	std::string to_string() const;

public: // Friends

	friend
	std::ostream &
	operator <<(
		std::ostream & os,
		NamedAtomID const & a
	);

	/// @brief input operator
	friend std::istream & operator >> ( std::istream & is, NamedAtomID & e );

	/// @brief a and b are the same atom
	friend
	inline
	bool
	operator ==(
		NamedAtomID const & a,
		NamedAtomID const & b
	) { return a.atom_ == b.atom_ && a.rsd_ == b.rsd_; }

	/// @brief a and b are different atom
	friend
	inline
	bool
	operator !=(
		NamedAtomID const & a,
		NamedAtomID const & b
	) { return a.atom_ != b.atom_ || a.rsd_ != b.rsd_; }

	/// @brief a is LOWER than b (e.g., first by smaller residue index number then by smaller atom index number)
	friend
	inline
	bool
	operator <(
		NamedAtomID const & a,
		NamedAtomID const & b
	) {
		return ( a.rsd_ <  b.rsd_ ||
			( a.rsd_ == b.rsd_ && a.atom_ < b.atom_ ) );
	}

private: // Fields
	/// @brief Atom number within the Residue
	std::string atom_;

	/// @brief Residue number within the complex
	Size rsd_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // NamedAtomID


/// @brief Globals
extern NamedAtomID const BOGUS_NAMED_ATOM_ID;
extern NamedAtomID const CHAINBREAK_BOGUS_NAMED_ATOM_ID;

} // id
} // core

#endif
