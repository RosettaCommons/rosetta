// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/id/AtomID.hh
/// @brief  Declarations for the TorsionID class
/// @author Phil Bradley


#ifndef INCLUDED_core_id_TorsionID_hh
#define INCLUDED_core_id_TorsionID_hh

// Unit headers
#include <core/id/TorsionID.fwd.hh>

// Package headers
#include <core/id/types.hh>

// C++ headers
#include <iostream>


namespace core {
namespace id {


/// @brief  Torsion identifier class
/// @note   Note the order of the elements in construction calls:
///
///        ( residue, type, torsion )
///
///        sort of from least to most specific
class TorsionID
{
public: // Creation

	/// @brief Default constructor
	inline
	TorsionID() :
		rsd_( 0 ),
		type_( BB ),
		torsion_( 0 )
	{};

	/// @brief Copy constructor
	inline
	TorsionID( TorsionID const & src ) :
		rsd_( src.rsd_ ),
		type_( src.type_ ),
		torsion_( src.torsion_ )
	{}

	/// @brief Property constructor
	inline
	TorsionID(
		Size const rsd_in,
		TorsionType const & type_in,
		Size const torsion_in
	) :
		rsd_( rsd_in ),
		type_( type_in ),
		torsion_( torsion_in )
	{}

public: // Properties

	inline
	Size
	rsd() const { return rsd_; }

	inline
	Size &
	rsd() { return rsd_; }

	inline
	TorsionType
	type() const { return type_; }

	inline
	TorsionType &
	type() { return type_; }

	inline
	uint
	torsion() const { return torsion_; }

	inline
	Size &
	torsion() { return torsion_; }

	/// @brief Is this id valid?
	/// \note Must return false for BOGUS_TORSION_ID
	inline
	bool
	valid() const { return ( rsd_ > 0 && torsion_ > 0 ); }

public: // Friends

 	friend
	inline
 	std::ostream &
 	operator <<(
 		std::ostream & os,
 		TorsionID const & a
 	)
	{
		os << "TorsionID " << a.rsd_ << ' ' << a.type_ << ' ' << a.torsion_;
		return os;
	}

	friend
	inline
	bool
	operator ==(
		TorsionID const & a,
		TorsionID const & b
	)
	{
		return a.type_ == b.type_ && a.torsion_ == b.torsion_ && a.rsd_ == b.rsd_;
	}

	friend
	inline
	bool
	operator !=(
		TorsionID const & a,
		TorsionID const & b
	)
	{
		return a.type_ != b.type_ || a.torsion_ != b.torsion_ || a.rsd_ != b.rsd_;
	}

	friend
	inline
	bool
	operator <(
		TorsionID const & a,
		TorsionID const & b
	)
	{
		if      ( a.rsd_     < b.rsd_     ) return true;
		else if ( a.rsd_     > b.rsd_     ) return false;

		if      ( a.type_    < b.type_    ) return true;
		else if ( a.type_    > b.type_    ) return false;

		if      ( a.torsion_ < b.torsion_ ) return true;

		return false;
	}

private: // Fields


	/// @brief Residue number within the complex
	Size rsd_;

	/// @brief  The type (BB,CHI,JUMP) of this torsion
	TorsionType type_;

	/// @brief Torsion number of the given type within the residue
	Size torsion_;


}; // TorsionID


/// @brief Globals
extern TorsionID const BOGUS_TORSION_ID;


} // namespace id
} // namespace core


#endif // INCLUDED_core_id_TorsionID_HH
