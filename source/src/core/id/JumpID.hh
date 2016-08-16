// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/id/AtomID.hh
/// @author Phil Bradley


#ifndef INCLUDED_core_id_JumpID_hh
#define INCLUDED_core_id_JumpID_hh


// Unit headers
#include <core/id/JumpID.fwd.hh>

// Package headers


// C++ headers
#include <iostream>

#include <core/types.hh>


namespace core {
namespace id {

/// /brief  Torsion identifier class
///
/// Note the order of the elements in construction calls:
///
///        ( residue, type, torsion )
///
///        sort of from least to most specific

class JumpID
{

public: // Creation

	/// @brief Default constructor
	inline
	JumpID() :
		rsd1_( 0 ),
		rsd2_( 0 )
	{};

	/// @brief Property constructor
	inline
	JumpID(
		Size const rsd1_in,
		Size const rsd2_in
	)
	{
		rsd1_ = rsd1_in < rsd2_in ? rsd1_in : rsd2_in;
		rsd2_ = rsd1_in >= rsd2_in ? rsd1_in : rsd2_in;
	}

public: // Properties

	inline
	Size
	rsd1() const { return rsd1_; }

	inline
	Size
	rsd2() const { return rsd2_; }

	//fd why is this returning a nonconst reference?
	inline
	Size &
	rsd2() { return rsd2_; }

	/// @brief Is this id valid?
	/// \note Must return false for BOGUS_TORSION_ID
	inline
	bool
	valid() const { return ( rsd1_ > 0 && rsd2_ > rsd1_ ); }

public: // Friends

	friend
	inline
	std::ostream &
	operator <<(
		std::ostream & os,
		JumpID const & a
	)
	{
		os << "JumpID " << a.rsd1_ << ' ' << a.rsd2_;
		return os;
	}

	friend
	inline
	bool
	operator ==(
		JumpID const & a,
		JumpID const & b
	)
	{
		return a.rsd1_ == b.rsd1_ && a.rsd2_ == b.rsd2_;
	}

	friend
	inline
	bool
	operator !=(
		JumpID const & a,
		JumpID const & b
	)
	{
		return a.rsd1_ != b.rsd1_ || a.rsd2_ != b.rsd2_;
	}

	friend
	inline
	bool
	operator <(
		JumpID const & a,
		JumpID const & b
	)
	{
		if      ( a.rsd1_ < b.rsd1_     ) return true;
		else if ( a.rsd1_ ==  b.rsd1_ && a.rsd2_ < b.rsd2_  ) return true;

		return false;
	}

private: // Fields


	/// @brief Residue number within the complex
	Size rsd1_;
	Size rsd2_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // JumpID

} // namespace id
} // namespace core


#endif // INCLUDED_core_id_JumpID_HH
