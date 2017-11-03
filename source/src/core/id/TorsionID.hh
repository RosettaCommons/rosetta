// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
/// @details Consider a few examples to get a better picture for how torsions
/// are uniquely identified:
///
/// @code
/// #include <core/id/types.hh>
/// using core::id::BB;
///
/// TorsionID(253, BB, 1)  // Phi backbone torsion of residue 253.
/// TorsionID(253, BB, 2)  // Psi backbone torsion of residue 253.
/// TorsionID(253, BB, 3)  // Omega backbone torsion of residue 253.
/// @endcode
///
/// Note the order of the elements in construction calls (residue, type,
/// torsion) go from least to most specific.
///
/// TorsionIDs are very different for JUMP TorsionTypes. In such a case, they are interpreted as follows:
/// TorsionID(1, JUMP, 2)  // RB2 of jump #1 for the Pose.
class TorsionID
{
public: // Creation

	/// @brief Default constructor
	inline
	constexpr
	TorsionID() :
		rsd_( 0 ),
		type_( BB ),
		torsion_( 0 )
	{};

	/// @brief Copy constructor
	inline
	TorsionID( TorsionID const & ) = default;

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

	static constexpr TorsionID BOGUS_TORSION_ID() { return TorsionID(); }

public: // Properties

	/// @brief For this TorsionID, return the Residue number within the complex OR the Jump number for the Pose.
	inline
	Size
	rsd() const { return rsd_; }

	/// @brief For this TorsionID, set the Residue number within the complex OR the Jump number for the Pose.
	inline
	Size &
	rsd() { return rsd_; }

	/// @brief Return the type (BB, CHI,NU, BRANCH, JUMP) of this torsion.
	inline
	TorsionType
	type() const { return type_; }

	/// @brief Set the type (BB, CHI,NU, BRANCH, JUMP) of this torsion.
	inline
	TorsionType &
	type() { return type_; }

	/// @brief Return the torsion number of the given type within the residue OR
	/// the rigid-body identifier for a JUMP TorsionType.
	inline
	uint
	torsion() const { return torsion_; }

	/// @brief Set the torsion number of the given type within the residue OR
	/// the rigid-body identifier for a JUMP TorsionType.
	inline
	Size &
	torsion() { return torsion_; }

	/// @brief Is this id valid?
	/// @note Must return false for BOGUS_TORSION_ID
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
		os << "TorsionID " << a.rsd_ << ' ' << to_string( a.type_ ) << ' ' << a.torsion_;
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

	// Residue number within the complex OR the Jump number for the Pose
	Size rsd_;

	// The type (BB, CHI,NU, BRANCH, JUMP) of this torsion
	TorsionType type_;

	// Torsion number of the given type within the residue OR the rigid-body identifier for a JUMP TorsionType
	Size torsion_;


#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // TorsionID


/// @brief Global -- may not be used until after core::init is called.
extern TorsionID const GLOBAL_BOGUS_TORSION_ID;

} // namespace id
} // namespace core


#endif // INCLUDED_core_id_TorsionID_HH
