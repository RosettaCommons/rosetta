// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/id/DOF_ID.hh
/// @author Phil Bradley


#ifndef INCLUDED_core_id_DOF_ID_hh
#define INCLUDED_core_id_DOF_ID_hh


// Unit headers
#include <core/id/DOF_ID.fwd.hh>

// Package headers
#include <core/id/types.hh>
#include <core/id/AtomID.hh>

// C++ header
#include <cassert>


namespace core {
namespace id {


/// @brief Kinematics DOF identifier class
class DOF_ID
{

public: // Creation

	DOF_ID(){};

	DOF_ID(
		AtomID const & atom_id_in,
		DOF_Type const & type_in
	):
		atom_id_( atom_id_in ),
		type_( type_in )
	{};

public: // Properties

	inline
	AtomID const &
	atom_id() const { return atom_id_; }

	inline
	Size
	rsd() const { return atom_id_.rsd(); }

	inline
	Size
	atomno() const { return atom_id_.atomno(); }

	inline
	DOF_Type
	type() const { return type_; }


	/// @brief Is this DOF_ID valid?
	/// @note Must return false for BOGUS_TORSION_ID.
	inline
	bool
	valid() const { return atom_id_.valid(); }

public: // Friends

	friend
	std::ostream &
	operator <<(
		std::ostream & os,
		DOF_ID const & a
	);

	friend
	inline
	bool
	operator ==(
		DOF_ID const & a,
		DOF_ID const & b
	) { return a.atom_id_ == b.atom_id_ && a.type_ == b.type_; }

	friend
	inline
	bool
	operator !=(
		DOF_ID const & a,
		DOF_ID const & b
	) { return a.atom_id_ != b.atom_id_ || a.type_ != b.type_; }

	friend
	inline
	bool
	operator <(
		DOF_ID const & a,
		DOF_ID const & b
	)
	{
		return ( a.atom_id_ < b.atom_id_ ||
				 ( a.atom_id_ == b.atom_id_ && a.type_ < b.type_ ) );
	}

private: // Fields

	/// @brief Atom identifier
	AtomID atom_id_;

	/// @brief DOF type
	DOF_Type type_;

}; // DOF_ID


/// @brief Globals
extern DOF_ID const BOGUS_DOF_ID;


/////////////////////////////////////////////////////////////////////////////
// @brief helpful conversion
inline
DOF_Type
get_rb_type( Size const k ) {
	assert( k>=1 && k<=6 );
//	return ( k == 1 ? RB1 : ( k==2 ? RB2 : ( k == 3 ? RB3 :
//				 ( k == 4 ? RB4 : ( k==5 ? RB5 : RB6 ) ) ) ) );
	return DOF_Type( RB1 + k - 1 ); //SGM I think this should work and be a little faster: Requires RB's to be contiguous but that seem safe
}


/////////////////////////////////////////////////////////////////////////////
inline
Size
get_rb_number( DOF_Type const t ) {
//	if ( t == PHI || t == THETA || t == D ) return 0;
//	else {
//		return ( t == RB1 ? 1 : ( t == RB2 ? 2 : ( t == RB3 ? 3 :
//					 ( t == RB4 ? 4 : t == RB5 ? 5 : 6 ) ) ) );
//	}
	return ( ( t >= RB1 ) && ( t <= RB6 ) ? t - RB1 + 1 : 0 ); //SGM I think this should work and be a little faster: Requires RB's to be contiguous but that seem safe
}


/////////////////////////////////////////////////////////////////////////////
inline
bool
DOF_type_is_rb( DOF_Type const t )
{
	return ( t >= RB1 && t <= RB6 );
}

} // namespace id
} // namespace core

#endif // INCLUDED_core_id_DOF_ID_HH
