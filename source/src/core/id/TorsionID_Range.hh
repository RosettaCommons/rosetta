// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/id/TorsionID_Range.hh
/// @author Colin A. Smith


#ifndef INCLUDED_core_id_TorsionID_Range_hh
#define INCLUDED_core_id_TorsionID_Range_hh


// Unit headers
#include <core/id/TorsionID_Range.fwd.hh>

// Package headers
#include <core/id/types.hh>
#include <core/id/TorsionID.hh>

// Numeric headers
//#include <numeric/constants.hh>
//#include <numeric/numeric.functions.hh>

// Utility headers
//#include <utility/exit.hh>

// C++ header
// AUTO-REMOVED #include <utility/assert.hh>
#include <iostream>


namespace core {
namespace id {


/// @brief Kinematics DOF identifier (with range) class
class TorsionID_Range
{

public: // Creation

	TorsionID_Range(){};

	TorsionID_Range(
		TorsionID const & torsion_id,
		Real const & min,
		Real const & max
	):
		torsion_id_( torsion_id ),
		min_( min ),
		max_( max )
	{};

public: // Properties

	inline
	TorsionID const &
	torsion_id() const { return torsion_id_; }

	inline
	Real
	min() const { return min_; }

	inline
	Real
	max() const { return max_; }

public: // Friends

	friend
	std::ostream &
	operator <<(
		std::ostream & os,
		TorsionID_Range const & a
	);

	friend
	inline
	bool
	operator ==(
		TorsionID_Range const & a,
		TorsionID_Range const & b
	) { return a.torsion_id_ == b.torsion_id_ && a.min_ == b.min_ && a.max_ == b.max_; }

	friend
	inline
	bool
	operator !=(
		TorsionID_Range const & a,
		TorsionID_Range const & b
	) { return a.torsion_id_ != b.torsion_id_ || a.min_ != b.min_ || a.max_ != b.max_; }

	friend
	inline
	bool
	operator <(
		TorsionID_Range const & a,
		TorsionID_Range const & b
	)
	{
		return ( a.torsion_id_ < b.torsion_id_ ||
				 ( ( a.torsion_id_ == b.torsion_id_ && a.min_ < b.min_ ) ||
				 ( a.min_ == b.min_ && a.max_ < b.max_ ) ) );
	}

private: // Fields

	/// @brief DOF identifier
	TorsionID torsion_id_;

	/// @brief minimum value
	core::Real min_;

	/// @brief maximum value
	core::Real max_;

}; // TorsionID_Range


} // namespace id
} // namespace core


#endif // INCLUDED_core_id_TorsionID_Range_HH
