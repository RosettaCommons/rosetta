// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/id/DOF_ID_Range.hh
/// @author Colin A. Smith


#ifndef INCLUDED_core_id_DOF_ID_Range_hh
#define INCLUDED_core_id_DOF_ID_Range_hh


// Unit headers
#include <core/id/DOF_ID_Range.fwd.hh>

// Package headers
#include <core/id/types.hh>
#include <core/id/DOF_ID.hh>

// Numeric headers
//#include <numeric/constants.hh>
//#include <numeric/numeric.functions.hh>

// Utility headers
//#include <utility/exit.hh>

// C++ header
#include <utility/assert.hh>

#if (defined min) && (defined WIN32)  // Workaround for MSVC and windows.h include which used #define min
	#undef min
#endif

#if (defined max) && (defined WIN32) // Workaround for MSVC and windows.h include which used #define max
	#undef max
#endif

namespace core {
namespace id {


/// @brief Kinematics DOF identifier (with range) class
class DOF_ID_Range
{

public: // Creation

	DOF_ID_Range(){};

	DOF_ID_Range(
		DOF_ID const & dof_id,
		Real const & min,
		Real const & max
	):
		dof_id_( dof_id ),
		min_( min ),
		max_( max )
	{};

public: // Properties

	inline
	DOF_ID const &
	dof_id() const { return dof_id_; }

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
		DOF_ID_Range const & a
	);

	friend
	inline
	bool
	operator ==(
		DOF_ID_Range const & a,
		DOF_ID_Range const & b
	) { return a.dof_id_ == b.dof_id_ && a.min_ == b.min_ && a.max_ == b.max_; }

	friend
	inline
	bool
	operator !=(
		DOF_ID_Range const & a,
		DOF_ID_Range const & b
	) { return a.dof_id_ != b.dof_id_ || a.min_ != b.min_ || a.max_ != b.max_; }

	friend
	inline
	bool
	operator <(
		DOF_ID_Range const & a,
		DOF_ID_Range const & b
	)
	{
		return ( a.dof_id_ < b.dof_id_ ||
			( ( a.dof_id_ == b.dof_id_ && a.min_ < b.min_ ) ||
			( a.min_ == b.min_ && a.max_ < b.max_ ) ) );
	}

private: // Fields

	/// @brief DOF identifier
	DOF_ID dof_id_;

	/// @brief minimum value
	core::Real min_;

	/// @brief maximum value
	core::Real max_;

}; // DOF_ID_Range


} // namespace id
} // namespace core


#endif // INCLUDED_core_id_DOF_ID_Range_HH
