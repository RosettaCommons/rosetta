// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/signals/RestypeDestructionEvent.hh
/// @brief  signal that the Pose is getting destroyed
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_pose_signals_RestypeDestructionEvent_hh
#define INCLUDED_core_pose_signals_RestypeDestructionEvent_hh


// unit headers
#include <core/chemical/RestypeDestructionEvent.fwd.hh>

// package headers
#include <core/chemical/ResidueType.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace chemical {


/// @brief special signal that the ResidueType is getting destroyed
struct RestypeDestructionEvent {

	/// @brief constructor
	inline
	RestypeDestructionEvent( ResidueType const * rt = nullptr ) :
		restype( rt )
	{}

	/// @brief copy constructor
	inline
	RestypeDestructionEvent( RestypeDestructionEvent const & ) = default;

	/// @brief default destructor
	inline
	virtual
	~RestypeDestructionEvent() = default;


	/// @brief copy assignment
	inline
	RestypeDestructionEvent &
	operator =( RestypeDestructionEvent const & /*rval*/ ) = default;

	/// @brief the ResidueType firing the signal
	ResidueType const * restype;

};


} // namespace chemical
} // namespace core


#endif
