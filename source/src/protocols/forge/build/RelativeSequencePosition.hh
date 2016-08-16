// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/forge/build/RelativeSequencePosition.hh
/// @brief function objects to compute relative positions for RelativeConnectRight
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_protocols_forge_build_RelativeSequencePosition_hh
#define INCLUDED_protocols_forge_build_RelativeSequencePosition_hh

// unit headers
#include <protocols/forge/build/RelativeSequencePosition.fwd.hh>

// Package headers
#include <protocols/forge/build/BuildInstruction.hh>

/// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace forge {
namespace build {


/// @brief computes a position wrt values in given BuildInstruction
struct RelativeSequencePosition : public utility::pointer::ReferenceCount {
	virtual ~RelativeSequencePosition() {}
	virtual RelativeSequencePositionOP clone() const = 0;
	virtual core::Size operator ()( BuildInstructionCAP i ) const = 0;
};


/// @brief count starting from interval().left in a BuildInstruction
struct CountFromLeft : public RelativeSequencePosition {
	inline
	virtual RelativeSequencePositionOP clone() const {
		return RelativeSequencePositionOP( new CountFromLeft( *this ) );
	}

	inline
	virtual core::Size operator ()( BuildInstructionCAP i ) const {
		return i.lock()->interval().left + left_skip + p - 1;
	}

	/// @brief the position to re-compute using the instruction's interval: 1-based indexing
	core::Size p;

	/// @brief additionally skip 'n' positions prior to computing the position
	core::Size left_skip;
};


} // namespace build
} // namespace forge
} // namespace protocols


#endif /* INCLUDED_protocols_forge_build_RelativeSequencePosition_HH */
