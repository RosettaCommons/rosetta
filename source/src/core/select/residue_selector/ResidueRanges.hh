// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/ResidueRanges.hh
/// @brief  The ResidueRanges class identifies a subset of residues from a Pose
/// @author Tom Linsky (tlinsky at uw dot edu)

#ifndef INCLUDED_core_select_residue_selector_ResidueRanges_HH
#define INCLUDED_core_select_residue_selector_ResidueRanges_HH

// Unit headers
#include <core/select/residue_selector/ResidueRanges.fwd.hh>

// Core headers
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/types.hh>

// Basic headers

// Utility Headers
#include <utility/vector1.hh>

// C++ headers

namespace core {
namespace select {
namespace residue_selector {

/// @brief Class for storing a range of resids
class ResidueRange {
public:
	ResidueRange():
		start_( 0 ), stop_( 0 ) {}

	ResidueRange( Size const startval, Size const stopval ):
		start_( startval ), stop_( stopval ) {}

	/// @brief less than operator that can be used for sorting
	/// @details If start < other.start, return true
	///          If start > other.start, return false
	///          if start == other.start, return (stop < other.stop)
	bool
	operator<( ResidueRange const & other ) const;

	/// @brief Returns the starting residue of the range
	Size
	start() const;

	/// @brief Sets starting residue of the range
	void
	set_start( Size const start_res );

	/// @brief Returns the stopping residue of the range
	Size
	stop() const;

	/// @brief Sets stopping residue of the range
	void
	set_stop( Size const stop_res );

private:
	Size start_;
	Size stop_;
};

class ResidueRanges : public utility::vector1< ResidueRange > {
public:

	/// @brief Constructs an empty vector of ResidueRanges
	ResidueRanges();

	/// @brief Constructs a set of contiguous ranges of residues from a residue subset
	/// @param subset : residue subset from which contiguous ranges of residues will be derived
	/// @details ResidueRanges created via this constructor are guaranteed
	///          to be ordered in ascending order by start resid and
	///          contain no duplicates
	ResidueRanges( ResidueSubset const & subset );

	/// @brief Destructor.
	virtual ~ResidueRanges();

	/// @brief Clears the ranges and uses the provided ResidueSubset to create new ranges
	/// @param subset : residue subset from which contiguous ranges of residues will be derived
	/// @details ResidueRanges created via this method are guaranteed
	///          to be ordered in ascending order by start resid and
	///          contain no duplicates
	void
	from_subset( ResidueSubset const & subset );

private:

};


} //namespace residue_selector
} //namespace select
} //namespace core

#endif
