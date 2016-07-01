// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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

struct ResidueRange {
	ResidueRange():
		start( 0 ), stop( 0 ) {}

	ResidueRange( core::Size const startval, core::Size const stopval ):
		start( startval ), stop( stopval ) {}

	core::Size start;
	core::Size stop;
};

class ResidueRanges : public utility::vector1< ResidueRange > {
public:

	/// @brief Constructs an empty vector of ResidueRanges
	ResidueRanges();

	/// @brief Constructs a set of contiguous ranges of residues from a residue subset
	/// @param subset : residue subset from which contiguous ranges of residues will be derived
	ResidueRanges( ResidueSubset const & subset );

	/// @brief Destructor.
	virtual ~ResidueRanges();

	/// @brief Clears the ranges and uses the provided ResidueSubset to create new ranges
	/// @param subset : residue subset from which contiguous ranges of residues will be derived
	void
	from_subset( ResidueSubset const & subset );

private:

};


} //namespace residue_selector
} //namespace select
} //namespace core

#endif
