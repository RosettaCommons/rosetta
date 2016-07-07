// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/select/residue_selector/ResidueRanges.cc
/// @brief  The ResidueRanges class identifies contiguous segments of a pose from a subset
/// @author Tom Linsky (tlinsky at uw dot edu)

// Unit headers
#include <core/select/residue_selector/ResidueRanges.hh>

// Core headers
#include <core/select/residue_selector/ResidueSelector.hh>

// Utility headers

namespace core {
namespace select {
namespace residue_selector {

/// @brief Returns the starting residue of the range
Size
ResidueRange::start() const
{
	return start_;
}

/// @brief Returns the stopping residue of the range
Size
ResidueRange::stop() const
{
	return stop_;
}

/// @brief less than operator that can be used for sorting
/// @details If start < other.start, return true
///          If start > other.start, return false
///          if start == other.start, return (stop < other.stop)
bool
ResidueRange::operator<( ResidueRange const & other ) const
{
	if ( start_ == other.start() ) return ( stop_ < other.stop() );
	return ( start_ < other.start() );
}

/// @brief Constructs an empty vector of ResidueRanges
ResidueRanges::ResidueRanges():
	utility::vector1< ResidueRange >()
{}

/// @brief Constructs a set of contiguous ranges of residues from a residue subset
/// @param subset : residue subset from which contiguous ranges of residues will be derived
ResidueRanges::ResidueRanges( ResidueSubset const & subset ):
	utility::vector1< ResidueRange >()
{
	from_subset( subset );
}

ResidueRanges::~ResidueRanges() {}

/// @brief Clears the ranges and uses the provided ResidueSubset to create new ranges
/// @param subset : residue subset from which contiguous ranges of residues will be derived
void
ResidueRanges::from_subset( ResidueSubset const & subset )
{
	this->clear();
	Size start_interval = 0;
	for ( Size resid=1; resid<=subset.size(); ++resid ) {
		if ( ( start_interval == 0 ) && ( subset[ resid ] ) ) {
			start_interval = resid;
		} else if ( ( start_interval == 0 ) && ( ! subset[ resid ] ) ) {
			continue;
		} else if ( ( start_interval != 0 ) && ( subset[ resid ] ) ) {
			continue;
		} else if ( ( start_interval != 0 ) && ( ! subset[ resid ] ) ) {
			push_back( ResidueRange( start_interval, resid - 1 ) );
			start_interval = 0;
		} else {
			// We should never be here!
			debug_assert( false );
		}
	}
	if ( start_interval != 0 ) {
		push_back( ResidueRange( start_interval, subset.size() ) );
	}
	std::sort( begin(), end() );
}

} //namespace residue_selector
} //namespace select
} //namespace core

