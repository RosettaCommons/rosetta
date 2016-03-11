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

ResidueRanges::ResidueRanges( ResidueSubset const & subset )
: utility::vector1< ResidueRange >()
{
	this->clear();
	core::Size start_interval = 0;
	for ( core::Size resid=1; resid<=subset.size(); ++resid ) {
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
		}
	}
	if ( start_interval != 0 ) {
		push_back( ResidueRange( start_interval, subset.size() ) );
	}
}

ResidueRanges::~ResidueRanges() {}

} //namespace residue_selector
} //namespace select
} //namespace core

