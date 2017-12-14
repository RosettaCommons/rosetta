// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/ResidueVector.cc
/// @brief  The ResidueVector class identifies contiguous segments of a pose from a subset
/// @author Tom Linsky (tlinsky at uw dot edu)

// Unit headers
#include <core/select/residue_selector/ResidueVector.hh>

// Core headers
#include <core/select/residue_selector/ResidueSelector.hh>

// Utility headers

namespace core {
namespace select {
namespace residue_selector {

ResidueVector::ResidueVector():
	utility::vector1< core::Size >()
{}

ResidueVector::ResidueVector( ResidueSubset const & subset ):
	utility::vector1< core::Size >()
{
	from_subset( subset );
}

ResidueVector::ResidueVector( utility::vector1< core::Size > const & vec ):
	utility::vector1< core::Size >( vec )
{}

void
ResidueVector::from_subset( ResidueSubset const & subset )
{
	this->clear();
	for ( core::Size resid=1; resid<=subset.size(); ++resid ) {
		if ( subset[resid] ) {
			push_back( resid );
		}
	}
}

ResidueVector::~ResidueVector() = default;

} //namespace residue_selector
} //namespace select
} //namespace core

