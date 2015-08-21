// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author

#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>

#include <core/pack/rotamer_set/WaterAnchorInfo.hh>

// utility headers
#include <string>

#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace rotamer_set {

/// @details Auto-generated virtual destructor
WaterAnchorInfo::~WaterAnchorInfo() {}

Size
WaterAnchorInfo::anchor_residue() const {
	return anchor_residue_;
}

void
WaterAnchorInfo::anchor_residue( Size const rsd ) {
	anchor_residue_ = rsd;
}

bool
WaterAnchorInfo::attaches_to_residue_type( ResidueType
	const & type ) const {
	return type.aa() == aa_ && type.has( anchor_atom_name_ );
}

Size
WaterAnchorInfo::anchor_atom( ResidueType const & type
) const {
	debug_assert( attaches_to_residue_type( type ) );
	return type.atom_index( anchor_atom_name_ );
}

void
WaterAnchorInfo::anchor_atom( std::string const & name
) {
	anchor_atom_name_ = name;
}

void
WaterAnchorInfo::aa( AA const & aa_in ) {
	aa_ = aa_in;
}

Size
WaterAnchorInfo::nstep() const {
	return nstep_;
}

void
WaterAnchorInfo::nstep( Size const nstep_in ) {
	nstep_ = nstep_in;
}

} // namespace rotamer_set
} // namespace pack
} // namespace core
