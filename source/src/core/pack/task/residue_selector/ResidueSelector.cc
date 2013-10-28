// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/residue_selector/ResidueSelector.cc
/// @brief  The ResidueSelector class identifies a subset of residues from a Pose
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

// Unit headers
#include <core/pack/task/residue_selector/ResidueSelector.hh>

// Utility headers
#include <utility/tag/Tag.hh>

namespace core {
namespace pack {
namespace task {
namespace residue_selector {

ResidueSelector::ResidueSelector() {}

ResidueSelector::~ResidueSelector() {}

/// @details Noop implementation in the base class in the case that a derived
/// class has no need to read data from an input tag
void ResidueSelector::parse_my_tag(
	utility::tag::TagCOP ,
	basic::datacache::DataMap &
)
{}

} //namespace residue_selector
} //namespace task
} //namespace pack
} //namespace core
