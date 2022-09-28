// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/modifications/ChemistryBase.hh
/// @brief  Base class for chemical modifications such as add atoms, delete atoms, replace atom, hydrogen manipulations, etc.
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <core/chemical/modifications/ChemistryBase.hh>
#include <core/chemical/AtomRefMapping.hh>

namespace core {
namespace chemical {
namespace modifications {

/// @brief Get the vertex mapping that was used for the last apply()
/// The base class implementation defaults to an identity mapping.
VDVDMapping
ChemistryBase::get_mapping() const {
	return VDVDMapping(true);
}


}
}
}
