// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/chemical/AtomProperties.fwd.hh
/// @brief   Forward declarations for AtomProperties.
/// @author  Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_core_chemical_AtomProperties_FWD_HH
#define INCLUDED_core_chemical_AtomProperties_FWD_HH

// Utility header
#include <utility/pointer/owning_ptr.hh>


namespace core {
namespace chemical {

/// @brief  A class for storing atom properties, stereochemistries, and classifications.
class AtomProperties;

typedef utility::pointer::shared_ptr< AtomProperties > AtomPropertiesOP;
typedef utility::pointer::shared_ptr< AtomProperties const > AtomPropertiesCOP;

}  // namespace chemical
}  // namespace core

#endif  // INCLUDED_core_chemical_AtomProperties_FWD_HH
