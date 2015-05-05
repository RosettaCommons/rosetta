// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/conformation/Atom.fwd.hh
/// @brief  Class declaration for conformation::Atom
/// @note   not to be confused with chemical::Atom
/// @author Phil Bradley


#ifndef INCLUDED_core_conformation_Atom_FWD_HH
#define INCLUDED_core_conformation_Atom_FWD_HH

// Utility header
#include <utility/vector1.fwd.hh>


namespace core {
namespace conformation {

/// @brief A simple class with an atom's position and its chemical type.
class Atom;

typedef utility::vector1< Atom > Atoms;

}  // conformation
}  // core

#endif  // INCLUDED_core_conformation_Atom_FWD_HH
