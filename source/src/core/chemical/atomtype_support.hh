// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/chemical/atomtype_support.hh
/// @author Rocco Moretti (rmorettiase@gmail.com)


#ifndef INCLUDED_core_chemical_atomtype_support_hh
#define INCLUDED_core_chemical_atomtype_support_hh

#include <core/chemical/ResidueType.hh>

#include <core/chemical/ResidueGraphTypes.hh>

namespace core {
namespace chemical {

bool retype_is_aromatic(VD const & atom, ResidueGraph const & graph);

/// @brief Reassign Rosetta atom types based on the current fullatom heuristics.
///
/// If preserve is true, only retype those atoms which have an atom_type_index of zero.
void rosetta_retype_fullatom(ResidueType & restype, bool preserve=false);

} // chemical
} // core

#endif
