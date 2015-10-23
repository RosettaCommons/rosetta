// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/chemical/AtomProperty.hh
/// @brief   Enumeration definitions for AtomProperties.
/// @author  Labonte <JWLabonte@jhu.edu>
/// @note    TODO: Auto-generate this file in the same way ResiduePropertys are generated.


#ifndef INCLUDED_core_chemical_AtomProperty_HH
#define INCLUDED_core_chemical_AtomProperty_HH

namespace core {
namespace chemical {

/// @brief   Enumerators for all the properties that can be assigned to a chemical::Atom.
enum AtomProperty {
	NO_ATOM_PROPERTY = 0,
	FIRST_ATOM_PROPERTY = 1,
	H_DONOR,
	H_ACCEPTOR,
	POLAR_HYDROGEN,
	AROMATIC_HYDROGEN,
	HAS_ORBITALS,
	VIRTUAL_ATOM,
	AROMATIC_CARBON_WITH_FREE_VALENCE,
	N_ATOM_PROPERTIES = AROMATIC_CARBON_WITH_FREE_VALENCE
};

}  // namespace chemical
}  // namespace core

#endif  // INCLUDED_core_chemical_AtomProperty_HH
