// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/chemical/ChemicalManager.hh
/// @brief  Chemical manager class header
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_chemical_ChemicalManager_fwd_hh
#define INCLUDED_core_chemical_ChemicalManager_fwd_hh

// Unit headers

// Package headers

// C++
#include <string>

namespace core {
namespace chemical {

//singleton class

class ChemicalManager;

/// @brief A type set category is the "compatibility class" of a type set.
/// That is, all e.g. ResidueTypes of a given TypeSetCategory should be
/// "compatible" with the scale of modeling resolution, indepenent of if
/// they're in the same ResidueTypeSet.
enum TypeSetCategory {
	/// @brief Don't use this as a category - it exists as an error signal only.
	INVALID_t = 0,
	/// @brief Full atom modeling - all atoms are represented
	FULL_ATOM_t = 1,
	/// @brief For type sets where there really isn't a clear distinction between
	/// fullatom and centroid.
	DEFAULT_t,
	/// @brief Standard Rosetta Centroid: amino acid sidechains are represented by
	/// a "superatom", and other atoms are represented by a unified-atom model.
	/// (No apolar hydrogens.)
	CENTROID_t,
	/// @brief A modification of the Centroid mode, where the amino acid sidechains
	/// have additional degrees of freedom.
	CENTROID_ROT_t,
	/// @brief ???
	HYBRID_FA_STANDARD_CENTROID_t,
	/// @brief ???
	COARSE_RNA_t,
	/// @brief MIXED is not an actual category, but is instead used for a a mixed TypeSetCategory situation
	MIXED_t,
	TYPE_SET_CATEGORIES_LENGTH = MIXED_t // keep at end
};

extern std::string const FA_STANDARD;
extern std::string const CENTROID;
extern std::string const CENTROID_ROT;
extern std::string const HYBRID_FA_STANDARD_CENTROID;
extern std::string const COARSE_RNA;

} // namespace core
} // namespace chemical


#endif
