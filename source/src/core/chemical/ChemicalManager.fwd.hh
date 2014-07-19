// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/ChemicalManager.hh
/// @brief  Chemical manager class header
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_chemical_ChemicalManager_fwd_hh
#define INCLUDED_core_chemical_ChemicalManager_fwd_hh

// Unit headers

// Package headers

// C++
// Commented by inclean daemon #include <map>
// Commented by inclean daemon #include <string>
#include <string>

namespace core {
namespace chemical {

//singleton class

class ChemicalManager;

extern std::string const FA_STANDARD;
extern std::string const CENTROID;
extern std::string const CENTROID_ROT;
extern std::string const COARSE_TWO_BEAD;
extern std::string const HYBRID_FA_STANDARD_CENTROID;
extern std::string const FA_RNA;
extern std::string const COARSE_RNA;

} // namespace core
} // namespace chemical


#endif
