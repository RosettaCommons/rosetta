// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/symmetry/SymmetricEnergies.fwd.hh
/// @brief  Symmetric Energies class to store cached energies and track the residue
/// neighbor relationships
/// @author Ingemar Andre

#ifndef INCLUDED_core_scoring_symmetry_SymmetricEnergies_fwd_hh
#define INCLUDED_core_scoring_symmetry_SymmetricEnergies_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace symmetry {

// Forward
class SymmetricEnergies;

typedef utility::pointer::shared_ptr< SymmetricEnergies       > SymmetricEnergiesOP;
typedef utility::pointer::shared_ptr< SymmetricEnergies const > SymmetricEnergiesCOP;

}	// namespace symmetry
} // namespace pose
} // namespace core

#endif // INCLUDED_core_pose_symmetry_SymmetricEnergies_FWD_HH
