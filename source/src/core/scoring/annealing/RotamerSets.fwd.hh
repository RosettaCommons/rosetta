// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/annealing/RotamerSets.fwd.hh
/// @brief  Residue Sets class foward declarations.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
/// @note This forward declaration is one library level down from where it should be.  This is because the ResidueArrayAnnealableEnergy base class, which needs
/// to be in core.3, must define the set_up_residuearrayannealableenergy_for_packing() method, which has a RotamerSets const reference as a parameter, but the
/// RotamersSets class isn't defined until core.4.  This is a minimal, necessary violation of the usual coding conventions that should not appreciably affect
/// compile times.  Since the file core/pack/rotamer_set/RotamerSets.fwd.hh includes this file, usual inclusion of forward declarations ought not to be affected,
/// either.


#ifndef INCLUDED_core_scoring_annealing_RotamerSets_fwd_hh
#define INCLUDED_core_scoring_annealing_RotamerSets_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pack {
namespace rotamer_set {

class RotamerSets;

typedef utility::pointer::shared_ptr< RotamerSets > RotamerSetsOP;
typedef utility::pointer::shared_ptr< RotamerSets const > RotamerSetsCOP;

} // namespace rotamer_set
} // namespace pack
} // namespace core


#endif // INCLUDED_core_scoring_annealing_RotamerSets_fwd_hh
