// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rna/AlignmentEnergy.fwd.hh
/// @brief An RMSD-type energy to a reference pose, complete with derivatives hacked out of coordinate constraints
/// @author Andy Watkins (amw579@nyu.edu)

#ifndef INCLUDED_protocols_rna_AlignmentEnergy_fwd_hh
#define INCLUDED_protocols_rna_AlignmentEnergy_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace rna {

class AlignmentEnergy;

typedef utility::pointer::shared_ptr< AlignmentEnergy > AlignmentEnergyOP;
typedef utility::pointer::shared_ptr< AlignmentEnergy const > AlignmentEnergyCOP;

} //protocols
} //rna

#endif //INCLUDED_protocols_rna_AlignmentEnergy_fwd_hh
