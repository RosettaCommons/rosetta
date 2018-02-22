// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/farnesyl/InstallFarnesylMover.fwd.hh
/// @brief Modifies a free cysteine residue with a branch of 3 DMA residues (the terpene monomer) to create farnesyl-cysteine
/// @author Andy Watkins (andy.watkins2@gmail.com)

#ifndef INCLUDED_protocols_farnesyl_InstallFarnesylMover_fwd_hh
#define INCLUDED_protocols_farnesyl_InstallFarnesylMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace farnesyl {

class InstallFarnesylMover;

typedef utility::pointer::shared_ptr< InstallFarnesylMover > InstallFarnesylMoverOP;
typedef utility::pointer::shared_ptr< InstallFarnesylMover const > InstallFarnesylMoverCOP;

} //protocols
} //farnesyl

#endif //INCLUDED_protocols_farnesyl_InstallFarnesylMover_fwd_hh
