// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/protein_mpnn/ProteinMPNNMover.fwd.hh
/// @brief  Mover for the ProteinMPNN PyTorch model developed by Dauparas et al.
/// @author Frederick Chan (fredchan@uw.edu)

#ifndef INCLUDED_protocols_protein_mpnn_ProteinMPNN_fwd_hh
#define INCLUDED_protocols_protein_mpnn_ProteinMPNN_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

// Forward
namespace protocols {
namespace protein_mpnn {

class ProteinMPNNMover;

typedef utility::pointer::shared_ptr< ProteinMPNNMover > ProteinMPNNMoverOP;
typedef utility::pointer::shared_ptr< ProteinMPNNMover const > ProteinMPNNMoverCOP;

} // protein_mpnn
} // protocols

#endif //INCLUDED_protocols_protein_mpnn_ProteinMPNNMover_fwd_hh
