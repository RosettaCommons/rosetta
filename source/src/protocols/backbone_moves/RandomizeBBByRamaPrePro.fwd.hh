// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/backbone_moves/RandomizeBBByRamaPrePro.fwd.hh
/// @brief A simple mover to randomize a backbone, or a portion of a backbone, biased by the rama_prepro score of each residue.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_protocols_backbone_moves_RandomizeBBByRamaPrePro_fwd_hh
#define INCLUDED_protocols_backbone_moves_RandomizeBBByRamaPrePro_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace backbone_moves {

class RandomizeBBByRamaPrePro;

typedef utility::pointer::shared_ptr< RandomizeBBByRamaPrePro > RandomizeBBByRamaPreProOP;
typedef utility::pointer::shared_ptr< RandomizeBBByRamaPrePro const > RandomizeBBByRamaPreProCOP;

} //protocols
} //backbone_moves

#endif //INCLUDED_protocols_backbone_moves_RandomizeBBByRamaPrePro_fwd_hh
