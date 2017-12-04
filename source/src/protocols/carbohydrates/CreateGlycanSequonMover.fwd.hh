// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/carbohydrates/CreateGlycanSequonMover.fwd.hh
/// @brief Mutates residues to create a potential glycosylation site using known sequence motifs of N- or C- linked glycans.  Includes options for Enhanced Sequons for N-linked glycans that have been shown to have higher rates of glycosylation as well as other positions that have been shown to influence the glycosylation chemistry.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_carbohydrates_CreateGlycanSequonMover_fwd_hh
#define INCLUDED_protocols_carbohydrates_CreateGlycanSequonMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace carbohydrates {

class CreateGlycanSequonMover;

typedef utility::pointer::shared_ptr< CreateGlycanSequonMover > CreateGlycanSequonMoverOP;
typedef utility::pointer::shared_ptr< CreateGlycanSequonMover const > CreateGlycanSequonMoverCOP;

} //protocols
} //carbohydrates

#endif //INCLUDED_protocols_carbohydrates_CreateGlycanSequonMover_fwd_hh
