// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/drug_design/ApplyChemistryMover.fwd.hh
/// @brief Apply a given Chemistry modifier to the ResidueType at a given position, then replace the ResidueType at that position.
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_protocols_drug_design_ApplyChemistryMover_fwd_hh
#define INCLUDED_protocols_drug_design_ApplyChemistryMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace drug_design {

class ApplyChemistryMover;

typedef utility::pointer::shared_ptr< ApplyChemistryMover > ApplyChemistryMoverOP;
typedef utility::pointer::shared_ptr< ApplyChemistryMover const > ApplyChemistryMoverCOP;

} //protocols
} //drug_design

#endif //INCLUDED_protocols_drug_design_ApplyChemistryMover_fwd_hh
