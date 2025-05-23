// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/constraint_generator/TerminiConstraintGenerator.fwd.hh
/// @brief Generates distance constraints between the upper and lower termini
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_constraint_generator_TerminiConstraintGenerator_fwd_hh
#define INCLUDED_protocols_constraint_generator_TerminiConstraintGenerator_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

// Forward
namespace protocols {
namespace constraint_generator {

class TerminiConstraintGenerator;

typedef utility::pointer::shared_ptr< TerminiConstraintGenerator > TerminiConstraintGeneratorOP;
typedef utility::pointer::shared_ptr< TerminiConstraintGenerator const > TerminiConstraintGeneratorCOP;

} //protocols
} //constraint_generator


#endif //INCLUDED_protocols_constraint_generator_TerminiConstraintGenerator_fwd_hh
