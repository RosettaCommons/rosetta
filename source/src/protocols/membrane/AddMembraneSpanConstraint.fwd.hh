// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       protocols/membrane/AddMembraneSpanConstraintCreator.hh
/// @brief      add membrane span constraint
/// @author     Jonathan Weinstein (jonahtan.weinstein@weizmann.ac.il)

#ifndef INCLUDED_protocols_membrane_AddMembraneSpanConstraint_fwd_hh
#define INCLUDED_protocols_membrane_AddMembraneSpanConstraint_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace membrane {

class AddMembraneSpanConstraint;
typedef utility::pointer::shared_ptr< AddMembraneSpanConstraint > AddMembraneSpanConstraintOP;
typedef utility::pointer::shared_ptr< AddMembraneSpanConstraint const > AddMembraneSpanConstraintCOP;

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_AddMembraneSpanConstraint_fwd_hh
