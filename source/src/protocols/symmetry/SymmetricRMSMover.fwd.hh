// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Ingemar Andre


#ifndef INCLUDED_protocols_simple_moves_symmetry_SymmetricRMSMover_fwd_hh
#define INCLUDED_protocols_simple_moves_symmetry_SymmetricRMSMover_fwd_hh

// Unit headers

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace symmetry {

class SymmetricRMSMover;
typedef utility::pointer::shared_ptr< SymmetricRMSMover > SymmetricRMSMoverOP;
typedef utility::pointer::shared_ptr< SymmetricRMSMover const > SymmetricRMSMoverCOP;

} // moves
} // rosetta
#endif
