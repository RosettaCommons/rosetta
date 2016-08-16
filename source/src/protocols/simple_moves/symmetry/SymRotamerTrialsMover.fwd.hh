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

#ifndef INCLUDED_protocols_simple_moves_symmetry_SymRotamerTrialsMover_fwd_hh
#define INCLUDED_protocols_simple_moves_symmetry_SymRotamerTrialsMover_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace simple_moves {
namespace symmetry {

class SymRotamerTrialsMover;
typedef utility::pointer::shared_ptr< SymRotamerTrialsMover > SymRotamerTrialsMoverOP;
typedef utility::pointer::shared_ptr< SymRotamerTrialsMover const > SymRotamerTrialsMoverCOP;

class SymEnergyCutRotamerTrialsMover;
typedef utility::pointer::shared_ptr< SymEnergyCutRotamerTrialsMover > SymEnergyCutRotamerTrialsMoverOP;
typedef utility::pointer::shared_ptr< SymEnergyCutRotamerTrialsMover const > SymEnergyCutRotamerTrialsMoverCOP;

} // symmetry
} // moves
} // protocols

#endif
