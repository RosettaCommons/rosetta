// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/design_opt/GreedyOptMutationMover.fwd.hh
/// @brief
/// @author Chris King (chrisk1@uw.edu)

#ifndef INCLUDED_protocols_design_opt_GreedyOptMutationMover_fwd_hh
#define INCLUDED_protocols_design_opt_GreedyOptMutationMover_fwd_hh


// Utility headers
#include <utility/pointer/owning_ptr.fwd.hh>

namespace protocols {
namespace design_opt {

class GreedyOptMutationMover;
typedef utility::pointer::shared_ptr< GreedyOptMutationMover >  GreedyOptMutationMoverOP;
typedef utility::pointer::shared_ptr< GreedyOptMutationMover const >  GreedyOptMutationMoverCOP;


} // namespace design_opt
} // namespace protocols

#endif
