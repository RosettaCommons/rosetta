// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/rational_mc/RationalMonteCarlo.fwd.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_protocols_simple_moves_rational_mc_RationalMonteCarlo_FWD_HH
#define INCLUDED_protocols_simple_moves_rational_mc_RationalMonteCarlo_FWD_HH

#include <utility/pointer/owning_ptr.fwd.hh>

namespace protocols {
namespace simple_moves {
namespace rational_mc {

class RationalMonteCarlo;
typedef utility::pointer::shared_ptr<RationalMonteCarlo> RationalMonteCarloOP;
typedef utility::pointer::shared_ptr<RationalMonteCarlo const> RationalMonteCarloCOP;

}  // namespace rational_mc
}  // namespace simple_moves
}  // namespace protocols

#endif  // protocols_simple_moves_rational_mc_RationalMonteCarlo_FWD_HH_
