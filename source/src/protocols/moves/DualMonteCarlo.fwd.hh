// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/moves/DualMonteCarlo.fwd.hh
/// @brief "dual" MonteCarlo header - wraps MonteCarlo to allow for centroid MonteCarlo scoring of a fullatom pose (or other similar situations)
/// @author Steven Lewis

#ifndef INCLUDED_protocols_moves_DualMonteCarlo_fwd_hh
#define INCLUDED_protocols_moves_DualMonteCarlo_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace moves {

class DualMonteCarlo;

typedef utility::pointer::shared_ptr< DualMonteCarlo > DualMonteCarloOP;
typedef utility::pointer::shared_ptr< DualMonteCarlo const > DualMonteCarloCOP;

} // moves
} // protocols

#endif //INCLUDED_protocols_moves_DualMonteCarlo_FWD_HH
