// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/moves/MixedMonteCarlo.fwd.hh
/// @brief A hybrid monte carlo mover
/// @author AmeyaHarmalkar (harmalkar.ameya24@gmail.com)

#ifndef INCLUDED_protocols_moves_MixedMonteCarlo_fwd_hh
#define INCLUDED_protocols_moves_MixedMonteCarlo_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace moves {

class MixedMonteCarlo;

using MixedMonteCarloOP = utility::pointer::shared_ptr< MixedMonteCarlo >;
using MixedMonteCarloCOP = utility::pointer::shared_ptr< MixedMonteCarlo const >;

} //moves
} //protocols

#endif //INCLUDED_protocols_moves_MixedMonteCarlo_fwd_hh
