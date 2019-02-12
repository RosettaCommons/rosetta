// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


/// @file protocols/tcr/TCRloopRefine.fwd.hh
/// @brief TCR cdr loop modeling/refinement
/// @author Ragul Gowthaman (ragul@umd.edu)

// utility headers
#ifndef INCLUDED_protocols_tcr_TCRloopRefine_fwd_hh
#define INCLUDED_protocols_tcr_TCRloopRefine_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace protocols {
namespace tcr {

class TCRloopRefine;

typedef utility::pointer::shared_ptr< TCRloopRefine > TCRloopRefineOP;
typedef utility::pointer::shared_ptr< TCRloopRefine const > TCRloopRefineCOP;
typedef utility::pointer::weak_ptr< TCRloopRefine > TCRloopRefineAP;
typedef utility::pointer::weak_ptr< TCRloopRefine const > TCRloopRefineCAP;

} // tcr
} // protocols

#endif
