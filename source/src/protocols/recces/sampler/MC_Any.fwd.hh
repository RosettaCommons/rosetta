// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/recces/sampler/MC_Any.fwd.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_recces_sampler_MC_Any_FWD_HH
#define INCLUDED_protocols_recces_sampler_MC_Any_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace recces {
namespace sampler {

class MC_Any;
typedef utility::pointer::shared_ptr< MC_Any > MC_AnyOP;
typedef utility::pointer::shared_ptr< MC_Any const > MC_AnyCOP;

} //sampler
} //recces
} //protocols

#endif
