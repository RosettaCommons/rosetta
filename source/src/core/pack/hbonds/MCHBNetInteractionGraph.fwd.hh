// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/hbonds/MCHBNetInteractionGraph.fwd.hh
/// @brief Dervied class of PDInteractionGraph that does not save twobody energy calculations but rather passes them directly to a HBondGraph
/// @details This is a HBondGraph creator that is wearing an InteractionGraph disguise so that monte carlo HBNet can collect energy information without having to create custom interfaces in many other classes. This class should not be used as an InteractionGraph because it does not store all of the information that InteractionGraphs need to store. There are a few utility_exit_with_message() calls sprinkled within this class to make sure it is not being misused, but there really is not any need to use it for anything other than HBondGraph creation.
/// @author Jack Maguire, jack@med.unc.edu


#ifndef INCLUDED_core_pack_hbonds_MCHBNetInteractionGraph_FWD_HH
#define INCLUDED_core_pack_hbonds_MCHBNetInteractionGraph_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pack {
namespace hbonds {

class MCHBNetInteractionGraph;
typedef utility::pointer::shared_ptr< MCHBNetInteractionGraph > MCHBNetInteractionGraphOP;
typedef utility::pointer::shared_ptr< MCHBNetInteractionGraph const > MCHBNetInteractionGraphCOP;

} //hbonds
} //pack
} //core

#endif
