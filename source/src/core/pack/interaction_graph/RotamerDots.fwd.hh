// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/interaction_graph/RotamerDots.fwd.hh
/// @brief  RotamerDots classes, forward declaration - ported from rosetta++
/// @author Ron Jacak

#ifndef INCLUDED_core_pack_interaction_graph_RotamerDots_fwd_hh
#define INCLUDED_core_pack_interaction_graph_RotamerDots_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pack {
namespace interaction_graph {

class DotSphere;
class RotamerDots;
class RotamerDotsCache;
class InvRotamerDots;

typedef utility::pointer::shared_ptr< RotamerDots > RotamerDotsOP;
typedef utility::pointer::shared_ptr< RotamerDots const > RotamerDotsCOP;

typedef utility::pointer::shared_ptr< InvRotamerDots > InvRotamerDotsOP;
typedef utility::pointer::shared_ptr< InvRotamerDots const > InvRotamerDotsCOP;


} // interaction_graph
} // pack
} // core


#endif // INCLUDED_core_pack_interaction_graph_RotamerDots_FWD_HH
