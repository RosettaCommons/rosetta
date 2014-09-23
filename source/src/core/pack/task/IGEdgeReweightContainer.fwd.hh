// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file Interaction graph edge reweighting
/// @brief
/// @author Florian Richter, floric@u.washington.edu

#ifndef INCLUDED_core_pack_task_IGEdgeReweightContainer_fwd_hh
#define INCLUDED_core_pack_task_IGEdgeReweightContainer_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pack {
namespace task {

class IGEdgeReweighter;
class IGEdgeReweightContainer;

typedef utility::pointer::shared_ptr< IGEdgeReweighter > IGEdgeReweighterOP;
typedef utility::pointer::shared_ptr< IGEdgeReweighter const > IGEdgeReweighterCOP;

typedef utility::pointer::shared_ptr< IGEdgeReweightContainer > IGEdgeReweightContainerOP;
typedef utility::pointer::shared_ptr< IGEdgeReweightContainer const > IGEdgeReweightContainerCOP;

} //task
} //pack
} //core

#endif
