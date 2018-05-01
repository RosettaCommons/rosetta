// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
/// @file   core/pack/interaction_graph/AnnealableGraphBase.cc
/// @brief  Base interface for annealable graphs.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Unit headers
#include <core/pack/interaction_graph/AnnealableGraphBase.hh>

// Package headers

// Project headers

// Utility headers
#include <basic/Tracer.hh>

#include <utility/vector1.hh>


static basic::Tracer TR( "core.pack.interaction_graph.AnnealableGraphBase" );


namespace core {
namespace pack {
namespace interaction_graph {

/// @brief Constructor.
///
AnnealableGraphBase::AnnealableGraphBase() = default;
//TODO -- initialize here

/// @brief Copy constructor.
///
AnnealableGraphBase::AnnealableGraphBase( AnnealableGraphBase const &/*src*/ ) {}
//TODO -- initialize here

/// @brief Destructor.
///
AnnealableGraphBase::~AnnealableGraphBase() = default;

/// @brief Provide the opportunity for an AnnealableGraph to clean up cached data in the pose or inside itself after packing.
/// @details Base class function does nothing; may be overridden in derived classes.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void
AnnealableGraphBase::clean_up_after_packing(
	core::pose::Pose &//pose
) {
	//GNDN
}

} // namespace interaction_graph
} // pack
} // core
