// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/ContextGraphFactory.cc
/// @brief  Context graph class factory implementation
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/scoring/ContextGraphFactory.hh>

// Package Headers
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/TwelveANeighborGraph.hh>

// STL Headers

// Utility Headers
#include <utility/exit.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {

ContextGraphOP
ContextGraphFactory::create_context_graph( ContextGraphType type ) {
	switch ( type ) {
	case ten_A_neighbor_graph :
		return ContextGraphOP( new TenANeighborGraph() );
	case twelve_A_neighbor_graph :
		return ContextGraphOP( new TwelveANeighborGraph() );
		break;
	default :
		utility_exit_with_message( "Error in ContextGraphFactory.cc.  Unsupported context graph requested" );
		break;
	}
	return nullptr;
}

} // scoring
} // core


