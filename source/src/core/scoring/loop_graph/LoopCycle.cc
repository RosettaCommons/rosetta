// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/loop_graph/LoopCycle.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <core/scoring/loop_graph/LoopCycle.hh>
#include <core/scoring/loop_graph/Loop.hh>

namespace core {
namespace scoring {
namespace loop_graph {

//Constructor
LoopCycle::LoopCycle()
{}

LoopCycle::LoopCycle( utility::vector1< Loop > const & loops ):
	loops_( loops )
{}

//Destructor
LoopCycle::~LoopCycle()
{}

//////////////////////////////////////////////
Loop const &
LoopCycle::loop( Size const n ) const{
	return loops_[ n ];
}

//////////////////////////////////////////////
Size
LoopCycle::find_index_for_loop_landing_at_domain( Size const & domain ) const
{
	for ( Size n = 1; n <= loops_.size(); n++ ) {
		if ( loops_[ n ].landing_domain() == domain ) return n;
	}
	return 0;
}

/// @brief Test IO operator for debug and Python bindings
std::ostream & operator << ( std::ostream & os, LoopCycle const & loop_cycle){
	for ( Size n = 1; n <= loop_cycle.loops_.size(); n++  ) {
		os << " " << loop_cycle.loops_[ n ];
	}
	return os;
}

} //loop_graph
} //scoring
} //core
