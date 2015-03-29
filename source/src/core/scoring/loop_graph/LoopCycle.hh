// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/scoring/loop_graph/LoopCycle.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_scoring_loop_graph_LoopCycle_HH
#define INCLUDED_core_scoring_loop_graph_LoopCycle_HH

#include <core/scoring/loop_graph/LoopCycle.fwd.hh>
#include <core/scoring/loop_graph/Loop.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

namespace core {
namespace scoring {
namespace loop_graph {

	class LoopCycle: public utility::pointer::ReferenceCount {

	public:

	//constructor
		LoopCycle();

		LoopCycle( utility::vector1< Loop > const & loops );

	//destructor
	~LoopCycle();

	public:

		Loop const & loop( Size const n ) const;

		utility::vector1< Loop > const & loops() const { return loops_; }

		Size size() const { return loops_.size(); }

		Size find_index_for_loop_landing_at_domain( Size const & takeoff_domain );

		friend
		/// @brief Test IO operator for debug and Python bindings
		std::ostream & operator << ( std::ostream & os, LoopCycle const & loop_cycle);

	private:

		utility::vector1< Loop > loops_;

	};


} //loop_graph
} //scoring
} //core

#endif
