// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/scoring/loop_graph/LoopScoreInfo.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_scoring_loop_graph_LoopScoreInfo_HH
#define INCLUDED_core_scoring_loop_graph_LoopScoreInfo_HH

#include <utility/pointer/ReferenceCount.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/loop_graph/LoopScoreInfo.fwd.hh>
#include <core/scoring/func/Func.fwd.hh>
#include <core/id/AtomID.hh>

using core::scoring::constraints::FuncOP;
using namespace core;

namespace core {
namespace scoring {
namespace loop_graph {

	class LoopScoreInfo: public utility::pointer::ReferenceCount {

	public:

	//constructor
		LoopScoreInfo();

	//destructor
	~LoopScoreInfo();

	public:

		void set_length( core::Size const & setting ){ length_ = setting; }
		core::Size length() const{ return length_; }

		void set_takeoff_atom( id::AtomID const & setting ){ takeoff_atom_ = setting; }
		id::AtomID takeoff_atom() const{ return takeoff_atom_; }

		void set_landing_atom( id::AtomID const & setting ){ landing_atom_ = setting; }
		id::AtomID landing_atom() const{ return landing_atom_; }

		void set_current_distance( core::Real const & setting ){ current_distance_ = setting; }
		core::Real current_distance() const{ return current_distance_; }

		void set_func( FuncOP const & setting );
		FuncOP func() const;

	private:

		core::Size length_;
		id::AtomID takeoff_atom_;
		id::AtomID landing_atom_;
		core::Real current_distance_;

		FuncOP func_;

	};

} //loop_graph
} //scoring
} //core

#endif
