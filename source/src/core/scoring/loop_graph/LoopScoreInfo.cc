// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/scoring/loop_graph/LoopScoreInfo.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <core/scoring/loop_graph/LoopScoreInfo.hh>
#include <core/scoring/func/Func.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/pose/Pose.fwd.hh>

using core::scoring::func::FuncOP;
using namespace core;

namespace core {
namespace scoring {
namespace loop_graph {

	//Constructor
	LoopScoreInfo::LoopScoreInfo():
		length_( 0 ),
		takeoff_atom_( id::BOGUS_ATOM_ID ),
		landing_atom_( id::BOGUS_ATOM_ID ),
		current_distance_( 0.0 )
	{}

	//Destructor
	LoopScoreInfo::~LoopScoreInfo()
	{}

	//////////////////////////
	void
	LoopScoreInfo::set_func( FuncOP const & setting ){
		func_ = setting;
	}

	//////////////////////////
	FuncOP
	LoopScoreInfo::func() const{
		return func_;
	}

} //loop_graph
} //scoring
} //core
