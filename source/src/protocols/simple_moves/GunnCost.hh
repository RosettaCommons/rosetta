// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief  metric for Gunn moves
/// @author Oliver Lange ( olange@u.washington.edu )
/// @date   Wed Aug 22 12:08:31 2007


#ifndef INCLUDED_protocols_simple_moves_GunnCost_HH
#define INCLUDED_protocols_simple_moves_GunnCost_HH

// Unit Headers
//#include <protocols/simple_moves/GunnTuple.fwd.hh>


// Package Headers
#include <core/fragment/Frame.fwd.hh>
#include <core/fragment/FragCache.hh>
#include <protocols/simple_moves/SmoothFragmentMover.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

#ifdef WIN32
#include <core/pose/Pose.hh> // WIN32 INCLUDE
#endif

// Utility headers
#include <utility/vector1.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <map>

#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves {


struct GunnTuple : public utility::pointer::ReferenceCount {
public:
	GunnTuple() {
		q1 = 0;
		q2 = 0;
		q3 = 0;
		q4 = 0;
		q5 = 0;
		q6 = 0;
	}
	core::Real q1, q2, q3, q4, q5, q6;
};

class GunnCost : public FragmentCost {
	friend class GunnCostTest;
public:
	GunnCost();
	GunnCost( core::Real cutoff );
	~GunnCost();

	void score( core::fragment::Frame const&, core::pose::Pose const& pose, ScoreList& scores );

	//private:
	void compute_gunn( core::fragment::Frame const& frame, core::Size frag_num, GunnTuple &data );
	void compute_gunn( core::pose::Pose const& pose, core::Size begin, core::Size end, GunnTuple &data);
	core::Real score_tuple( GunnTuple const& g1, GunnTuple const& g2 );
	core::fragment::FragCache<GunnTuple> frag_cache_;

private:
	utility::vector1< core::pose::PoseOP > various_length_poses_;

};


} //simple_moves
} //protocols
#endif
