// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief  metric for Gunn moves
/// @author Oliver Lange ( olange@u.washington.edu )
/// @date   Wed Aug 22 12:08:31 2007
///

#ifndef INCLUDED_protocols_frag_picker_scores_GunnCost_HH
#define INCLUDED_protocols_frag_picker_scores_GunnCost_HH


// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

#ifdef WIN32
#include <core/pose/Pose.hh> // WIN32 INCLUDE
#endif

// Utility headers
#include <utility/vector1.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

using namespace core;

struct GunnTuple {
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

class GunnCost  {
public:
	GunnCost() /*: cutoff_(7.0)*/ {}
	GunnCost( core::Real /*cutoff*/ ) /*: cutoff_( cutoff)*/ {}
	~GunnCost() {}

	void compute_gunn( core::pose::Pose const& pose, core::Size begin, core::Size end, GunnTuple &data);
	core::Real score_tuple( GunnTuple const& g1, GunnTuple const& g2 );

private:
	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// Real cutoff_;
//	utility::vector1< core::pose::PoseOP > various_length_poses_;
};


} // scores
} // frag_picker
} // protocols

#endif
