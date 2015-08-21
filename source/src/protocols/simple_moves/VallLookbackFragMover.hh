// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


/// @brief
/// @author TJ brunette


#ifndef INCLUDED_protocols_simple_moves_VallLookbackFragMover_HH
#define INCLUDED_protocols_simple_moves_VallLookbackFragMover_HH

// Unit Headers
#include <protocols/simple_moves/VallLookbackFragMover.fwd.hh>

// Package Headers
#include <protocols/simple_moves/FragmentMover.hh>

// Project Headers
#include <core/fragment/Frame.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/FrameList.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/indexed_structure_store/FragmentStore.hh>

// Utility headers

#include <utility/vector1.hh>

#include <map>

namespace protocols {
namespace simple_moves {

using namespace core;
using namespace fragment;

class VallLookbackFragMover : virtual public ClassicFragmentMover {
public:
	VallLookbackFragMover(core::fragment::FragSetCOP fragset, core::kinematics::MoveMapCOP movemap);

	~VallLookbackFragMover();

	Real lookbackMaxRmsd(pose::Pose & pose, Size resStart,Size resEnd) const;

	void markChangedResid(pose::Pose & pose,Size resStart, Size resEnd) const;

	virtual std::string get_name() const;

	virtual bool apply_frames( pose::Pose &pose, FrameList const& frames ) const;

};

} //simple_moves
} //protocols

#endif

