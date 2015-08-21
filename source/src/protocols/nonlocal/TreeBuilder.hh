// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/TreeBuilder.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_PROTOCOLS_NONLOCAL_TREEBUILDER_HH
#define INCLUDED_PROTOCOLS_NONLOCAL_TREEBUILDER_HH

// Unit headers
#include <protocols/nonlocal/TreeBuilder.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/loops/Loops.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace nonlocal {

class TreeBuilder : public utility::pointer::ReferenceCount {
public:
	/// @brief Programmatically constructs a FoldTree, updating <pose>.
	virtual void set_up(const protocols::loops::Loops& chunks, core::pose::Pose* pose) = 0;

	/// @brief Reverts any modifications to <pose> introduced in preceding calls
	/// to set_up(). Only subclasses that introduce modifications are responsible
	/// for overriding this method.
	virtual void tear_down(core::pose::Pose*) {}
};

}  // namespace nonlocal
}  // namespace protocols

#endif  // PROTOCOLS_NONLOCAL_TREEBUILDER_HH_
