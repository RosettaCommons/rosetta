// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/nonlocal/SimpleTreeBuilder.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_PROTOCOLS_NONLOCAL_SIMPLETREEBUILDER_HH
#define INCLUDED_PROTOCOLS_NONLOCAL_SIMPLETREEBUILDER_HH

// Unit headers
#include <protocols/nonlocal/SimpleTreeBuilder.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/loops/Loops.fwd.hh>

// Package headers
#include <protocols/nonlocal/TreeBuilder.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace nonlocal {

class SimpleTreeBuilder : public TreeBuilder {
public:
	void set_up(const protocols::loops::Loops& chunks, core::pose::Pose* pose) override;
};

}  // namespace nonlocal
}  // namespace protocols

#endif  // PROTOCOLS_NONLOCAL_SIMPLE_TREE_BUILDER_HH_
