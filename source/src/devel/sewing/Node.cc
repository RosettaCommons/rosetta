// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file Node.cc
///
/// @brief
/// @author Tim Jacobs

//Unit headers
#include <devel/sewing/Node.hh>

//Core headers
#include <core/types.hh>

//Utility headers
#include <utility/string_util.hh>

//Boost


//C++ headers
#include <map>
#include <string>

namespace devel {
namespace sewing {

Node::Node():
	pair_id_(0),
	addition_direction_(0)
{}

Node::Node(
	protocols::features::StructureID struct_id,
	core::Size bundle_id,
	core::Size pair_id,
	core::Size helix_1_id,
	core::Size helix_1_begin,
	core::Size helix_1_end,
	core::Size helix_2_id,
	core::Size helix_2_begin,
	core::Size helix_2_end
):
	struct_id_(struct_id),
	bundle_id_(bundle_id), //change to substructure_id
	pair_id_(pair_id),
	helix_1_id_(helix_1_id), //change to list of substructures (substructure class has id, beginning, end)
	helix_1_begin_(helix_1_begin),
	helix_1_end_(helix_1_end),
	helix_2_id_(helix_2_id),
	helix_2_begin_(helix_2_begin),
	helix_2_end_(helix_2_end),
	addition_direction_(0)
{}

core::Size Node::addition_direction() const { return addition_direction_; }
protocols::features::StructureID Node::struct_id() const{ return struct_id_; }
core::Size Node::bundle_id() const{ return bundle_id_; }
core::Size Node::pair_id() const{ return pair_id_; }

core::Size Node::helix_1_id() const{ return helix_1_id_; }
core::Size Node::helix_1_begin() const{ return helix_1_begin_; }
core::Size Node::helix_1_end() const{ return helix_1_end_; }

core::Size Node::helix_2_id() const{ return helix_2_id_; }
core::Size Node::helix_2_begin() const{ return helix_2_begin_; }
core::Size Node::helix_2_end() const{ return helix_2_end_; }

std::string Node::print() const{
	return utility::to_string(pair_id_);
}

core::Size Node::reversed() const{ return helix_1_begin() > helix_2_begin(); }

bool Node::operator<(const Node & rhs) const{
	return (pair_id_ < rhs.pair_id());
}

bool Node::operator==(const Node & rhs) const{
	return (pair_id_ == rhs.pair_id());
}

} //sewing namespace
} //devel namespace
