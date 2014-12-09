// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file Node.hh
///
/// @brief
/// @author Tim Jacobs

#ifndef INCLUDED_devel_sewing_Node_HH
#define INCLUDED_devel_sewing_Node_HH

//Unit headers
#include <devel/sewing/Node.fwd.hh>

//Core headers
#include <core/types.hh>

//Utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <protocols/features/FeaturesReporter.fwd.hh>

//Boost


//C++ headers
#include <map>
#include <string>

namespace devel {
namespace sewing {

class Node : public utility::pointer::ReferenceCount
{
	
public:
	
	Node();
	
	Node(
		 protocols::features::StructureID struct_id,
		 core::Size bundle_id,
		 core::Size pair_id,
		 core::Size helix_1_id,
		 core::Size helix_1_begin,
		 core::Size helix_1_end,
		 core::Size helix_2_id,
		 core::Size helix_2_begin,
		 core::Size helix_2_end
	);
	
	protocols::features::StructureID struct_id() const;
	core::Size bundle_id() const;
	core::Size pair_id() const;
	core::Size helix_1_id() const;
	core::Size helix_1_begin() const;
	core::Size helix_1_end() const;
	core::Size helix_2_id() const;
	core::Size helix_2_begin() const;
	core::Size helix_2_end() const;
	core::Size addition_direction() const;
	std::string print() const;
	
	core::Size reversed() const;
	
	bool operator<(const Node & rhs) const;
	bool operator==(const Node & rhs) const;
	
private:
	
	protocols::features::StructureID struct_id_;
	core::Size bundle_id_;
	core::Size pair_id_;
	
	core::Size helix_1_id_;
	core::Size helix_1_begin_;
	core::Size helix_1_end_;
	
	core::Size helix_2_id_;
	core::Size helix_2_begin_;
	core::Size helix_2_end_;
	
	core::Size addition_direction_;//this doesn't work
	
public:
	std::map<core::Size, core::Size> helix_1_positions_; //map residue numbers from this node to positions in final bundle
	std::map<core::Size, core::Size> helix_2_positions_;
	
};

} //sewing namespace
} //devel namespace

#endif
