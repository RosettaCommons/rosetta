// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/residue_selector/ResidueSelector.hh
/// @brief  The ResidueSelector class identifies a subset of residues from a Pose
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_pack_task_residue_selector_ResidueSelector_HH
#define INCLUDED_core_pack_task_residue_selector_ResidueSelector_HH

// Unit headers
#include <core/pack/task/residue_selector/ResidueSelector.fwd.hh>

// Package headers
#include <core/pose/Pose.fwd.hh>

// Basic headers
#include <basic/datacache/DataMap.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>

// C++ headers
#include <string>

namespace core {
namespace pack {
namespace task {
namespace residue_selector {

class ResidueSelector : public utility::pointer::ReferenceCount {
public:
	ResidueSelector();
	virtual ~ResidueSelector();

	/// @brief Make modifications to the subset vector given an input pose.
	/// The input ResidueSubset should already be of size pose.total_residue()
	/// Unless explicitly combining the subset such as AND OR, should clear the given subset.
	virtual void apply(
		core::pose::Pose const & pose,
		ResidueSubset & subset
	) const = 0;

	/// @brief Initialize any data members of this instance from an input tag
	/// and a DataMap object
	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datacache
	);

	virtual
	std::string
	get_name() const = 0;

};


} //namespace residue_selector
} //namespace task
} //namespace pack
} //namespace core


#endif
