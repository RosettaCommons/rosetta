// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/features/strand_assembly/SandwichFragment.hh
/// @brief Small helper class that stores the start and end of a strand secondary structure
/// @author Doo Nam Kim (started from Tim jacobs' code)

#ifndef INCLUDED_protocols_features_strand_assembly_SANDWICHFRAGMENT_HH
#define INCLUDED_protocols_features_strand_assembly_SANDWICHFRAGMENT_HH

//C++ Headers
#include <string>
#include <vector>
#include <map>

//Core
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>

//Numeric
#include <numeric/xyzVector.hh>
#include <numeric/HomogeneousTransform.hh>

//Utility
#include <utility/string_util.hh>
#include <utility/vector1.fwd.hh>

namespace protocols {
namespace features {
namespace strand_assembly {

class SandwichFragment{

public:
	SandwichFragment(
				   core::Size chain_B_resNum);

	SandwichFragment(
				   core::Size start,
				   core::Size end);

	SandwichFragment(
				   core::Size sheet_id,
				   core::Size start,
				   core::Size end);

	SandwichFragment(
				   core::Size sheet_id,
				   core::Size strand_id,
				   core::Size start,
				   core::Size end);

	SandwichFragment(
				   core::Size sw_can_by_sh_id,
				   core::Size sheet_id,
				   core::Size strand_id,
				   core::Size start,
				   core::Size end);

	//This really shouldn't exist, but boost uses it for serialization (I couldn't figure out the way around this).
	//Either way, start_ and end_ are read-only members, so a Fragment assembled in this way won't be very useful....
	SandwichFragment();
	~SandwichFragment();

	// Undefined, commenting out to fix PyRosetta build  std::string get_pdb_source() const;

	core::Size get_sw_can_by_sh_id() const;

	core::Size get_sheet_id() const;
	core::Size get_strand_id() const;

	core::Size get_start() const;
	core::Size get_end() const;

	core::Size get_size() const;

	core::Size get_resNum() const;

	// Undefined, commenting out to fix PyRosetta build  void set_pdb_source(std::string pdb_source_);

private:
	core::Size	sw_can_by_sh_id_;

	core::Size sheet_id_;
	core::Size strand_id_;
	core::Size start_;
	core::Size end_;

	core::Size chain_B_resNum_;

	std::string pdb_source_;
};

} //namespace strand_assembly
} //namespace features
} //namespace protocols

#endif /* SANDWICHFRAGMENT_HH_ */
