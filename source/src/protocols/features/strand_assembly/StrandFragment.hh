// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StrandFragment.hh
/// @brief Small helper class that stores the start and end of a strand secondary structure
/// @author Doo Nam Kim (based on Tim jacobs' code)

#ifndef INCLUDED_protocols_features_strand_assembly_STRANDFRAGMENT_HH
#define INCLUDED_protocols_features_strand_assembly_STRANDFRAGMENT_HH

//Core
#include <core/types.hh>
#include <core/pose/Pose.hh>

//Utility
#include <utility/vector1.fwd.hh>

//C++ Headers
#include <string>
#include <vector>
#include <map>

namespace protocols {
namespace features {
namespace strand_assembly {

class StrandFragment{

public:

	StrandFragment(core::Size start, core::Size end);
	StrandFragment(core::Size beta_selected_segments_id, core::Size start, core::Size end);
	
	StrandFragment();
	~StrandFragment();
	
	std::string get_pdb_source() const;
	core::Size get_id() const;
	
	core::Size get_start() const;
	core::Size get_end() const;
	
	core::Size get_i() const;
	core::Size get_j() const;
	
	core::Size get_size() const;
	void set_pdb_source(std::string pdb_source_);

private:

	core::Size id_;
	core::Size start_;
	core::Size end_;
	
	core::Size beta_id_i_;
	core::Size beta_id_j_;
	
	std::string pdb_source_;
};

} //namespace strand_assembly
} //namespace features
} //namespace protocols

#endif /* STRANDFRAGMENT_HH_ */
