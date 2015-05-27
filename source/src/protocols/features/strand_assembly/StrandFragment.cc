// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /rosetta/main/source/src/protocols/features/strand_assembly/StrandFragment.cc
/// @brief Small helper class that stores the start and end of a strand secondary structure
/// @author Doo Nam Kim (doonam.kim@gmail.com) based on Tim jacobs' helix assembly

//Unit Headers
#include <protocols/features/strand_assembly/StrandFragment.hh>

//Core
#include <core/conformation/Residue.hh>

//Utility
#include <utility/string_util.hh>

//Numeric
#include <numeric/xyzVector.hh>
#include <numeric/HomogeneousTransform.hh>

namespace protocols {
namespace features {
namespace strand_assembly {

StrandFragment::StrandFragment()
{}

// 2 parameters
StrandFragment::StrandFragment(core::Size start, core::Size end):
start_(start),
end_(end),
pdb_source_("")
{}

// 3 parameters	arguments
StrandFragment::StrandFragment(core::Size beta_selected_segments_id, core::Size start, core::Size end):
id_(beta_selected_segments_id),
start_(start),
end_(end),
beta_id_i_(start), //beta_selected_segments_id_i
beta_id_j_(end), //beta_selected_segments_id_j
pdb_source_("")
{}

StrandFragment::~StrandFragment(){}
	
core::Size StrandFragment::get_id() const
{
	return id_;
}

core::Size StrandFragment::get_start() const
{
	return start_;
}

core::Size StrandFragment::get_end() const
{
	return end_;
}

core::Size StrandFragment::get_i() const
{
	return beta_id_i_;
}
	
core::Size StrandFragment::get_j() const
{
	return beta_id_j_;
}
	
core::Size StrandFragment::get_size() const
{
	if(end_ == start_){return 0;}
	return end_-start_+1;
}

std::string StrandFragment::get_pdb_source() const
{
	return pdb_source_;
}

void StrandFragment::set_pdb_source(std::string pdb_source_)
{
	this->pdb_source_ = pdb_source_;
}

} //namespace strand_assembly
} //namespace features
} //namespace protocols

