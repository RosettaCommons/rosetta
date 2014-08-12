// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /rosetta/main/source/src/protocols/features/strand_assembly/SandwichFragment.cc
/// @brief
/// @author Doo Nam Kim based on Tim jacobs' code

//Unit Headers
#include <protocols/features/strand_assembly/SandwichFragment.hh>

namespace protocols {
namespace features {
namespace strand_assembly {

SandwichFragment::SandwichFragment()
{}

// 1 parameter
SandwichFragment::SandwichFragment(
							   core::Size chain_B_resNum):
	chain_B_resNum_(chain_B_resNum),
	pdb_source_("")
{}


// 2 parameters
SandwichFragment::SandwichFragment(
							   core::Size start,
							   core::Size end):
	start_(start),
	end_(end),
	pdb_source_("")
{}

// 3 parameters
	// used in get_current_strands_in_sheet
SandwichFragment::SandwichFragment(
							   core::Size sheet_id,
							   core::Size start,
							   core::Size end):
	sheet_id_(sheet_id),
	start_(start),
	end_(end),
	pdb_source_("")
{}

// 4 parameters
SandwichFragment::SandwichFragment(
							   core::Size sheet_id,
							   core::Size strand_id,
							   core::Size start,
							   core::Size end):
	sheet_id_(sheet_id),
	strand_id_(strand_id),
	start_(start),
	end_(end),
	pdb_source_("")
{}

// 5 parameters
SandwichFragment::SandwichFragment(
								   core::Size sw_can_by_sh_id,
								   core::Size sheet_id,
								   core::Size strand_id,
								   core::Size start,
								   core::Size end):
	sw_can_by_sh_id_(sw_can_by_sh_id),
	sheet_id_(sheet_id),
	strand_id_(strand_id),
	start_(start),
	end_(end),
	pdb_source_("")
{}

SandwichFragment::~SandwichFragment()
{}

core::Size
SandwichFragment::get_start() const
{
	return start_;
}

core::Size
SandwichFragment::get_end() const
{
	return end_;
}


core::Size
SandwichFragment::get_sw_can_by_sh_id() const
{
	return sw_can_by_sh_id_;
}

core::Size
SandwichFragment::get_sheet_id() const
{
	return sheet_id_;
}

core::Size
SandwichFragment::get_strand_id() const
{
	return strand_id_;
}


core::Size
SandwichFragment::get_resNum() const
{
	return chain_B_resNum_;
}

core::Size
SandwichFragment::get_size() const
{
	if(end_ == start_){return 0;}
	return end_-start_+1;
}


} //namespace strand_assembly
} //namespace features
} //namespace protocols
