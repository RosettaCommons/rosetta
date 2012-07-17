// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/loops/LoopsFileOptions.hh
/// @brief
/// @author

//unit headers
#include <protocols/loops/LoopsFileOptions.hh>
#include <protocols/loops/LoopsFileOptionsCreator.hh>

//utility headers
#include <utility/tag/Tag.hh>

//C++ headers

namespace protocols {
namespace loops {

std::string
LoopsFileOptionsCreator::options_type() const { return "LoopsFileOptions"; }

basic::resource_manager::ResourceOptionsOP
LoopsFileOptionsCreator::create_options() const { return new LoopsFileOptions; }


LoopsFileOptions::LoopsFileOptions() : prohibit_single_residue_loops_( true ) {}
LoopsFileOptions::~LoopsFileOptions() {}

void
LoopsFileOptions::parse_my_tag(
	utility::tag::TagPtr tag
)
{
	prohibit_single_residue_loops( tag->getOption< bool >( "prohibit_single_residue_loops", true ));
}

std::string
LoopsFileOptions::type() const
{
	return "LoopsFileOptions";
}

bool LoopsFileOptions::prohibit_single_residue_loops() const { return prohibit_single_residue_loops_; }
void LoopsFileOptions::prohibit_single_residue_loops( bool setting ) { prohibit_single_residue_loops_ = setting; }

} // namespace loops
} // namespace protocols
