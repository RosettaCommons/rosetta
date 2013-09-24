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

/// @details Return a string specifying the type of %ResourceOptions to create (LoopsFileOptions).
std::string
LoopsFileOptionsCreator::options_type() const { return "LoopsFileOptions"; }

/// @details Return an owning pointer to a newly constructed default instance of LoopsFileOptions.
basic::resource_manager::ResourceOptionsOP
LoopsFileOptionsCreator::create_options() const { return new LoopsFileOptions; }

/// @details The prohibit_single_residue_loops property is set to true by default.
LoopsFileOptions::LoopsFileOptions() : prohibit_single_residue_loops_( true ) {}
LoopsFileOptions::~LoopsFileOptions() {}

/// @details Read the resource definitions file's tag and set the value for prohibit_single_residue_loops appropriately.
/// If this option is omitted it is set to true by default.
void
LoopsFileOptions::parse_my_tag(
	utility::tag::TagPtr tag
)
{
	prohibit_single_residue_loops( tag->getOption< bool >( "prohibit_single_residue_loops", true ));
}

/// @details Return the string value for the name of this class (LoopsFileOptions).
std::string
LoopsFileOptions::type() const
{
	return "LoopsFileOptions";
}

// No details necessary - the @brief's for these methods (in the header) are sufficient.
bool LoopsFileOptions::prohibit_single_residue_loops() const { return prohibit_single_residue_loops_; }
void LoopsFileOptions::prohibit_single_residue_loops( bool setting ) { prohibit_single_residue_loops_ = setting; }

} // namespace loops
} // namespace protocols
