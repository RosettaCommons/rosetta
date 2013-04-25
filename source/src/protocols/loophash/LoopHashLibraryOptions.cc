// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loophash/LoopHashLibraryOptions.cc
/// @brief Load the Loop Hash library using the resource manager
/// @author Tim Jacobs


//unit headers
#include <protocols/loophash/LoopHashLibraryOptions.hh>
#include <protocols/loophash/LoopHashLibraryOptionsCreator.hh>

//utility headers
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <utility/exit.hh>

//C++ headers

namespace protocols {
namespace loophash {

///// LoopHashLibraryOptionsCreator /////
LoopHashLibraryOptionsCreator::LoopHashLibraryOptionsCreator() {}

LoopHashLibraryOptionsCreator::~LoopHashLibraryOptionsCreator() {}

basic::resource_manager::ResourceOptionsOP
LoopHashLibraryOptionsCreator::create_options() const {
	return new LoopHashLibraryOptions;
}

std::string
LoopHashLibraryOptionsCreator::options_type() const {
	return "LoopHashLibraryOptions";
}

///// LoopHashLibraryOptions /////
LoopHashLibraryOptions::LoopHashLibraryOptions() {}
LoopHashLibraryOptions::~LoopHashLibraryOptions() {}

utility::vector1<core::Size>
LoopHashLibraryOptions::loop_sizes() const{
	return loop_sizes_;
}

void
LoopHashLibraryOptions::loop_sizes(
	utility::vector1<core::Size> loop_sizes
){
	loop_sizes_=loop_sizes;
}

void
LoopHashLibraryOptions::parse_my_tag(
	utility::tag::TagPtr tag
)
{
	if(tag->hasOption("loop_sizes")){
		loop_sizes_.clear();
		utility::vector1<std::string> loop_sizes_strings =
			utility::string_split(tag->getOption< std::string >("loop_sizes"), ',');
		for(core::Size i=1; i<=loop_sizes_strings.size(); ++i)
		{
			loop_sizes_.push_back(utility::string2int(loop_sizes_strings[i]));
		}
	}
	else{
		utility_exit_with_message("You must provide a 'loop_sizes' option to the LoopHashLibrary tag");
	}
}

std::string
LoopHashLibraryOptions::type() const
{
	return "LoopHashLibraryOptions";
}

} // namespace loophash
} // namespace protocols
