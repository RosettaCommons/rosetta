// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/loophash/LoopHashLibraryOptionsCreator.hh
/// @brief Load the Loop Hash library using the resource manager
/// @author Tim Jacobs

#ifndef INCLUDED_protocols_loophash_LoopHashLibraryOptionsCreator_hh
#define INCLUDED_protocols_loophash_LoopHashLibraryOptionsCreator_hh

//unit headers

//package headers
#include <basic/resource_manager/ResourceOptionsCreator.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>

//C++ headers

namespace protocols {
namespace loophash {

/// @brief %LoopHashLibraryOptionsCreator allows the ResourceLoaderFactory to create a LoopHashLibraryOptions instance.
/// @details The LoopHashLibraryOptions class can be constructed from the string "LoopHashLibraryOptions", which enables
/// a user to configure a LoopHashLibrary %resource in his/her resource definitions file.
class LoopHashLibraryOptionsCreator : public basic::resource_manager::ResourceOptionsCreator
{
public:
	LoopHashLibraryOptionsCreator();
	~LoopHashLibraryOptionsCreator() override;

	/// @brief Return the string identifier for the associated ResourceOptions (LoopHashLibraryOptions).

	std::string
	options_type() const override;

	/// @brief Return a up-casted owning pointer (ResourceOptionsOP) to the resource options.

	basic::resource_manager::ResourceOptionsOP
	create_options() const override;

};

} // namespace
} // namespace

#endif //INCLUDED_protocols_loophash_LoopHashLibraryOptionsCreator_hh
