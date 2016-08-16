// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/loophash/LoopHashLibraryLoaderCreator.hh
/// @brief Load the Loop Hash library using the resource manager
/// @author Tim Jacobs

#ifndef INCLUDED_protocols_loophash_LoopHashLibraryLoaderCreator_hh
#define INCLUDED_protocols_loophash_LoopHashLibraryLoaderCreator_hh

//unit headers
#include <basic/resource_manager/ResourceLoaderCreator.hh>

namespace protocols {
namespace loophash {

/// @brief %LoopHashLibraryLoaderCreator allows the ResourceLoaderFactory to create a LoopHashLibraryLoader instance.
/// @details The LoopHashLibraryLoader class can be constructed from the string "LoopHashLibrary", which enables a user
/// to specify that this type of %resource is required for a particular %job in their XML input file.
class LoopHashLibraryLoaderCreator : public basic::resource_manager::ResourceLoaderCreator
{
public:
	/// @brief Return a up-casted owning pointer (ResourceLoaderOP) to the resource loader.
	virtual
	basic::resource_manager::ResourceLoaderOP
	create_resource_loader() const;

	/// @brief Return the string identifier for the associated ResourceLoader (LoopHashLibrary).
	virtual
	std::string loader_type() const;

};

} // namespace
} // namespace

#endif // include guard
