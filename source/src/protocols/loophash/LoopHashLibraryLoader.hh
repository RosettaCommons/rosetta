// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/loophash/LoopHashLibraryLoader.hh
/// @brief Load the Loop Hash library using the resource manager
/// @author Tim Jacobs

#ifndef INCLUDED_protocols_loophash_LoopHashLibraryLoader_hh
#define INCLUDED_protocols_loophash_LoopHashLibraryLoader_hh

//unit headers
#include <protocols/loophash/LoopHashLibraryLoader.fwd.hh>
#include <protocols/loophash/LoopHashLibraryOptions.hh>
#include <basic/resource_manager/ResourceLoader.hh>

//package headers
#include <basic/resource_manager/ResourceOptions.hh>
#include <basic/resource_manager/types.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>

//C++ headers
#include <istream>

namespace protocols {
namespace loophash {

/// @brief %LoopHashLibraryLoader constructs a LoopHashLibrary instance from data provided by the %ResourceManager.
/// @details The %LoopHashLibraryLoader is given a LoopHashLibraryOptions containing a %vector of loop lengths from the
/// ResourceManager.  This information is then used to produce a LoopHashLibraryOP to return to the protocol.
class LoopHashLibraryLoader : public basic::resource_manager::ResourceLoader
{
public:
	/// @brief Construct the %LoopHashLibraryLoader.
	LoopHashLibraryLoader();

	/// @brief Destructor.
	virtual ~LoopHashLibraryLoader() {}

	/// @brief Return a LoopHashLibraryOP constructed from the given ResourceOptions.
	virtual
	basic::resource_manager::ResourceOP
	create_resource(
		basic::resource_manager::ResourceOptions const & options,
		basic::resource_manager::LocatorID const &,
		std::istream &
	) const;

	/// @brief Return a ResourceOptionsOP with the default set of options.
	virtual
	basic::resource_manager::ResourceOptionsOP
	default_options() const { return basic::resource_manager::ResourceOptionsOP( new LoopHashLibraryOptions() );}
};

} // namespace
} // namespace

#endif
