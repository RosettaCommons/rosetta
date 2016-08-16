// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/loops/LoopsFileOptionsCreator.hh
/// @brief
/// @author

#ifndef INCLUDED_protocols_loops_LoopsFileOptionsCreator_hh
#define INCLUDED_protocols_loops_LoopsFileOptionsCreator_hh

//unit headers

//package headers
#include <basic/resource_manager/ResourceOptionsCreator.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>

//C++ headers

namespace protocols {
namespace loops {

/// @brief %LoopsFileOptionsCreator allows the ResourceLoaderFactory to create a LoopsFileOptions instance.
/// @details The LoopsFileOptions class can be constructed from the string "LoopsFileOptions", which enables a user to
/// configure a LoopsFile %resource in his/her resource definitions file.
class LoopsFileOptionsCreator : public basic::resource_manager::ResourceOptionsCreator
{
public:

	/// @brief Return the string identifier for the associated ResourceOptions (LoopsFileOptions).
	virtual std::string options_type() const;

	/// @brief Return a up-casted owning pointer (ResourceOptionsOP) to the resource options.
	virtual basic::resource_manager::ResourceOptionsOP create_options() const;

};

} // namespace loops
} // namespace protocols

#endif //INCLUDED_protocols_loops_LoopsFileOptionsCreator_hh
