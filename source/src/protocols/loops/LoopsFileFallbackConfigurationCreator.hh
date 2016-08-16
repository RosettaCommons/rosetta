// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/resource_manager/LoopsFileFallbackConfigurationCreator.hh
/// @brief
/// @author Brian D. Weitzner brian.weitzner@gmail.com

#ifndef INCLUDED_protocols_loops_loops_file_fallback_configuration_creator_HH
#define INCLUDED_protocols_loops_loops_file_fallback_configuration_creator_HH

//unit headers
#include <basic/resource_manager/FallbackConfigurationCreator.hh>

// package headers
#include <basic/resource_manager/FallbackConfiguration.fwd.hh>

//C++ headers
#include <string>

namespace protocols {
namespace loops {

/// @brief %LoopsFileFallbackConfigurationCreator allows the ResourceManager to create a LoopsFileIO instance.
/// @details The LoopsFileIO class can be constructed from the string "loops_file", which provides backwards
/// compatibility with the options system.
class LoopsFileFallbackConfigurationCreator : public basic::resource_manager::FallbackConfigurationCreator
{
public:
	/// @brief Return a up-casted owning pointer (FallbackConfigurationOP) to the resource.
	virtual
	basic::resource_manager::FallbackConfigurationOP
	create_fallback_configuration() const;

	/// @brief Return the string identifier for the associated Resource (loops_file).
	virtual
	std::string resource_description() const;

};


} // namespace loops
} // namespace protocols

#endif // INCLUDED_protocols_loops_loops_file_fallback_configuration_creator_HH
