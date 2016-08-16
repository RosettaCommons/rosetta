// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/io/silent/SilentFileLoader.hh
/// @brief  Load a silent file from an input stream and a SilentFileOptions object
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_core_io_silent_SilentLoader_hh
#define INCLUDED_core_io_silent_SilentLoader_hh


// Unit Headers
#include <core/io/silent/SilentFileLoader.fwd.hh>

// Project Headers
#include <basic/resource_manager/ResourceLoader.hh>
#include <basic/resource_manager/ResourceOptions.fwd.hh>
#include <basic/resource_manager/types.hh>

// Platform Headers
#include <core/types.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
#include <string>
#include <istream>

namespace core {
namespace io {
namespace silent {

class SilentFileLoader : public basic::resource_manager::ResourceLoader
{
public:
	SilentFileLoader();

	~SilentFileLoader();

	SilentFileLoader(
		SilentFileLoader const & src);

	basic::resource_manager::ResourceOP
	create_resource(
		basic::resource_manager::ResourceOptions const & options,
		basic::resource_manager::LocatorID const & locator_id,
		std::istream & istream) const;

	basic::resource_manager::ResourceOptionsOP
	default_options() const;

};


} // namespace
} // namespace
} // namespace


#endif
