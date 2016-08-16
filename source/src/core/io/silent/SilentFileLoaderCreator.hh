// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


/// @file   core/io/silent/SilentFileLoaderCreator.fwd.hh
/// @brief  Creator for loader for constructing a pose from a silent file
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_core_io_silent_SilentFileLoaderCreator_hh
#define INCLUDED_core_io_silent_SilentFileLoaderCreator_hh

// Unit Headers
#include <basic/resource_manager/ResourceLoaderCreator.hh>
#include <core/io/silent/SilentFileLoader.fwd.hh>

#include <core/types.hh>
#include <utility/vector1.hh>
#include <string>

namespace core {
namespace io {
namespace silent {

/// @brief creator for the SilentFileLoader class
class SilentFileLoaderCreator : public basic::resource_manager::ResourceLoaderCreator
{
public:
	SilentFileLoaderCreator();
	virtual ~SilentFileLoaderCreator();

	virtual basic::resource_manager::ResourceLoaderOP create_resource_loader() const;
	virtual std::string loader_type() const;
};

} //namespace
} //namespace
} //namespace

#endif
