// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/database/DatabaseSessionLoader.hh
/// @brief  load the utility database session
/// @author Tim Jacobs

#ifndef INCLUDED_basic_database_DatabaseSessionLoader_hh
#define INCLUDED_basic_database_DatabaseSessionLoader_hh

//unit headers
#include <basic/database/DatabaseSessionOptions.hh>
#include <basic/resource_manager/ResourceLoader.hh>

//package headers
#include <basic/resource_manager/ResourceOptions.hh>
#include <basic/resource_manager/types.hh>

//C++ headers
#include <istream>

namespace basic {
namespace database {

class DatabaseSessionLoader : public basic::resource_manager::ResourceLoader
{
public:
	~DatabaseSessionLoader() override = default;


	basic::resource_manager::ResourceOP
	create_resource(
		basic::resource_manager::ResourceOptions const &,
		basic::resource_manager::LocatorID const &,
		std::istream & istream
	) const override;


	basic::resource_manager::ResourceOptionsOP
	default_options() const override { return basic::resource_manager::ResourceOptionsOP( new DatabaseSessionOptions() );}

};

} // database
} // basic

#endif // include guard
