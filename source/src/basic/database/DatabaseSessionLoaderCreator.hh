// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/database/DatabaseSessionLoaderCreator.hh
/// @brief  load the database session
/// @author Tim Jacobs

#ifndef INCLUDED_basic_database_DatabaseSessionLoaderCreator_hh
#define INCLUDED_basic_database_DatabaseSessionLoaderCreator_hh

//unit headers
#include <basic/resource_manager/ResourceLoaderCreator.hh>


namespace basic {
namespace database {

class DatabaseSessionLoaderCreator : public basic::resource_manager::ResourceLoaderCreator
{
public:
	virtual
	basic::resource_manager::ResourceLoaderOP
	create_resource_loader() const;

	virtual
	std::string loader_type() const;

};

} // database
} // basic

#endif // include guard
