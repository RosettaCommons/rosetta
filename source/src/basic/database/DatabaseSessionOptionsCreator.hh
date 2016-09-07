// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/database/DatabaseSessionOptionsCreator.hh
/// @brief  load the database session
/// @author Tim Jacobs

#ifndef INCLUDED_basic_database_DatabaseSessionOptionsCreator_hh
#define INCLUDED_basic_database_DatabaseSessionOptionsCreator_hh

//unit headers
#include <basic/resource_manager/ResourceOptionsCreator.hh>
#include <basic/database/DatabaseSessionOptions.fwd.hh>

namespace basic {
namespace database {

class DatabaseSessionOptionsCreator : public basic::resource_manager::ResourceOptionsCreator
{
public:
	DatabaseSessionOptionsCreator();
	~DatabaseSessionOptionsCreator() override;


	basic::resource_manager::ResourceOptionsOP
	create_options() const override;


	std::string options_type() const override;

};

} // database
} // basic

#endif // include guard
