// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file basic/database/schema_generator/DatabaseSchemaGeneratorListener.hh
///
/// @brief manage the one-time initializiation of database schemas in a database
/// @author Matthew O'Meara (mattjomeara@gmail.com)


#ifndef INCLUDED_basic_message_listening_DatabaseSchemaGeneratorListener_HH
#define INCLUDED_basic_message_listening_DatabaseSchemaGeneratorListener_HH


#include <basic/message_listening/DatabaseSchemaGeneratorListener.fwd.hh>
#include <basic/message_listening/MessageListener.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <string>
#include <set>

namespace basic {
namespace message_listening {

class DatabaseSchemaGeneratorListener : public MessageListener
{
public:

	//recieve data from a slave node
	void
	receive(
		std::string const & data);

	//Given the identifier for a slave node, fill the return data to be sent
	//back and specify whether the slave should give more data. The use case for
	//this would be when you have some logic in the slave that you want executed
	//only once, or only once for each unique identifier.
	bool request(
		std::string const & table_name,
		std::string & return_data);

private:
	std::set<std::string> table_names_;

};

} //namespace
} //namespace

#endif
