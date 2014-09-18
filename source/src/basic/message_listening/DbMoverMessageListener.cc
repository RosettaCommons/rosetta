// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file MessageListener.cc
///
/// @brief
/// @author Tim Jacobs

#include <basic/message_listening/DbMoverMessageListener.hh>

#include <utility/string_util.hh>

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>


namespace basic {
namespace message_listening {

static thread_local basic::Tracer TR( "basic.message_listening.DbMoverMessageListener" );

DbMoverMessageListener::DbMoverMessageListener():
	protocol_id_(0)
{}

bool
DbMoverMessageListener::request(
	std::string const & identifier,
	std::string & return_data
){

	bool need_slave_data=false;

	//Has a protocol id been set yet? If not, return 0 and request a protocol id from the slave
	if(protocol_id_==0){
		need_slave_data=true;
	}
	if(!batch_ids_.count(identifier)){
		batch_ids_[identifier]=0;
		need_slave_data=true;
	}
	return_data = utility::to_string(protocol_id_) + " " + utility::to_string(batch_ids_[identifier]);
	return need_slave_data;
}

void DbMoverMessageListener::receive(std::string const & data) {

	TR << "Received this message from slave: " << data << std::endl;

	deserialize_data(data);
}

void
DbMoverMessageListener::deserialize_data(
	std::string const & data){
	utility::vector1< std::string > tokens = utility::split(data);
	if(tokens.size() != 3){
		utility_exit_with_message("failed to deserialize the database message from slave node.");
	}
	//Received a protocol id from the slave node, and thus it has been written to the database
	protocol_id_=utility::string2int(tokens[1]);

	TR << "Setting batch id identifier to: '" << tokens[2] << "'" << std::endl;
	batch_ids_[tokens[2]]=utility::string2int(tokens[3]);
}


} //namespace
} //namespace
