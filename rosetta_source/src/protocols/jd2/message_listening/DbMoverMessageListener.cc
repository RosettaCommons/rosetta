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

#include <protocols/jd2/message_listening/DbMoverMessageListener.hh>

#include <utility/string_util.hh>

#include <basic/Tracer.hh>

namespace protocols {
namespace jd2 {
namespace message_listening { 

static basic::Tracer TR("protocols.jd2.message_listening.DbMoverMessageListener");
    
DbMoverMessageListener::DbMoverMessageListener():
protocol_id_(0)
{}
    
bool DbMoverMessageListener::request(std::string identifier, std::string & return_data) {
    
    TR << "slave request with identifier: " << identifier << std::endl;
    
    bool need_slave_data=false;
    
    if(protocol_id_==0){need_slave_data=true;}
    
    if(!batch_ids_.count( identifier )){
        TR << "new batch id recieved" << std::endl;
        //haven't seen this mover before, initialize to 0
        batch_ids_[identifier]=0;
        need_slave_data=true;
    }
    return_data = utility::to_string(protocol_id_) + " " + utility::to_string(batch_ids_[identifier]);
    return need_slave_data;
}

void DbMoverMessageListener::recieve(std::string data) {
    
    TR << "recieved this message from slave: " << data << std::endl;
    
    deserialize_data(data);
}
    
void DbMoverMessageListener::deserialize_data(std::string data){
    std::vector< std::string > tokens = utility::split(data);
    if(tokens.size() != 3){
        utility_exit_with_message("failed to deserialize the database message from slave node.");
    }
    protocol_id_=utility::string2int(tokens[0]);    
    batch_ids_[tokens[1]]=utility::string2int(tokens[2]);
}

    
} //namespace message_listening
} //namespace jd2
} //namespace protocols
