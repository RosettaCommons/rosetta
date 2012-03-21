// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file util.cc
///
/// @brief

/// @author tim

#include <core/types.hh>

#include <protocols/features/util.hh>
#include <protocols/jd2/util.hh>
#include <protocols/features/ProtocolFeatures.hh>
#include <protocols/features/BatchFeatures.hh>

#include <basic/Tracer.hh>
#include <basic/database/sql_utils.hh>

#include <utility/string_util.hh>
#include <utility/exit.hh>

#include <cppdb/frontend.h>

#include <string>
#include <utility>
#include <iostream>


namespace protocols{
namespace features{

// Static data for the serial case
core::Size static_protocol_id_ = 0;
bool protocol_table_initialized_ = false;
std::map<std::string, core::Size> static_batch_id_map_;    

static basic::Tracer TR("protocols.features.util");
// End static data
    
std::pair<core::Size, core::Size> get_protocol_and_batch_id(std::string identifier, utility::sql_database::sessionOP db_session) {
    
    int protocol_id, batch_id;
    ProtocolFeaturesOP protocol_features = new ProtocolFeatures();
    BatchFeaturesOP batch_features = new BatchFeatures();
    
#ifdef USEMPI
    
    //Send an identifier to the head node, along with a message_listening tag so that the 
    //messageListenerFactory of the job distributor knows who to give the data to.
    std::string listener_data = protocols::jd2::request_data_from_head_node(protocols::jd2::message_listening::DB_TAG, identifier);
    
    TR << "recieved data from head node: " << listener_data << std::endl;
    
    //deserialize data into protocol_id and batch_id
    std::pair<core::Size, core::Size> ids = deserialize_db_listener_data(listener_data);
    protocol_id = ids.first;
    batch_id = ids.second;
    
    //no protocol id set yet - create protocol and first batch
    if(protocol_id==0){
        protocol_id = protocol_features->report_features(protocol_id, db_session);
        
        //currently batches don't have descriptions
        batch_id = batch_features->report_features(protocol_id, identifier, "", db_session);
        
        protocols::jd2::send_data_to_head_node(protocols::jd2::message_listening::DB_TAG, serialize_ids(protocol_id, identifier, batch_id));
    }
    //protocol is set, but this is a new batch
    else if(batch_id==0){
        
        batch_id = batch_features->report_features(protocol_id, identifier, "", db_session);
        
        protocols::jd2::send_data_to_head_node(protocols::jd2::message_listening::DB_TAG, serialize_ids(protocol_id, identifier, batch_id));
    }
    
#endif
#ifndef USEMPI
    //serial case
    if(!protocol_table_initialized_){
        TR << "Initializing protocol table" << std::endl;
        static_protocol_id_ = protocol_features->report_features(protocol_id, db_session);
        protocol_table_initialized_ = true;
    }
    protocol_id = static_protocol_id_;
    
    if(!static_batch_id_map_.count(identifier)){
        TR << "Initializing batch table" << std::endl;
        batch_id = batch_features->report_features(protocol_id, identifier, "", db_session);
        static_batch_id_map_[identifier] = batch_id;
    }
    else{
        batch_id=static_batch_id_map_[identifier];
    }
#endif
    
    return std::pair<core::Size, core::Size>(protocol_id, batch_id);
}
    
std::pair<core::Size, core::Size> deserialize_db_listener_data(std::string data){
    std::vector< std::string > tokens = utility::split(data);
    if(tokens.size() != 2){
        utility_exit_with_message("failed to deserialize the message from master node. Message was: " + data + " You will get this message if trying to run ReportToDb mover in MPI mode with only on processor.");
    }
    int protocol_id=utility::string2int(tokens[0]);
    int batch_id=utility::string2int(tokens[1]);
    return std::pair<core::Size, core::Size>(protocol_id, batch_id);
}

std::string serialize_ids(int protocol_id, std::string identifier, core::Size batch_id){
    return utility::to_string(protocol_id) + " " + identifier + " " + utility::to_string(batch_id);
}
    
} //namespace protocols
} //namespace features

