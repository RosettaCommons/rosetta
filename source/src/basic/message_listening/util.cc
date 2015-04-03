// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file basic/message_listening/util.cc
/// @brief Utility functions for the message listening framework
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifdef USEMPI
#include <mpi.h>
#endif

#include <basic/message_listening/MessageListener.fwd.hh>
#include <basic/message_listening/util.hh>
#include <basic/Tracer.hh>

#include <utility/mpi_util.hh>
#include <utility/assert.hh>
#include <utility/exit.hh>

#include <boost/lexical_cast.hpp>

#include <string>

namespace basic {
namespace message_listening {

using std::string;
using std::endl;
using utility::send_string_to_node;
using utility::receive_string_from_node;


static thread_local basic::Tracer TR( "basic.message_listening" );

/// @brief used for message passing to the
///MPIWorkPoolJobDistributor. This function will ask the head node for
///data.  The type of data returned is based on the type of listener
///created based on the listener_tags of the MessageListenerFactory

string
request_data_from_head_node(
  listener_tags MPI_ONLY( listener_tag ) ,
  string const & MPI_ONLY( data )
){

#ifdef USEMPI

  //send a message to the head node that tells jd2 to create a message listener
  TR.Debug
    << "Requesting data from head node for tag "
    << "'" << listener_tag_to_name(listener_tag) << "'" << endl;
  TR.flush();
  MPI_Send( &listener_tag, 1, MPI_INT, 0/*head node*/, 40 /*REQUEST_MESSAGE_TAG*/, MPI_COMM_WORLD );

  //send a string to be processed by the listener
  TR.Debug << "Sending " << data << " to head node for tag '" << listener_tag_to_name(listener_tag) << "'" << std::endl;
  TR.flush();
  send_string_to_node(0/*head node*/, data);

  //receive a response from the head node listener
  return utility::receive_string_from_node(0/*head node*/);
#endif
#ifndef USEMPI
  utility_exit_with_message(
    "ERROR: You have tried to request a message from the head node but you are not in mpi mode (compile with extras=mpi)");
#endif

  return "";  // required for compilation on Windows
}

void
send_data_to_head_node(
  std::string const & MPI_ONLY( data )
){
#ifdef USEMPI
  //send a string to be processed by the listener
  utility::send_string_to_node(0/*head node*/, data);
#endif
#ifndef USEMPI
  utility_exit_with_message("ERROR: You have tried to send a message to the head node but you are not in mpi mode (compile with extras=mpi)");
#endif
}

std::string
listener_tag_to_name(listener_tags tag)
{
	switch(tag){
		case DATABASE_PROTOCOL_AND_BATCH_ID_TAG:
			return "DATABAASE_PROTOCOL_AND_BATCH_ID_TAG";
		default:
			return "Unrecognized listener_tag: " + boost::lexical_cast<std::string>(tag);
	}
}


listener_tags
name_to_listener_tag(std::string const & listener_tag_name){
  if(listener_tag_name == "DATABASE_PROTOCOL_AND_BATCH_ID_TAG"){
    return DATABASE_PROTOCOL_AND_BATCH_ID_TAG;
  } else {
    utility_exit_with_message("Unknown listener tag name '" + listener_tag_name + "'");
  }
}


} // namespace
} // namespace
