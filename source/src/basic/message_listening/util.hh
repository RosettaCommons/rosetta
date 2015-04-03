// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/message_listening/util.hh
/// @author Tim Jacobs
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_basic_message_listening_util_hh
#define INCLUDED_basic_message_listening_util_hh

#include <basic/message_listening/MessageListener.fwd.hh>
#include <string>

namespace basic {
namespace message_listening {

/// @brief send mpi message to head node in order to request data. The
/// data sent back will be a string formatted by the listener
/// specified in the listener_tags enum of MessageListenerFactory
std::string
request_data_from_head_node(
	listener_tags listener_tag,
	std::string const & data);


void
send_data_to_head_node(
	std::string const & data);

std::string
listener_tag_to_name(listener_tags tag);

listener_tags
name_to_listener_tag(std::string const & listener_tag_name);


} // namespace
} // namespace

#endif
