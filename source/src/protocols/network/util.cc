// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/network/util.cc
/// @brief: various supplemental declarations for network layer
///
/// @author Sergey Lyskov


#ifdef ZEROMQ

#include <protocols/network/util.hh>

namespace protocols {
namespace network {


//  Receive ZeroMQ message as a string
std::string receive_message(zmq::socket_t & socket)
{
	zmq::message_t message;
	socket.recv(&message);

	return std::string( static_cast<char*>( message.data() ), message.size() );
}

//  Send ZeroMQ message as string
bool send_message(zmq::socket_t & socket, std::string const & string_message, int flags)
{
	zmq::message_t message( string_message.size() );
	memcpy(message.data(), string_message.data(), string_message.size());

	return socket.send(message, flags);
}


} // namespace network
} // namespace protocols

#endif // ZEROMQ
