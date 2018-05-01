// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/network/util.hh
/// @brief: various supplemental declarations for network layer
///
/// @author Sergey Lyskov

#ifdef ZEROMQ

#pragma once

#include <libzmq/include/zmq.hpp>
#include <libzmq/include/zmq_addon.hpp>

#include <memory>
#include <vector>
#include <string>

namespace protocols {
namespace network {

// message types
auto const _m_ping_          = "ping";
auto const _m_quit_          = "quit";
auto const _m_specification_ = "specification";


// socket addresses
auto const _bus_address_    = "inproc://bus";
auto const _server_address_ = "ipc:///tmp/rosetta-ui";


using ContextSP = std::shared_ptr<zmq::context_t>;
using SocketUP  = std::unique_ptr<zmq::socket_t>;

using Bytes = std::vector<std::uint8_t>;
using BytesUP = std::unique_ptr<Bytes>;


//  Receive ZeroMQ message as a string
std::string receive_message(zmq::socket_t & socket);


//  Send ZeroMQ message as string
bool send_message(zmq::socket_t & socket, std::string const & string_message, int flags = 0);



} // namespace network
} // namespace protocols

#endif // ZEROMQ
