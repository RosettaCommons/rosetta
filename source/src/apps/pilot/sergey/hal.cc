// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
///
/// @brief
/// @author Sergey Lyskov

#if defined(SERIALIZATION)

#include <devel/init.hh>

#include <basic/options/option.hh>
#include <basic/Tracer.hh>

#include <libzmq/include/zmq.hpp>

#include <thread>
#include <chrono>

static basic::Tracer TR( "HAL" );

auto const _bus_address_    = "inproc://bus";
auto const _server_address_ = "ipc:///tmp/rosetta-ui";

auto const _greetings_   = "READY PLAYER ONE";

using ContextSP = std::shared_ptr<zmq::context_t>;
using SocketUP  = std::unique_ptr<zmq::socket_t>;

//  Receive ZeroMQ message as a string
std::string receive_message(zmq::socket_t & socket)
{
	zmq::message_t message;
	socket.recv(&message);

	return std::string( static_cast<char*>( message.data() ), message.size() );
}

//  Send ZeroMQ message as string
bool send_message(zmq::socket_t & socket, std::string const & string_message) {

	zmq::message_t message( string_message.size() );
	memcpy(message.data(), string_message.data(), string_message.size());

	return socket.send(message);
}


SocketUP create_ui_socket(zmq::context_t & context)
{
	std::cout << "Connecting to UI server " << _server_address_ << "..." << std::endl;

	SocketUP ui( new zmq::socket_t(context, ZMQ_DEALER) );
	//socket.connect ("tcp://localhost:62055");
	ui->connect(_server_address_);

	return ui;
}

void hal_client(ContextSP context) // note: pass-by-value here is intentional
{
	zmq::socket_t bus(*context, ZMQ_PAIR);
	bus.connect(_bus_address_);

	SocketUP ui = create_ui_socket(*context);

	send_message(bus, _greetings_);

	std::vector<zmq_pollitem_t> items = { {*ui, 0, ZMQ_POLLIN, 0 }, {bus, 0, ZMQ_POLLIN, 0 } };

	int const max_ping_count = 5;
	int ping_count = max_ping_count;
	while ( true ) {
		zmq::poll( items, std::chrono::milliseconds(500) );

		if ( items [0].revents & ZMQ_POLLIN ) { // messages from ui server
			std::cout << "got message from ui: " << receive_message(*ui) << std::endl;
			ping_count = max_ping_count;
		}

		if ( items [1].revents & ZMQ_POLLIN ) { // messages from main Rosetta thread
		}

		if ( --ping_count == 0 ) {
			std::cout << "UI server seems to be dead, re-connecting..." << std::endl;

			ui.release();
			ui = create_ui_socket(*context);
			ping_count = max_ping_count;
		}
		send_message(*ui, "ping");
		std::cout << ping_count;
	}

	//  Do 10 requests, waiting each time for a response
	// for (int request_nbr = 0; request_nbr != 10; request_nbr++) {
	//     zmq::message_t request (5);
	//     memcpy (request.data (), "Hello", 5);
	//     std::cout << "Sending Hello " << request_nbr << "..." << std::endl;
	//     ui->send(request);
	//     //  Get the reply.
	//     zmq::message_t reply;
	//     ui->recv (&reply);
	//     std::cout << "Received World " << request_nbr << std::endl;
	// }
}

void hal(int argc, char * argv [])
{
	ContextSP context = std::make_shared<zmq::context_t>(1); // zmq::context_t context(1);

	zmq::socket_t bus(*context, ZMQ_PAIR);
	bus.bind(_bus_address_);

	std::thread client(hal_client, context);

	devel::init(argc, argv);

	auto greetings = receive_message(bus);
	if ( greetings == _greetings_ ) {
		std::cout << TR.bgRed << TR.Black << greetings << TR.Reset << std::endl;

		{
			// std::cout << "sizeof(long) = " << sizeof(long) << std::endl;
			// std::vector<char> v;
			// for(long i=0; i<1024*1024*1024*1L; ++i) v.push_back('a');

			TR << "Sleeping then restating..." << std::endl;
			std::this_thread::sleep_for(std::chrono::seconds(1*1024));
			execvp(argv[0], argv);
		}
	} else {
		utility_exit_with_message("Error: got unexpected greetings from server process, exiting...\nGreetings:" + greetings);
	}

	client.join();
}

//int main(int, char * [] )
int main(int argc, char * argv [])
{

	try {
		hal(argc, argv);
		//TR << "TTest ended. #@!  --------------------------------" << std::endl;
		return 0;

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}

#else // !defined(SERIALIZATION)

#include <utility/excn/Exceptions.hh>
#include <devel/init.hh>
#include <iostream>

int main(int argc, char * argv [])
{
	try {
		devel::init(argc, argv);

		std::cerr << "HAL app need to be build with extras=serialization! Aborting..." << std::endl;
		return 1;
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

#endif // defined(SERIALIZATION)
