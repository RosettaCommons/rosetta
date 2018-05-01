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

#include <utility/json_utilities.hh>


#if defined(ZEROMQ)  and  defined(_NLOHMANN_JSON_ENABLED_)

#include <devel/init.hh>

#include <basic/options/option.hh>
#include <basic/Tracer.hh>

#include <utility/CSI_Sequence.hh>

#include <protocols/network/util.hh>

#include <json.hpp>

#include <thread>
#include <chrono>

using std::string;
using namespace utility;
using namespace protocols::network;

static basic::Tracer TR( "HAL" );

auto const _greetings_   = "READY PLAYER ONE";


SocketUP create_ui_socket(zmq::context_t & context)
{
	//std::cout << "Connecting to UI server " << _server_address_ << "..." << std::endl;

	SocketUP ui( new zmq::socket_t(context, ZMQ_DEALER) );
	//socket.connect ("tcp://localhost:62055");
	ui->connect(_server_address_);

	return ui;
}


// receive specification or terminate if thats not possible
std::string receive_specification(zmq::socket_t & bus)
{
	zmq::multipart_t message(bus);

	if( message.size() == 2 ) {
		string type = message.popstr();
		if( type == _m_specification_ ) {
			return message.popstr();
		}
		else std::cerr << "ERROR: first message from Rosetta should be specification! Instead received: " << type << std::endl;

	}
	std::cerr << "ERROR: first message from Rosetta should be specification! Terminating..." << std::endl;
	std::exit(1);

	return ""; // dummy, should never be reached
}

void hal_client(ContextSP context) // note: pass-by-value here is intentional
{
	zmq::socket_t bus(*context, ZMQ_PAIR);
	bus.connect(_bus_address_);

	SocketUP ui = create_ui_socket(*context);

	send_message(bus, _greetings_);

	std::string specification = receive_specification(bus);
	std::cout << "Got specification from Rosetta: " << nlohmann::json::from_msgpack(specification) << std::endl;

	auto create_pool = [&ui, &bus]() -> std::vector<zmq_pollitem_t> { return { {*ui, 0, ZMQ_POLLIN, 0 }, {bus, 0, ZMQ_POLLIN, 0 } }; };

	std::vector<zmq_pollitem_t> items = create_pool();  //{ {*ui, 0, ZMQ_POLLIN, 0 }, {bus, 0, ZMQ_POLLIN, 0 } };

	int const max_ping_count = 10; // maximum number of pings we sent before re-creating a socket
	int ping_count = max_ping_count;

	auto const timeout = std::chrono::milliseconds(2500);
	auto const ping_interval = std::chrono::milliseconds(500);
	std::chrono::steady_clock::time_point last_ping_time = std::chrono::steady_clock::now() - ping_interval;

	bool alive = false;
	std::chrono::steady_clock::time_point time_point_of_ui_last_message = std::chrono::steady_clock::now();

	while ( true ) {
		if( zmq::poll( items, std::chrono::milliseconds(250) ) ) {
			if ( items [0].revents & ZMQ_POLLIN ) { // messages from ui server
				auto message = receive_message(*ui);
				//std::cout << "got message from ui: " << message << std::endl;
				time_point_of_ui_last_message = std::chrono::steady_clock::now();

				if( not alive ) {
					std::cout << CSI_bgGreen() << CSI_Black() << "Connection to UI front-end established!" << CSI_Reset() << std::endl;
					alive = true;

					// sending our specification to UI
					send_message(*ui, _m_specification_, ZMQ_SNDMORE);
					send_message(*ui, specification);
				}
				ping_count = max_ping_count;
			}

			if ( items [1].revents & ZMQ_POLLIN ) { // messages from main Rosetta thread
				zmq::multipart_t message(bus);

				if( message.size() == 2 ) {
					string type = message.popstr();
					std::cout << "got message from Rosetta:" << type << std::endl;
				}

			}
		}

		std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
		//std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(now - time_point_of_ui_last_message).count()
		//		  << " _:" << (now - time_point_of_ui_last_message > timeout) <<  std::endl;

		if( alive  and  now - time_point_of_ui_last_message > timeout ) {
			std::cout << CSI_bgRed() << CSI_Black() << "UI front-end disconnected due to timeout!" << CSI_Reset() << std::endl;
			alive = false;
			ping_count = 1; // force re-creation of socket
		}

		if ( --ping_count == 0 ) {
			ui->setsockopt(ZMQ_LINGER, 0);
			ui->close(); // ui.release();
			ui = create_ui_socket(*context);

			items = create_pool();

			ping_count = max_ping_count;
		}

		if( now - last_ping_time > ping_interval ) {
			//std::cout << "sending ping..." << std::endl;
			send_message(*ui, _m_ping_);
			last_ping_time = now;
		}
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


/// Generate HAL specification
string specification()
{
	nlohmann::json j;
	j["functions"] = { {"name", "test"}, {"arguments", {"int", "float"} } };

	string r;
	nlohmann::json::basic_json::to_msgpack(j, r);
	return r;
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

		{ // sending HAL specification
			send_message(bus, _m_specification_, ZMQ_SNDMORE);
			send_message(bus, specification() );
		}

		{
			// std::cout << "sizeof(long) = " << sizeof(long) << std::endl;
			// std::vector<char> v;
			// for(long i=0; i<1024*1024*1024*1L; ++i) v.push_back('a');

			TR << "Sleeping then restating..." << std::endl;
			std::this_thread::sleep_for(std::chrono::seconds(1*1024));
			//std::this_thread::sleep_for(std::chrono::seconds(2));
			//execvp(argv[0], argv);
		}
	} else {
		utility_exit_with_message("Error: got unexpected greetings from server process, exiting...\nGreetings:" + greetings);
	}

	client.join();
}


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

#else // ! ( defined(ZEROMQ) and  defined(_NLOHMANN_JSON_ENABLED_) )

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

#endif // defined(ZEROMQ) and  defined(_NLOHMANN_JSON_ENABLED_)
