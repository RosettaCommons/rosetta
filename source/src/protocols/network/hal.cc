// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/network/hal.cc
/// @brief: HAL implementaion
///
/// @author Sergey Lyskov

#ifdef ZEROMQ

#include <protocols/network/hal.hh>
//#include <protocols/network/hal.util.hh>
#include <protocols/network/util.hh>

#include <basic/Tracer.hh>

#include <utility/CSI_Sequence.hh>


#include <thread>
#include <chrono>
#include <mutex>
#include <condition_variable>

#include <unistd.h>

#include <cstring>


namespace protocols {
namespace network {

using std::string;
using namespace utility;
using namespace protocols::network;

static basic::Tracer TR( "HAL" );


auto const _greetings_ = "READY PLAYER ONE";

struct HAL_pause_state {
	static std::mutex hal_pause_mutex;
	static std::condition_variable condition_variable;
	static bool paused;
};

std::mutex HAL_pause_state::hal_pause_mutex;
std::condition_variable HAL_pause_state::condition_variable;
bool HAL_pause_state::paused = false;

void sleep_if_paused()
{
	std::unique_lock<std::mutex> lock(HAL_pause_state::hal_pause_mutex);
	HAL_pause_state::condition_variable.wait(lock, []{ return not HAL_pause_state::paused;} );
}

void pause_rosetta_threads(bool state)
{
	{
		std::lock_guard<std::mutex> _(HAL_pause_state::hal_pause_mutex);
		HAL_pause_state::paused = state;
	}
	HAL_pause_state::condition_variable.notify_all();
}


///
/// actual implementation of hal_client
/// Think of this as not of data type but as way to implement single, a few-pages long function
///
class H
{

public:
	H(zmq::context_t &context, CommandLineArguments const &args);
	void hal_client();



private:

	SocketUP create_ui_socket();

	void receive_specification();

	void process_ui_messages();
	void process_bus_messages();
	void process_progress_messages();


private:
	int const _ui_socket_high_water_mark_  = 16;
	int const _bus_socket_high_water_mark_ = 16;

	int const max_ping_count = 10; // maximum number of pings we sent before re-creating a socket

	zmq::context_t &context;
	CommandLineArguments const &command_line_arguments;

	zmq::socket_t bus;
	zmq::socket_t progress;

	SocketUP ui;

	int ping_count = max_ping_count;
	std::chrono::milliseconds const timeout       = std::chrono::milliseconds(2500);
	std::chrono::milliseconds const ping_interval = std::chrono::milliseconds(500);

	std::chrono::steady_clock::time_point last_ping_time;
	std::chrono::steady_clock::time_point time_point_of_ui_last_message;

	bool ui_is_alive = false;

	std::string specification;

	HAL_Settings settings;

	std::vector<zmq_pollitem_t> poll;

	// WARNING when chaning make sure that `create_poll` (below) got adjusted!
	enum SocketIndexes {
						_SI_ui_       = 0,
						_SI_bus_      = 1,
						_SI_progress_ = 2,
	};
	std::vector<zmq_pollitem_t> create_poll() {
		//       _SI_ui_                   _SI_bus_                  _SI_progress_
		return { {*ui, 0, ZMQ_POLLIN, 0 }, {bus, 0, ZMQ_POLLIN, 0 }, { progress, 0, ZMQ_POLLIN, 0 } };
	}
};



H::H(zmq::context_t &context, CommandLineArguments const &args) :
	context(context),
	command_line_arguments(args),

	bus     (context, ZMQ_PAIR),
	progress(context, ZMQ_ROUTER),

	ui( create_ui_socket() ),

	poll( create_poll() )
{
	bus.connect(_bus_address_);

	bus.setsockopt(ZMQ_SNDHWM, _bus_socket_high_water_mark_);
	bus.setsockopt(ZMQ_RCVHWM, _bus_socket_high_water_mark_);


	progress.bind(_hal_address_);


	last_ping_time = std::chrono::steady_clock::now() - ping_interval;

	time_point_of_ui_last_message = std::chrono::steady_clock::now();
}

SocketUP H::create_ui_socket()
{
	//std::cout << "Connecting to UI server " << _server_address_ << "..." << std::endl;

	SocketUP ui( new zmq::socket_t(context, ZMQ_DEALER) );
	//socket.connect ("tcp://localhost:62055");
	ui->setsockopt(ZMQ_LINGER, 0);
	ui->connect(_server_address_);
	ui->setsockopt(ZMQ_SNDHWM, _ui_socket_high_water_mark_);
	ui->setsockopt(ZMQ_RCVHWM, _ui_socket_high_water_mark_);
	return ui;
}

// receive specification or terminate if thats not possible
void H::receive_specification()
{
	zmq::multipart_t message(bus);

	if( message.size() == 2 ) {
		string type = message.popstr();
		if( type == _m_specification_ ) {
			specification = message.popstr();
			return;
		}
		else std::cerr << "ERROR: first message from Rosetta should be specification! Instead received: " << type << std::endl;

	}
	std::cerr << "ERROR: first message from Rosetta should be specification! Terminating..." << std::endl;
	std::exit(1);
}


void H::process_ui_messages()
{
	//auto message = receive_message(*ui);
	zmq::multipart_t message(*ui);
	if( message.size() >= 1 ) {
		string type = message.popstr();

		//if( type != _m_ping_ ) std::cout << "got message from ui:" << type << " extra-parts:" << message.size() << std::endl;

		if( type == _m_abort_ ) {
			std::cout << CSI_bgRed() << CSI_Black() << "got " << type << " from UI, restarting..." << CSI_Reset() << std::endl;
			execvp(command_line_arguments.argv[0], command_line_arguments.argv);
		}
		else if( type == _m_execute_  and  message.size() == 1  ) {
			//std::cout << "___ _m_execute_ "  << std::endl;
			string command = message.popstr();

			send_message(bus, _m_execute_, ZMQ_SNDMORE);
			send_message(bus, command);

		}
		else if( type == _m_settings_  and  message.size() == 1 ) {
			zmq::message_t m = message.pop();
			if( m.size() == sizeof(HAL_Settings) ) {
				HAL_Settings o;
				memcpy( &o, m.data(), sizeof(HAL_Settings));

				if( o != settings ) {
					poll = create_poll();

					settings = o;
					//std::cout << CSI_bgBlue() << CSI_Black() << "HAL: got message `" << _m_options_ << "`, updating my settings..." << CSI_Reset() << std::endl;
					if( settings.pause ) std::cout << CSI_bgBlue() << CSI_Black() << "paused" << CSI_Reset() << std::endl;
					else std::cout << CSI_bgBlue() << CSI_Black() << "running" << CSI_Reset() << std::endl;

					pause_rosetta_threads(settings.pause);
				}

			}
			else std::cout << CSI_bgRed() << CSI_Black() << "HAL: message `" << _m_settings_ << "` have wrong payload length - ignoring..." << CSI_Reset() << std::endl;
		}
		else if( type != _m_ping_ ) {
			std::cout << "HAL: Got '" << type << "' size=" << message.size() << " message from UI, type is unknown - discarding..." <<  std::endl;
		}

	}

	{ // keep-alive machinery
		time_point_of_ui_last_message = std::chrono::steady_clock::now();

		if( not ui_is_alive ) {
			std::cout << CSI_bgGreen() << CSI_Black() << "Connection to UI front-end established!" << CSI_Reset() << std::endl;
			ui_is_alive = true;

			// sending our specification to UI
			send_message(*ui, _m_specification_, ZMQ_SNDMORE);
			send_message(*ui, specification);
		}
		ping_count = max_ping_count;
	}
}


void H::process_bus_messages()
{
	zmq::multipart_t message(bus);

	if( message.size() >= 1 ) {
		string type = message.popstr();
		//std::cout << "got message from Rosetta:" << type << " size:" << message.size() << std::endl;

		if( type == _m_result_  and  message.size() == 1  ) {
			string result = message.popstr();

			send_message(*ui, _m_result_, ZMQ_SNDMORE);
			send_message(*ui, result);
		}
	}
}


void H::process_progress_messages()
{
	//std::cout << "got message from Rosetta: UIMover..." << std::endl;
	zmq::multipart_t message(progress);
	/* string mover_id = */ message.popstr();
	message.send(*ui);

	// while( message.size() ) {
	// 	string m = message.popstr();
	// 	send_message(*ui, m, message.size() ? ZMQ_SNDMORE : 0);
	// }

	//send_message(*ui, _m_intermediate_, ZMQ_SNDMORE);
	//send_message(*ui, "m_intermediate_result");
}


void H::hal_client()
{
	send_message(bus, _greetings_);

	receive_specification();
	//std::cout << "Got specification from Rosetta: " << nlohmann::json::from_msgpack(specification) << std::endl;
	std::cout << "Got specification from Rosetta, ready to connect to UI client!" << std::endl;

	while ( true ) {
		if( zmq::poll( poll, std::chrono::milliseconds(250) ) ) {

			if ( poll[_SI_ui_].revents & ZMQ_POLLIN ) process_ui_messages(); // messages from ui server


			if ( poll[_SI_bus_].revents & ZMQ_POLLIN ) process_bus_messages();  // messages from main Rosetta thread, main()

			if ( poll[_SI_progress_].revents & ZMQ_POLLIN ) process_progress_messages(); // messages from main Rosetta thread, UIMover
		}

		std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
		//std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(now - time_point_of_ui_last_message).count()
		//		  << " _:" << (now - time_point_of_ui_last_message > timeout) <<  std::endl;

		if( ui_is_alive  and  now - time_point_of_ui_last_message > timeout ) {
			std::cout << CSI_bgRed() << CSI_Black() << "UI front-end disconnected due to timeout!" << CSI_Reset() << std::endl;
			ui_is_alive = false;
			ping_count = 1; // force re-creation of socket
		}

		if( now - last_ping_time > ping_interval ) {
			if ( --ping_count == 0 ) {
				//std::cout << "__________1" << std::endl;
				//ui->close(); // ui.release();
				ui = create_ui_socket();
				poll = create_poll();
				ping_count = max_ping_count;
				//std::cout << "__________2" << std::endl;
			}

			//std::cout << "sending ping..." << std::endl;
			send_message(*ui, _m_ping_);
			last_ping_time = now;
		}

		// {
		// 	std::cout << "testing..." << std::endl;
		// 	for(int i=0; ; ++i) {
		// 		//ui->setsockopt(ZMQ_LINGER, 0);
		// 		//ui->close();
		// 		//ui.release();
		// 		ui = create_ui_socket(*context);
		// 		if( not (i % 64) ) {
		// 			std::cout << i << std::endl;
		// 			std::this_thread::sleep_for( std::chrono::milliseconds(10) );
		// 			//std::this_thread::sleep_for( std::chrono::seconds(1) );
		// 		}
		// 	}
		// }
	}


}


void hal_client(zmq::context_t &context, CommandLineArguments args) // hal_client(ContextSP context) // note: pass-by-value here is intentional
{
	H h(context, args);
	h.hal_client();

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

/// Extract human-redable executable name from given path and extension: `/my/path/hal-demo.gccrelease.default` should yield `hal-demo`
std::string extract_executable_name(CommandLineArguments const &args)
{
	if (args.argc) {
		char const *s = args.argv[0];  //std::cout << s << std::endl;

		int begin = 0;
		int end = strlen(s);

		for(int i=end; i >= 0; --i ) {
			if( s[i] == '.' ) end = i;
			if( s[i] == '/' ) { begin = i + 1; break; }
		}

		if( end > begin ) return std::string(s + begin, end - begin);
	}
	return "unknown";
}

void hal(SpecificationCallBack const &specification_callback, ExecutionerCallBack const & executioner_callback, CommandLineArguments const &args)
{
	//ContextSP context = std::make_shared<zmq::context_t>(1); // zmq::context_t context(1);
	zmq::context_t & context = zmq_context();

	zmq::socket_t bus(context, ZMQ_PAIR);
	bus.bind(_bus_address_);

	std::thread client(hal_client, std::ref(context), args);

	auto greetings = receive_message(bus);
	if ( greetings == _greetings_ ) {
		std::cout << TR.bgRed << TR.Black << greetings << TR.Reset << std::endl;

		{ // sending HAL specification
			send_message(bus, _m_specification_, ZMQ_SNDMORE);

			auto s = specification_callback();

			if( not s.count(_f_name_) ) s[_f_name_] = extract_executable_name(args);

			string r;  nlohmann::json::basic_json::to_msgpack(s, r);

			send_message(bus, r);
		}

		{
			// std::cout << "sizeof(long) = " << sizeof(long) << std::endl;
			// std::vector<char> v;
			// for(long i=0; i<1024*1024*1024*1L; ++i) v.push_back('a');

			// TR << "Sleeping then restating..." << std::endl;
			// std::this_thread::sleep_for(std::chrono::seconds(1*1024));
			// //std::this_thread::sleep_for(std::chrono::seconds(2));
			// //execvp(argv[0], argv);

			for(;;) {
				zmq::multipart_t message(bus);

				if( message.size() >= 1 ) {
					string type = message.popstr();

					//std::cout << "Rosetta: got message from ui:" << type << " extra-parts:" << message.size() << std::endl;

					if( type == _m_execute_  and  message.size() == 1  ) {
						string binary_command = message.popstr();

						json command = nlohmann::json::from_msgpack(binary_command);
						json result = executioner_callback(command);

						string binary_result;
						nlohmann::json::basic_json::to_msgpack(result, binary_result);

						send_message(bus, _m_result_, ZMQ_SNDMORE);
						send_message(bus, binary_result);
					}
				}

			}
		}
	} else {
		utility_exit_with_message("Error: got unexpected greetings from server process, exiting...\nGreetings:" + greetings);
	}

	client.join();
}


void hal(HAs _actions, CommandLineArguments const &args)
{
	std::map<string, HA> actions( std::move(_actions) );

	hal([&actions]() {
			auto f = json::object();

			for(auto const & a : actions) f[a.first] = a.second.specification;

			json spec;
			spec[_f_functions_] = f;
			return spec;
		},
		[&actions](json const &command) -> json {
			std::string name;
			if ( extract_value_if_present(command, _f_name_, name) ) {
				auto it = actions.find(name);
				if( it != actions.end() ) return it->second.execute(command[_f_arguments_]);
			}
			return json();
		},
		args);

}


} // namespace network
} // namespace protocols

#endif // defined(ZEROMQ) and  defined(_NLOHMANN_JSON_ENABLED_)
