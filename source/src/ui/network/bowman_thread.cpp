#include <ui/network/bowman_thread.h>

#include <protocols/network/util.hh>

#include <utility/CSI_Sequence.hh>

#include <libzmq/include/zmq.hpp>

#include <iostream>

namespace ui {
namespace network {

using std::string;
using namespace utility;
using namespace protocols::network;

string as_hexadecimal(string const &s, bool as_bytes)
{
	std::ostringstream o;
	//o << s.size() << "_";

	if( as_bytes ) {
		for(unsigned char c : s) { // cast to `unsigned char` is intentional here
			o << std::hex << std::setw(2) << std::setfill('0') << static_cast<short>(c);// << ' ';
		}
	}
	else {
		for(unsigned char c : s) { // cast to `unsigned char` is intentional here
			if( c >= 33 && c < 127) o << c;
			else o << std::hex /*<< std::setw(2)*/ << std::setfill('0') << static_cast<short>(c);
		}
	}
	//o << std::resetiosflags(std::ios_base::basefield);
	return o.str();
}

BowmanThread::BowmanThread(QObject *parent) : QThread(parent)
{
	//context();
	//bus(); // making sure bus socket is created and bound _before_ thread start to execute

	context_ = std::make_shared<zmq::context_t>(1);

	bus_.reset( new zmq::socket_t(*context_, ZMQ_PAIR) );
	bus_->bind(_bus_address_);

	//connect(this, &BowmanThread::client_connected, this, &BowmanThread::test);
}


//BowmanThread::~BowmanThread() = default;
BowmanThread::~BowmanThread()
{
	send_message(*bus_, _m_quit_);
	wait();
}


void BowmanThread::run()
{
	std::cout << "BowmanThread::run() starting..." << std::endl;

	zmq::socket_t bus(*context_, ZMQ_PAIR);
	bus.connect(_bus_address_);

	zmq::socket_t hal(*context_, ZMQ_ROUTER);
	hal.bind(_server_address_);

	std::vector<zmq_pollitem_t> items = { {hal, 0, ZMQ_POLLIN, 0 }, {bus, 0, ZMQ_POLLIN, 0 } };

	auto const timeout = std::chrono::milliseconds(2500);
	auto const ping_interval = std::chrono::milliseconds(500);
	std::chrono::steady_clock::time_point last_ping_time = std::chrono::steady_clock::now() - ping_interval;

	std::map<string, std::chrono::steady_clock::time_point> clients; // client_id -> time of last message from that client

	while (true) {
		//zmq::message_t id;

		if( zmq::poll( items, std::chrono::milliseconds(500) ) > 0 ) {

			if ( items [0].revents & ZMQ_POLLIN ) { // messages from hal clients

				std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();

				//std::cout << "got message from hal: " << receive_message(hal) << std::endl;
				zmq::multipart_t message(hal);

				if( message.size() >= 2 ) {
					string hal_id = message.popstr();
					string type = message.popstr();
					//std::cout << "got message from hal with id:" << as_hexadecimal(id) << std::endl;

					auto it = clients.find(hal_id);

					if( it == clients.end() ) {
						std::cout << "BowmanThread::run(): " << CSI_bgBlue() << CSI_Black() << "Client " << as_hexadecimal(hal_id) << " connected!" << CSI_Reset() << std::endl;
						clients.emplace(hal_id, now);
						Q_EMIT client_connected(hal_id);  // NOTE: cross-thread emit

					} else it->second = now;

					//std::cout << "got message from hal: " << message.str() << std::endl;
					//id = std::move(message[0]);

					if( type == _m_specification_  and  message.size() == 1 ) {
						string specification = message.popstr();
						//std::cout << "BowmanThread::run(): " << CSI_bgBlue() << CSI_Black() << "Got specification from client " << as_hexadecimal(id) << ":" << CSI_Reset() << json::from_msgpack( Bytes(specification.begin(), specification.end() ) ) << std::endl;
						std::cout << "BowmanThread::run(): " << CSI_bgGreen() << CSI_Black() << "Got specification from client " << as_hexadecimal(hal_id) << ":" << CSI_Reset() << std::endl;

						JSON_SP j_up = std::make_shared<json>(json::from_msgpack(specification) );
						Q_EMIT specification_received(hal_id, j_up);  // NOTE: cross-thread emit

					} else if( type == _m_result_  and  message.size() == 1 ) {
						string binary_result = message.popstr();
						JSON_SP result = std::make_shared<json>(json::from_msgpack(binary_result) );
						//std::cout << "BowmanThread::run(): Got " << type << " message from hal " <<  as_hexadecimal(hal_id) /*<< " result: " << *result */<< std::endl;

						Q_EMIT result_received(result);  // NOTE: cross-thread emit

					} else if( type == _m_progress_  and  message.size() == 1 ) {
						string binary_result = message.popstr();
						JSON_SP result = std::make_shared<json>(json::from_msgpack(binary_result) );
						//std::cout << "BowmanThread::run(): Got " << type << " message from hal " <<  as_hexadecimal(hal_id) /*<< " result: " << *result */<< std::endl;

						Q_EMIT progress_data_received(result);  // NOTE: cross-thread emit

					} else if( type != _m_ping_) {
						std::cout << "BowmanThread::run(): Got '" << type << "' message from hal, type is uknown - discarding..." <<  std::endl;
					}


				} else std::cout << "BowmanThread::run(): ERROR Got unexpected number of parts in hal message:" << message.size() << "!" << std::endl;

			}

			if ( items [1].revents & ZMQ_POLLIN ) { // messages from main UI thread
				zmq::multipart_t message(bus);

				if( message.size() >= 1 ) {
					string type = message.popstr();

					if( type == _m_quit_ ) {
						std::cout << "got message from main thread: `" << type << "`, exiting..."  << std::endl;
						break;
					}
					else if( type == _m_abort_  and  message.size() == 1 ) {
						string hal_id = message.popstr();
						std::cout << "got message from main thread: `" << type << "`, sending `abort` signal to " << as_hexadecimal(hal_id) << "..."  << std::endl;
						send_message(hal, hal_id, ZMQ_SNDMORE);
						send_message(hal, _m_abort_);
					}
					else if( type == _m_execute_  and  message.size() == 2 ) {
						string hal_id = message.popstr();
						string command = message.popstr();
						//std::cout << "BowmanThread::run(): Got " << _m_execute_ << " message for hal " <<  as_hexadecimal(hal_id) << std::endl;

						send_message(hal, hal_id, ZMQ_SNDMORE);
						send_message(hal, _m_execute_, ZMQ_SNDMORE);
						send_message(hal, command);

					} else {
					}

				} else std::cout << "BowmanThread::run(): ERROR Got unexpected number of parts in bus message:" << message.size() << "!" << std::endl;
			}
		}

		// if ( --ping_count == 0 ) {
		// 	std::cout << "UI server seems to be dead, re-connecting..." << std::endl;

		// 	ui.release();
		// 	ui = create_ui_socket(*context);
		// 	ping_count = max_ping_count;
		// }

		std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
		if( clients.size()  and  now - last_ping_time > ping_interval ) { // std::chrono::duration_cast<std::chrono::milliseconds>(now - last_ping_time) > ping_interval
			last_ping_time = now;

			for (auto it = clients.begin(); it != clients.end(); ) {
				if( now - it->second > timeout ) { // std::chrono::duration_cast<std::chrono::milliseconds>(now - it->second) > timeout
					std::cout << "BowmanThread::run(): " << CSI_bgRed() << CSI_Black() << "Client " << as_hexadecimal(it->first) << " has disconnected!" << CSI_Reset() << std::endl;

					Q_EMIT client_disconnected(it->first); // NOTE: this is cross-thread emit

					it = clients.erase(it);
				} else ++it;
			}

			for(auto & it : clients) {
				//std::cout << "Sending ping to " << it.first << std::endl;
				send_message(hal, it.first, ZMQ_SNDMORE);
				send_message(hal, _m_ping_);
			}
		}

		//std::cout << ping_count;
	}

	std::cout << "BowmanThread::run() exiting..." << std::endl;
}


// BowmanThread::ContextSP BowmanThread::context()
// {
// 	static ContextSP context = std::make_shared<zmq::context_t>(1);
// 	return context;
// }
// BowmanThread * BowmanThread::get_instance()
// {
// 	static std::unique_ptr<BowmanThread> bowman( new BowmanThread() );
// 	if( not bowman->isRunning() ) bowman->start();
// 	return bowman.get();
// }

// zmq::socket_t & BowmanThread::bus()
// {
// 	static SocketUP bus;
// 	if( not bus) {
// 		bus.reset( new zmq::socket_t(*context(), ZMQ_PAIR) );
// 		bus->bind(_bus_address_);
// 	}
// 	return *bus;
// }


//
// below is functions executed in main (GUI) thread
//

void BowmanThread::abort(std::string const & hal_id)
{
	//std::cout << "BowmanThread::abort()..." << std::endl;
	send_message(*bus_, _m_abort_, ZMQ_SNDMORE);
	send_message(*bus_, hal_id);
}

void BowmanThread::execute(std::string const & hal_id, JSON_CSP const & command)
{
	send_message(*bus_, _m_execute_, ZMQ_SNDMORE);
	send_message(*bus_, hal_id, ZMQ_SNDMORE);

	string binary_command;
	nlohmann::json::basic_json::to_msgpack(*command, binary_command);
	send_message(*bus_, binary_command);
}

//#include <QCoreApplication>
// void BowmanThread::test()
// {
// 	bool is_main_thread = QThread::currentThread() == QCoreApplication::instance()->thread();
// 	bool is_bowman_thread = QThread::currentThread() == this;
// 	std::cout << "BowmanThread::test(): is_main_thread = " << std::boolalpha << is_main_thread << std::endl;
// 	std::cout << "BowmanThread::test(): is_bowman_thread = " << std::boolalpha << is_bowman_thread << std::endl;
// }

} // namespace network
} // namespace ui
