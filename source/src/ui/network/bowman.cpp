
#include <ui/network/bowman.h>

#include <protocols/network/util.hh>

#include <utility/CSI_Sequence.hh>

#include <libzmq/include/zmq.hpp>

#include <json.hpp>


// qDebug is not thread safe so we using raw std::cout instead
#include <iostream>

namespace ui {
namespace network {

using std::string;
using namespace utility;
using namespace protocols::network;

string as_hexadecimal(string const &s)
{
	std::ostringstream o;
	//o << s.size() << "_";

	for(unsigned char c : s) { // cast to `unsigned char` is intentional here
		if( c >= 33 && c < 127) o << c;
		else o << std::hex /*<< std::setw(2)*/ << std::setfill('0') << static_cast<short>(c);
	}
	//o << std::resetiosflags(std::ios_base::basefield);
	return o.str();
}

quint64 as_number(string const &s)
{
	quint64 r = 0;
	for(unsigned char c : s) r += c*256; // cast to `unsigned char` is intentional here
	return r;
}


BowmanThread::BowmanThread(QObject *parent) : QThread(parent)
{
	//context();
	//bus(); // making sure bus socket is created and bound _before_ thread start to execute

	context_ = std::make_shared<zmq::context_t>(1);

	bus_.reset( new zmq::socket_t(*context_, ZMQ_PAIR) );
	bus_->bind(_bus_address_);
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
					string id = message.popstr();
					string type = message.popstr();
					//std::cout << "got message from hal with id:" << as_hexadecimal(id) << std::endl;

					auto it = clients.find(id);

					if( it == clients.end() ) {
						std::cout << "BowmanThread::run(): " << CSI_bgGreen() << CSI_Black() << "Client " << as_hexadecimal(id) << " connected!" << CSI_Reset() << std::endl;
						clients.emplace(id, now);
						Q_EMIT client_connected( as_number(id) ); // NOTE: this is cross-thread emit
					} else it->second = now;

					//std::cout << "got message from hal: " << message.str() << std::endl;
					//id = std::move(message[0]);

					//ping_count = max_ping_count;
					if( message.size() == 1  and  type == _m_specification_ ) {
						string specification = message.popstr();
						std::cout << "BowmanThread::run(): " << CSI_bgBlue() << CSI_Black() << "Got specification from client " << as_hexadecimal(id) << ":" << CSI_Reset() << json::from_msgpack( Bytes(specification.begin(), specification.end() ) ) << std::endl;

					}

				} else std::cout << "BowmanThread::run(): ERROR Got unexpected number of parts in multipart message:" << message.size() << "!" << std::endl;

			}

			if ( items [1].revents & ZMQ_POLLIN ) { // messages from main UI thread
				std::cout << "got message from main thread: " << receive_message(bus) << std::endl;
				break;
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

					Q_EMIT client_disconnected( as_number(it->first) ); // NOTE: this is cross-thread emit

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



} // namespace network
} // namespace ui
