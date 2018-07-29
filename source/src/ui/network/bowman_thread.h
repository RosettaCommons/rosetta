#pragma once

#include <QThread>

#include <protocols/network/util.hh>

#include <json.hpp>
#include <libzmq/include/zmq_addon.hpp>

namespace ui {
namespace network {

class BowmanThread final : public QThread
{
    Q_OBJECT

private:
	using ContextSP = std::shared_ptr<zmq::context_t>;
	using SocketUP  = std::unique_ptr<zmq::socket_t>;
	using SocketSP  = std::shared_ptr<zmq::socket_t>;

private:

public:
	BowmanThread(QObject *parent = 0);
    ~BowmanThread();

	//static BowmanThread * get_instance();

	void abort(std::string const & hal_id);
	void execute(std::string const & hal_id, JSON_CSP const & command);

	void set_hal_settings(protocols::network::HAL_Settings const &);

public Q_SLOTS:

	//void test();

Q_SIGNALS:

	// cross-thread signals
	void client_connected(std::string const &);
	void client_disconnected(std::string const &);

	void specification_received(std::string const &, JSON_CSP const &);

	void result_received(JSON_SP const &);
	void progress_data_received(JSON_SP const &);

private:
	void run() override;

	//static zmq::socket_t & bus();
	//static ContextSP context();

	// shared data between main-thread and network-thread
	ContextSP context_;

	// main-thread data
	SocketUP bus_;
};

std::string as_hexadecimal(std::string const &s, bool as_bytes = true);


} // namespace network
} // namespace ui
