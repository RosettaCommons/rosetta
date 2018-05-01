#pragma once


#include <QThread>

#include <libzmq/include/zmq_addon.hpp>


namespace ui {
namespace network {


class BowmanThread : public QThread
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

Q_SIGNALS:
	void client_connected(quint64);
	void client_disconnected(quint64);


private:
	void run() override;

	//static zmq::socket_t & bus();
	//static ContextSP context();

	// shared data between main-thread and network-thread
	ContextSP context_;

	// main-thread data
	SocketUP bus_;
};

} // namespace network
} // namespace ui
