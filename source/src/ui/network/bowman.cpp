
#include <ui/network/bowman.h>

#include <protocols/network/util.hh>

//#include <json.hpp>


#include <QDebug>
//#include <QCoreApplication>

// qDebug is not thread safe so we using raw std::cout instead
//#include <iostream>

using namespace ui::network;
using namespace protocols::network;


namespace ui {
namespace network {

Bowman::Bowman(QObject *parent) : QObject(parent)
{
	bowman_thread_.start();
	connect(&bowman_thread_, &BowmanThread::client_connected,       this, &Bowman::client_connected);
	connect(&bowman_thread_, &BowmanThread::client_disconnected,    this, &Bowman::client_disconnected);
	//connect(&bowman_thread_, &BowmanThread::specification_received, this, &Bowman::specification_received);  ← we can't simply re-emit the signal because we want to update inner specification list first, so `specification_received` will be emited manually from on_bowman_thread_specification_received
	connect(&bowman_thread_, &BowmanThread::result_received,        this, &Bowman::on_bowman_thread_result_received);
	connect(&bowman_thread_, &BowmanThread::progress_data_received, this, &Bowman::on_bowman_thread_progress_data_received);

	//QMetaObject::connectSlotsByName(this);

	connect(&bowman_thread_, &BowmanThread::specification_received, this, &Bowman::on_bowman_thread_specification_received);
}

void Bowman::on_bowman_thread_specification_received(std::string const &client_id, JSON_CSP const &specification)
{
	back_ends_.emplace(client_id, specification);

	Q_EMIT specification_received(client_id, specification);

	// bool is_main_thread = QThread::currentThread() == QCoreApplication::instance()->thread();
	// bool is_bowman_thread = QThread::currentThread() == &bowman_thread_;
	// std::cout << "Bowman::on_bowman_thread_specification_received: is_main_thread = " << std::boolalpha << is_main_thread << std::endl;
	// std::cout << "Bowman::on_bowman_thread_specification_received: is_bowman_thread = " << std::boolalpha << is_bowman_thread << std::endl;
}

void Bowman::on_bowman_thread_client_disconnected(std::string const &id)
{
	auto it = back_ends_.find(id);
	if( it != back_ends_.end() ) back_ends_.erase(it);
}



/// Initiate execution of given command on `back_end`.
/// If no command is specified, - execute `null` command
/// If no back_end specified, - execute of first avalible back_end
void Bowman::execute(core::pose::Pose const &pose, std::string const & command_name, std::string back_end)
{
	auto pose_binary = protocols::network::pose_to_bytes(pose);

	JSON_SP j = std::make_shared<nlohmann::json>();
	(*j)[_f_name_] = command_name;
	(*j)[_f_arguments_] = {
						   {_f_pose_, pose_binary},
	};

	execute(j, back_end);
}

/// Initiate execution of given command on `back_end`.
/// If no back_end specified, - execute of first avalible back_end
void Bowman::execute(JSON_CSP const &command, std::string back_end)
{
	if( back_ends_.empty() ) return;

	std::string chosen_back_end;
	if( back_end.empty() ) chosen_back_end = back_ends_.begin()->first;
	else {
		auto it = back_ends_.find(back_end);
		if( it == back_ends_.end() ) return;
		else chosen_back_end = back_ends_.begin()->first;
	}

	bowman_thread_.execute(chosen_back_end, command);
}


void Bowman::on_bowman_thread_result_received(JSON_SP const &_result)
{
	json & result = *_result;

	core::pose::PoseOP pose;
	auto it_pose = result.find(_f_pose_);
	if( it_pose != result.end() ) {
		pose = protocols::network::bytes_to_pose( it_pose.value() );
		result.erase(it_pose);
	}
	Q_EMIT result_received(pose, _result);
}


void Bowman::on_bowman_thread_progress_data_received(JSON_SP const &_result)
{
	json & result = *_result;

	core::pose::PoseOP pose;
	auto it_pose = result.find(_f_pose_);
	if( it_pose != result.end() ) {
		pose = protocols::network::bytes_to_pose( it_pose.value() );
		result.erase(it_pose);
	}
	Q_EMIT progress_data_received(pose, _result);
}



// void Bowman::on_bowman_thread_client_connected(std::string const &)
// {
// 	//qDebug() << "Bowman::on_bowman_thread__client_connected(...)";
// }


} // namespace network
} // namespace ui
