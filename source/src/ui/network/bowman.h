#pragma once

#include <ui/network/bowman_thread.h>

#include <core/pose/Pose.fwd.hh>

#include <json.hpp>

namespace ui {
namespace network {

///
/// Bowman representation in main (GUI) thread
/// store list of currently connected hal's clients and they specifications
///
class Bowman final : public QObject
{
    Q_OBJECT

public:
	using BackEnds = std::map<std::string, JSON_CSP>;

	Bowman(QObject *parent = nullptr);
    //~Bowman();

	/// Initiate execution of given command on `back_end`.
	/// If no command is specified, - execute `null` command
	/// If no back_end specified, - execute of first avalible back_end
	void execute(core::pose::Pose const &, std::string const & command_name="", std::string back_end="");

	/// Initiate execution of given command on `back_end`.
	/// If no back_end specified, - execute of first avalible back_end
	void execute(JSON_CSP const &command, std::string back_end="");


	/// BackEnds access
	int size() const { return back_ends_.size(); }
	BackEnds::const_iterator begin() const { return back_ends_.begin(); }
	BackEnds::const_iterator end() const   { return back_ends_.end(); }


Q_SIGNALS:
	void client_connected(std::string const &);
	void client_disconnected(std::string const &);
	void specification_received(std::string const &, JSON_CSP const &);

	void result_received(core::pose::PoseOP const &, JSON_CSP const &);
	void progress_data_received(core::pose::PoseOP const &, JSON_CSP const &);

	//void result_pose_received(core::pose::PoseOP const &);
	//void progress_pose_received(core::pose::PoseOP const &);


private Q_SLOTS:
	void on_bowman_thread_specification_received(std::string const &, JSON_CSP const &);
	void on_bowman_thread_client_disconnected(std::string const &);

	void on_bowman_thread_result_received(JSON_SP const &);
	void on_bowman_thread_progress_data_received(JSON_SP const &);

	//void on_bowman_thread_client_connected(std::string const &);

private:
	BackEnds back_ends_;

	ui::network::BowmanThread bowman_thread_;
};


// enum class State {
// 	ok         = 0,
// 	error      = 1
// };


} // namespace network
} // namespace ui
