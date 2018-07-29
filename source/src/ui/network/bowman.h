#pragma once

#include <ui/network/bowman_thread.h>
#include <ui/network/bowman.h>

#include <ui/network/argument.fwd.h>


#include <core/pose/Pose.fwd.hh>

#include <utility/json_utilities.hh>

namespace ui {
namespace network {

class Bowman;

struct FunctionID {
	std::string name, hal_id;
};

/// lookup for function with given name and if found create Arguments for it, return nullptr if no function could be found
ArgumentsSP get_function_arguments(Bowman const &, FunctionID const &);


///
/// Bowman representation in main (GUI) thread
/// store list of currently connected hal's clients and they specifications
///
class Bowman final : public QObject
{
    Q_OBJECT

	Bowman(QObject *parent = nullptr);

public:
	using FunctionSignature = std::map<std::string, ArgumentCSP>;
	//using FunctionSignatureSP  = std::shared_ptr< FunctionSignature >;
	//using FunctionSignatureCSP = std::shared_ptr< FunctionSignature const >;

	using Functions = std::map<std::string, FunctionSignature>;

	struct BackEnd {
		std::string name;
		Functions functions;
	};

	using BackEnds  = std::map<std::string, BackEnd>;  // client_id --> back_end

    ~Bowman();

	/// Only one inctance of Bowman object should exist per GUI application, use bowman() to get reference to that instance
	/// Note: you should only call this _after_ QCoreApplication instance got constructed and _before_ it got destroyed
	static Bowman & bowman();


	/// Initiate execution of given command on `back_end`.
	/// If given PoseOP is empty - do not add `pose` as an argument to command
	/// If no command.name is specified, - execute `null` command
	/// If no command.back_end specified, - execute of first avalible back_end
	void execute(core::pose::PoseCOP const &, FunctionID const & command, Arguments const &args = Arguments());

	/// Initiate execution of given command on `back_end`.
	/// If no back_end specified, - execute of first avalible back_end
	void execute(JSON_CSP const &command, std::string back_end="");


	/// Send an `abort` signal to back-end, if no back-end specified then signal sent to _all_ connected clients
	void abort(std::string back_end="");


	/// BackEnds access
	BackEnds const & back_ends() const { return back_ends_; }
	int size() const { return back_ends_.size(); }
	BackEnds::const_iterator begin() const { return back_ends_.begin(); }
	BackEnds::const_iterator end() const   { return back_ends_.end(); }


	/// access to various HAL options
	void pause(bool);

Q_SIGNALS:
	void client_connected(std::string const &);
	void client_disconnected(std::string const &);
	void specification_received(std::string const &, JSON_CSP const &);

	void result_received(core::pose::PoseOP const &, JSON_CSP const &);
	void progress_data_received(core::pose::PoseOP const &, JSON_CSP const &);

	void back_ends_changed(Bowman const *);

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

	protocols::network::HAL_Settings hal_settings_;
	ui::network::BowmanThread bowman_thread_;
};



// enum class State {
// 	ok         = 0,
// 	error      = 1
// };


} // namespace network
} // namespace ui
