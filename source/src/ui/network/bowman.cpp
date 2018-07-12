
#include <ui/network/bowman.h>

#include <ui/network/argument.h>

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

/// lookup for function with given name and if found create Arguments for it, return nullptr if no function could be found
ArgumentsSP get_function_arguments(Bowman const &bowman, std::string const & function_name)
{
	for(auto const & back_end_it : bowman ) {
		Bowman::Functions const & functions(back_end_it.second);

		auto it = functions.find(function_name);
		if( it != functions.end() ) {
			ArgumentsSP r = std::make_shared<Arguments>();

			for(auto const & a : it->second) r->emplace( std::make_pair(a.first, a.second->clone()) );

			return r;
		}
	}

	return ArgumentsSP();
}


/// ----------------
/// Bowman
/// ----------------
Bowman & Bowman::bowman()
{
	static std::unique_ptr<Bowman> b( new Bowman() );

	return *b;
}

Bowman::~Bowman() = default;
// Bowman::~Bowman()
// {
// 	qDebug() << "Bowman::~Bowman()";
// }

Bowman::Bowman(QObject *parent) : QObject(parent)
{
	bowman_thread_.start();
	connect(&bowman_thread_, &BowmanThread::client_connected,       this, &Bowman::client_connected);

	// ← we can't simply re-emit the signal because we want to update inner specification list first, so `specification_received` and `client_disconnected` will be emited manually from on_bowman_thread_<event> handlers
	//connect(&bowman_thread_, &BowmanThread::client_disconnected,    this, &Bowman::client_disconnected);
	//connect(&bowman_thread_, &BowmanThread::specification_received, this, &Bowman::specification_received);

	connect(&bowman_thread_, &BowmanThread::result_received,        this, &Bowman::on_bowman_thread_result_received);
	connect(&bowman_thread_, &BowmanThread::progress_data_received, this, &Bowman::on_bowman_thread_progress_data_received);

	//QMetaObject::connectSlotsByName(this);

	connect(&bowman_thread_, &BowmanThread::specification_received, this, &Bowman::on_bowman_thread_specification_received);
	connect(&bowman_thread_, &BowmanThread::client_disconnected,    this, &Bowman::on_bowman_thread_client_disconnected);
}


void Bowman::on_bowman_thread_specification_received(std::string const &client_id, JSON_CSP const &specification)
{
	//qDebug() << "Bowman::on_bowman_thread_specification_received: id=" << client_id.c_str() << "\n" << specification->dump(2).c_str(); //QString::fromStdString(j->dump(2));

	auto it_functions = specification->find(_f_functions_);

	if( it_functions != specification->end()  and  it_functions->is_object() ) {
		auto it_and_bool = back_ends_.emplace(client_id, Functions());
		Functions & functions = it_and_bool.first->second;

		//qDebug() << "Bowman::on_bowman_thread_specification_received: it_functions:" << it_functions->dump(2).c_str(); //QString::fromStdString(j->dump(2));

		for (auto f = it_functions->begin(); f != it_functions->end(); ++f) {
			//qDebug() << "Bowman::on_bowman_thread_specification_received: adding function:" << f.key().c_str();

			json const & args = f.value();
			if( args.is_object() ) {
				auto it_and_bool = functions.emplace(f.key(), FunctionSignature());
				FunctionSignature &signature = it_and_bool.first->second;

				for( auto arg = args.begin(); arg != args.end(); ++arg) {
					ArgumentSP a = Argument::create( arg.value() );
					if(a) {
						//signature.insert( arg.key(), a );
						signature.emplace( std::make_pair(arg.key(), a) );
						//qDebug() << "Bowman::on_bowman_thread_specification_received:             arg:" << arg.key().c_str() << *a;
					}
				}
				// for(auto const & args : *f ) {
				// 	qDebug() << "Bowman::on_bowman_thread_specification_received: " << args.dump(2).c_str(); //QString::fromStdString(j->dump(2));
				// }
			}
		}

		//qDebug() << "Bowman::on_bowman_thread_specification_received:" << it.is_object();
		// if( it.count(_f_name_) ) {
		// 	std::string name = it[_f_name_];
		// 	functions_.emplace_back( FunctionIdentifier{ name, back_end_it->first } );
		// }
	}

	Q_EMIT specification_received(client_id, specification);
	Q_EMIT back_ends_changed(this);


	// for( auto const & b : back_ends_) {
	// 	for( auto const &f : b.second )
	// 		qDebug() << "Bowman::on_bowman_thread_specification_received: back-end: " << as_hexadecimal(b.first).c_str() <<  " function:" << f.first.c_str();
	// }

	// bool is_main_thread = QThread::currentThread() == QCoreApplication::instance()->thread();
	// bool is_bowman_thread = QThread::currentThread() == &bowman_thread_;
	// std::cout << "Bowman::on_bowman_thread_specification_received: is_main_thread = " << std::boolalpha << is_main_thread << std::endl;
	// std::cout << "Bowman::on_bowman_thread_specification_received: is_bowman_thread = " << std::boolalpha << is_bowman_thread << std::endl;
}

void Bowman::on_bowman_thread_client_disconnected(std::string const &id)
{
	auto it = back_ends_.find(id);
	if( it != back_ends_.end() ) back_ends_.erase(it);

	Q_EMIT client_disconnected(id);
	Q_EMIT back_ends_changed(this);
}



/// Initiate execution of given command on `back_end`.
/// If no command is specified, - execute `null` command
/// If no back_end specified, - execute of first avalible back_end
void Bowman::execute(core::pose::PoseCOP const &pose, std::string const & function_name, Arguments const &args, std::string back_end)
{
	JSON_SP j = std::make_shared<nlohmann::json>();
	(*j)[_f_name_] = function_name;
	// (*j)[_f_arguments_] = {
	// 					   {_f_pose_, pose_binary},
	// };

	(*j)[_f_arguments_] = json::object();

	for(auto const & a : args) (*j)[_f_arguments_][a.first] = a.second->encode();

	if(pose) (*j)[_f_arguments_][_f_pose_] = "..."; // this is only for make pretty-printing line below to display that `pose` will be indeed one of the arguments of this call

	qDebug() << "Bowman::execute: calling as:" << j->dump(2).c_str();

	if(pose) {
		auto pose_binary = protocols::network::pose_to_bytes(pose);
		(*j)[_f_arguments_][_f_pose_] = pose_binary;
	}

	if( back_end.empty() ) {
		for( auto const & b : back_ends_) {
			//qDebug() << "Bowman::execute: back-end: " << as_hexadecimal(b.first).c_str();
			if( b.second.find(function_name) != b.second.end() ) {
				//if( b.second.count(function_name) ) {
				back_end = b.first;
				//qDebug() << "Bowman::execute: selecting back-end" << as_hexadecimal(back_end).c_str() << " for executing" << function_name.c_str();
				break;
			}
		}
	}
	execute(j, back_end);
}

/// Initiate execution of given command on `back_end`.
/// If no back_end specified, - execute of first avalible back_end
void Bowman::execute(JSON_CSP const &command, std::string back_end)
{
	if( back_ends_.empty() ) return;

	if( back_end.empty() ) back_end = back_ends_.begin()->first;
	else {
		auto it = back_ends_.find(back_end);
		if( it == back_ends_.end() ) return;
	}

	bowman_thread_.execute(back_end, command);
}


/// Send an `abort` signal to back-end, if no back-end specified then signal sent to _all_ connected clients
void Bowman::abort(std::string back_end)
{
	for( auto const & b : back_ends_) {
		if( back_end.empty()  or  b.first == back_end ) bowman_thread_.abort(b.first);
	}
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
