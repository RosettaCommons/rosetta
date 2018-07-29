// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  apps/pilot/sergey/hal_demo.cc
/// @brief simple HAL protocol implementation demo, use this file as blueprint for your HAL-capable Rosetta application
/// @author Sergey Lyskov

#include <utility/json_utilities.hh>


#if defined(ZEROMQ)  and  defined(_NLOHMANN_JSON_ENABLED_)

#include <devel/init.hh>

#include <basic/options/option.hh>

#include <core/pose/Pose.hh>

#include <protocols/network/hal.hh>
#include <protocols/network/util.hh>

#include <json.hpp>

#include <basic/Tracer.hh>

#include <protocols/network/ui_mover.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>

#include <utility/pointer/memory.hh>

static basic::Tracer TR("hal-demo");


using namespace utility;
using namespace protocols::network;

using std::string;


auto const _n_set_all_phi_and_psi_ = "set all phi and psi";
auto const _n_change_phi_          = "change phi";
auto const _n_pose_from_sequence_  = "pose from sequence";

/// Generate HAL specification
json specification()
{
	json j;

	// optional, set name of your app (if no name provided then name of executables will be used): j["name"] = "hal-demo!";

	j["functions"] = json::object();


	json set_all_phi_and_psi = json::object();
	set_all_phi_and_psi[_f_pose_] = { {_f_type_, _t_pose_}, };
	set_all_phi_and_psi["phi"]    = { {_f_type_, _t_float_}, {_f_default_, 180}, {_f_description_, "value to which all phi angles will be set"}, };
	set_all_phi_and_psi["psi"]    = { {_f_type_, _t_float_}, {_f_default_, 180}, {_f_description_, "value to which all psi angles will be set"}, };
	j["functions"][_n_set_all_phi_and_psi_] = set_all_phi_and_psi;

	json change_phi = json::object();
	change_phi[_f_pose_]     = { {_f_type_, _t_pose_}, };
	change_phi["residue_n"]  = { {_f_type_, _t_integer_}, {_f_default_, 10},   {_f_min_, 1},    {_f_max_, 100}, {_f_description_, "residue number of phi angle"}, };
	change_phi["delta"]      = { {_f_type_, _t_float_},   {_f_default_, 1},    {_f_min_, -359}, {_f_max_, 359}, {_f_description_, "how much phi should be adjusted on each iteration"}, };
	change_phi["repeat"]     = { {_f_type_, _t_integer_}, {_f_default_, 2048}, {_f_min_, 1},                    {_f_description_, "how many times to repeat"}, };
	j["functions"][_n_change_phi_] = change_phi;

	json pose_from_sequence = json::object();
	pose_from_sequence["sequence"] = { {_f_type_, _t_string_}, {_f_default_, "aaaa"},  {_f_description_, "sequence to use"}, };
	j["functions"][_n_pose_from_sequence_] = pose_from_sequence;

	return j;
}

// json file_and_dir = json::object();
// file_and_dir["file"] = { {_f_type_, _t_file_}, };
// file_and_dir["dir"]  = { {_f_type_, _t_directory_}, };
// j["functions"]["file_and_dir"] = file_and_dir;

// json test_ui_mover = json::object();
// j["functions"]["test_ui_mover"] = test_ui_mover;



//j["functions"] = { {"foo_mover", foo}, {_n_zero_all_phi_, zero_all_phi} };
//change_phi.emplace("magnitude",   json::object( { {"type", "float"}, {"optional", true}, {"default", 1}, {"min", -10}, {"max", -2}, } ) );



void set_all_phi_and_psi(core::pose::Pose & pose, double phi, double psi)
{
	protocols::network::UIMover ui;

	for ( uint i=1; i <= pose.size(); ++i ) {
		pose.set_phi(i, phi);
		pose.set_psi(i, psi);
		ui.apply(pose);
	}
}

void change_phi(core::pose::Pose & pose, uint residue, double delta, int repeat)
{
	protocols::network::UIMover ui;
	//protocols::network::AddUIObserver( pose );

	if ( residue >=1  and residue <= pose.size() ) {
		for ( int i=0; i < repeat; ++i ) {

			pose.set_phi(residue, pose.phi(residue) + delta);
			ui.apply(pose);
		}
	}
}


core::pose::Pose pose_from_sequence(std::string sequence)
{
	auto pose = core::pose::Pose();

	std::transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);
	core::pose::make_pose_from_sequence(pose, sequence, core::chemical::FA_STANDARD );

	for ( uint i=1; i <= pose.size(); ++i ) pose.set_omega(i, 180);

	return pose;
}

void test_ui_mover()
{
	zmq_context().setctxopt(ZMQ_MAX_SOCKETS, 30999);

	protocols::network::UIMoverSP ui = std::make_shared<protocols::network::UIMover>();
	int i=0;
	core::pose::Pose pose;
	for ( ; ; ++i ) {
		TR << "test_ui_mover(): max_sockets:" << zmq_context().getctxopt(ZMQ_MAX_SOCKETS) << " i=" << i << std::endl; // ZMQ_SOCKET_LIMIT

		ui->apply(pose);

		ui = std::make_shared<protocols::network::UIMover>(*ui);
		//protocols::network::UIMoverSP ui = std::make_shared<protocols::network::UIMover>();

		//zmq::socket_t hal(zmq_context(), ZMQ_DEALER);
	}


	// zmq::context_t context;
	// for(;;) {
	//  zmq::socket_t s (context, ZMQ_DEALER);
	// }


	// zmq::context_t *context = new zmq::context_t();
	// zmq::socket_t  *s = new zmq::socket_t(*context, ZMQ_DEALER);

	// delete s;
	// TR << "deleting socket... Done!" << std::endl;

	// TR << "deleting context..." << std::endl;
	// delete context;
	// TR << "deleting context... Done!" << std::endl;
}


json hal_executioner(json const &command)
{
	TR << "Rosetta: executing command: " << command[_f_name_] << std::endl;

	core::pose::PoseOP pose = protocols::network::bytes_to_pose( command[_f_arguments_].value(_f_pose_, "") );

	if ( not pose ) pose = utility::pointer::make_shared<core::pose::Pose>();

	string name;
	if ( extract_value_if_present(command, _f_name_, name) ) {

		if ( pose ) {
			if ( name == _n_set_all_phi_and_psi_ ) {
				double phi = command[_f_arguments_].value("phi", 180.0);
				double psi = command[_f_arguments_].value("psi", 180.0);
				set_all_phi_and_psi(*pose, phi, psi);
			}

			if ( name == _n_change_phi_ ) {
				uint residue = command[_f_arguments_].value("residue_n", 1);
				double delta = command[_f_arguments_].value("delta", 1.0);
				int repeat   = command[_f_arguments_].value("repeat", 1);

				change_phi(*pose, residue, delta, repeat);
			}

		}
		if ( name == _n_pose_from_sequence_ ) {
			string sequence  = command[_f_arguments_].value("sequence", "");
			*pose = pose_from_sequence(sequence);
		}

		if ( name == "test_ui_mover" ) {
			test_ui_mover();
		}
	}


	auto pose_binary = protocols::network::pose_to_bytes(pose);

	nlohmann::json result;

	result[_f_pose_] = pose_binary;

	return result;
}

// protocols::network::UIMover ui;
// for ( int i=0; i<10000; ++i ) {
//  pose->set_phi(16, pose->phi(16) + 1);
//  ui.apply(*pose);
// }
//protocols::network::AddUIObserver(*pose);

int main(int argc, char * argv [])
{

	try {
		devel::init(argc, argv);

		{ // creating dummy pose object to trigger database load so later we can create Pose immeditaly
			core::pose::Pose p;
			core::import_pose::pose_from_pdbstring(p, "ATOM     17  N   ILE A   1      16.327  47.509  23.466  1.00  0.00\n");
		}

		/// deprecated, leaved as example for cases when by-hand implementation of `specification` and  `hal_executioner` is needed
		//protocols::network::hal(specification, hal_executioner, protocols::network::CommandLineArguments{argc, argv} );

		hal({
			// void set_all_phi_and_psi(core::pose::Pose & pose, double phi, double psi)
			{"set all phi and psi", {set_all_phi_and_psi, {
			{ {_f_name_, "phi"}, {_f_default_, 180}, {_f_description_, "value to which all phi angles will be set"}, }, // full specification except the type ( {_f_type_, _t_float_}, is ommitted)
			{ {_f_name_, "psi"}, {_f_default_, 180}, {_f_description_, "value to which all psi angles will be set"}, } // full specification except the type ( {_f_type_, _t_float_}, is ommitted)
			} } },

			/// Example: we can specify only name of the argument and letting compiler determent type for us
			// {"set all phi and psi 2", {set_all_phi_and_psi, {"phi", "psi"} } },

			/// Example: full specification including the type override, `psi` will be entered as interge on UI side
			// {"set all phi and psi 3", {set_all_phi_and_psi, {
			//               "phi",
			//               { {_f_name_, "psi"}, {_f_type_, _t_integer_}, {_f_default_, 180}, {_f_description_, "value to which all phi angles will be set"}, }
			//               } } },

			// {"set all phi and psi 4", {set_all_phi_and_psi, {"phi"} } }, // assign missed parameter names automatically to arg-0, arg-1, ...


			// void change_phi(core::pose::Pose & pose, uint residue, double delta, int repeat)
			{"change phi", {change_phi, {
			{ {_f_name_, "residue_n"}, {_f_default_, 10},   {_f_min_, 1},    {_f_max_, 100}, {_f_description_, "residue number where PHI will be chnaged"}, },
			{ {_f_name_, "delta"},     {_f_default_, 1},    {_f_min_, -359}, {_f_max_, 359}, {_f_description_, "how much phi should be adjusted on each iteration"}, },
			{ {_f_name_, "repeat"},    {_f_default_, 2048}, {_f_min_, 1},                    {_f_description_, "how many times to repeat"}, },
			} } },


			/// Note that `core::pose::Pose pose_from_sequence(std::string sequence)` will be automatically treated as Action that does not require Pose as an argument
			{"pose from sequence", {pose_from_sequence, {"sequence"} } },

			},
			CommandLineArguments{argc, argv} );


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
