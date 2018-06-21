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


#include <protocols/network/ui_mover.hh>


using std::string;

/// Generate HAL specification
string specification()
{
	nlohmann::json j;
	j["functions"] = {
		{{ "name", "my_mover" },
		{ "arguments",
		{
		{"param1", "int"},
		{"param2", "float"},
		},
		}},

		{{ "name", "test_mover2" },
		{ "arguments",
		{
		{"param1", "string"},
		{"param2", "float"},
		},
		}},

		};

	string r;
	nlohmann::json::basic_json::to_msgpack(j, r);
	return r;
}


json hal_executioner(json const &command)
{
	//std::cout << "Rosetta: executing command: " << command << std::endl;
	std::cout << "Rosetta: executing command: " << command[protocols::network::_f_name_] << std::endl;

	core::pose::PoseOP pose = protocols::network::bytes_to_pose(command[protocols::network::_f_arguments_][protocols::network::_f_pose_]);

	protocols::network::UIMover ui;

	for ( int i=0; i<10000; ++i ) {
		pose->set_phi(16, pose->phi(16) + 1);
		ui.apply(*pose);
	}

	auto pose_binary = protocols::network::pose_to_bytes(*pose);

	nlohmann::json result;

	result["pose"] = pose_binary;

	return result;
}


int main(int argc, char * argv [])
{

	try {
		devel::init(argc, argv);

		protocols::network::hal(specification, hal_executioner, argc, argv);
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
