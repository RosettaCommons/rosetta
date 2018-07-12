// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  apps/pilot/awatkins/hal_rna_denovo.cc
/// @brief A first pass at running RNA fragment assembly in a HAL-aware fashion. Can only do
/// single PDB input for now.
/// @author Andy Watkins

#include <utility/json_utilities.hh>


#if defined(ZEROMQ)  and  defined(_NLOHMANN_JSON_ENABLED_)

#include <devel/init.hh>

#include <basic/options/option.hh>

#include <core/pose/Pose.hh>

#include <protocols/network/hal.hh>
#include <protocols/network/util.hh>

#include <json.hpp>

#include <basic/Tracer.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include <protocols/network/ui_mover.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/rna/denovo/RNA_DeNovoProtocol.hh>
#include <core/import_pose/RNA_DeNovoSetup.hh>
#include <core/import_pose/options/RNA_DeNovoProtocolOptions.hh>

static basic::Tracer TR("hal-demo");

using namespace utility;
using namespace protocols::network;

using std::string;


auto const _n_single_loop_modeling_ = "single_loop_modeling";


/// Generate HAL specification
string specification()
{
	nlohmann::json j;
	j["functions"] = json::object();


	json single_loop_modeling = json::object();
	single_loop_modeling["pdb_fname"] = { {"type", _t_string_},   {"default", "1ycr.pdb"}, };
	single_loop_modeling["full_target_sequence"] = { {"type", _t_string_},   {"default", "gcaa"}, };
	j["functions"][_n_single_loop_modeling_] = single_loop_modeling;

	string r;
	nlohmann::json::basic_json::to_msgpack(j, r);
	return r;
}



//j["functions"] = { {"foo_mover", foo}, {_n_zero_all_phi_, zero_all_phi} };
//change_phi.emplace("magnitude",   json::object( { {"type", "float"}, {"optional", true}, {"default", 1}, {"min", -10}, {"max", -2}, } ) );



void single_loop_modeling(core::pose::Pose & pose, string const & pdb_fname, string const & full_target_sequence)
{
	//protocols::network::UIObserver ui;

	// Do an RNA_denovo run with the ui mover set up to this sequence...
	{
		using namespace protocols::rna::denovo;
		utility::options::OptionKeyList opts;
		core::import_pose::options::RNA_DeNovoProtocolOptions::list_options_read( opts );

		using namespace basic::options::OptionKeys;
		using namespace core::import_pose;
		// TODO: enable provision of sequence and secstruct file paths.
		utility::options::OptionCollection faux_cl_opts = basic::options::option;
		faux_cl_opts[ in::file::s ].set_cl_value( pdb_fname );
		faux_cl_opts[ rna::denovo::minimize_rna ].set_cl_value( "true" );
		faux_cl_opts[ rna::denovo::sequence ].set_cl_value( full_target_sequence );
		faux_cl_opts[ rna::denovo::secstruct ].set_cl_value( string( full_target_sequence.size(), '.' ) );
		RNA_DeNovoSetup rna_de_novo_setup;
		rna_de_novo_setup.initialize_from_options( faux_cl_opts );
		pose = *rna_de_novo_setup.pose();

		protocols::network::AddUIObserver( pose );
		//ui.attach(pose);//rna_de_novo_setup.pose());
		RNA_DeNovoProtocol protocol( rna_de_novo_setup.options(), rna_de_novo_setup.rna_params() );

		protocol.apply( pose );

	}

}

json hal_executioner(json const &command)
{
	TR << "Rosetta: executing command: " << command[_f_name_] << std::endl;

	//core::pose::PoseOP pose = protocols::network::bytes_to_pose(command[_f_arguments_][_f_pose_]);
	core::pose::PoseOP pose = core::pose::PoseOP( new core::pose::Pose );

	string name;
	if ( extract_value_if_present(command, _f_name_, name) ) {
		if ( name == _n_single_loop_modeling_ ) {
			// This isn't actually compatible with no-default.
			string const & pdb_fname = command[_f_arguments_].value("pdb_fname", "1ycr.pdb");
			string const & full_target_sequence = command[_f_arguments_].value("full_target_sequence", "gcaa");
			// Still have to pass a pose to the function.
			single_loop_modeling(*pose, pdb_fname, full_target_sequence);
		}
	}

	auto pose_binary = protocols::network::pose_to_bytes(*pose);

	nlohmann::json result;

	result["pose"] = pose_binary;

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

		protocols::network::hal(specification, hal_executioner, protocols::network::CommandLineArguments{argc, argv} );

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
