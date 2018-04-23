// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/task/util.cc
/// @brief Utility functions for main task classes.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <core/pack/task/util.hh>


#include <core/sequence/sequence_motif.hh>
#include <core/pack/task/ResfileReader.hh>

#include <basic/Tracer.hh>

#include <utility/string_util.hh>

static basic::Tracer TR( "core.pack.task.util" );


namespace core {
namespace pack {
namespace task {

using namespace core::sequence;
using namespace utility;


utility::vector1< ResfileCommandOP >
parse_res_agnostic_commands(
	std::string const & line,
	std::map< std::string, ResfileCommandOP > const & command_map)
{
	bool found_commands(false);
	vector1< ResfileCommandOP > commands;

	//std::string final_line = "1 "+line; //Add bogus residue number.
	vector1< std::string > tokens = tokenize_line( line);

	Size which_token = 1, ntokens(tokens.size());
	while ( which_token <= ntokens ) {

		//Locate the command
		std::map< std::string, ResfileCommandOP >::const_iterator command(
			command_map.find(get_token(which_token, tokens)));

		if ( command == command_map.end() ) {
			utility_exit_with_message(get_token(which_token, tokens) + " Not recognized in parse_res_agnostic_commands!");
		} else { }

		//Create the actual command from tokens
		ResfileCommandOP new_command = command->second->clone();

		new_command->initialize_from_tokens( tokens, which_token, 1 );
		commands.push_back( new_command );
		found_commands = true;
	}

	if ( !found_commands ) {
		utility_exit_with_message("Could not convert command: "+line);
	}

	return commands;
}

vector1< vector1< ResfileCommandOP> >
get_resfile_commands( std::string const & motif ){

	vector1< vector1< ResfileCommandOP > > per_position_commands;
	vector1< std::string > final_split_motif = split_sequence_motif(motif);


	TR.Debug << "Final Split Motif: " << utility::to_string(final_split_motif) << std::endl;
	std::map< std::string, ResfileCommandOP > command_map = create_command_map();

	//First, we take the raw command and turn it into a resfile command if needed.
	for ( std::string const & raw_cmd : final_split_motif ) {
		utility::vector1< ResfileCommandOP > commands;
		if ( startswith(raw_cmd, "%") ) {
			std::string final_cmd = raw_cmd.substr(1, raw_cmd.size() - 1); //Remove %
			commands = parse_res_agnostic_commands(final_cmd, command_map);
			TR.Debug << "Parsing Resfile Command" << std::endl;
		} else if ( raw_cmd == "-" ) {
			ResfileCommandOP cmd = ResfileCommandOP( new NATAA );
			commands.push_back(cmd);
			TR.Debug << "Parsing - Command" << std::endl;
		} else if ( raw_cmd == "X" ) {
			ResfileCommandOP cmd = ResfileCommandOP( new ALLAA );
			commands.push_back(cmd);
			TR.Debug << "Parsing X Command" << std::endl;
		} else if ( startswith(raw_cmd, "^") ) {
			std::string final_cmd = "NOTAA "+raw_cmd.substr(1, raw_cmd.size() - 1);
			commands = parse_res_agnostic_commands(final_cmd, command_map);
			TR.Debug << "Parsing Not AA ^ Comamnd " << std::endl;
		} else {
			//List of residues as PIKAA
			std::string final_cmd = "PIKAA " + raw_cmd;
			commands = parse_res_agnostic_commands(final_cmd, command_map);
			TR.Debug << "Parsing a standard set of amino acids " << std::endl;
		}
		per_position_commands.push_back(commands);
	}

	TR.Debug << "Created commands for " << per_position_commands.size() << " positions" << std::endl;
	return per_position_commands;
}



} //core
} //pack
} //task


