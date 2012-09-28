// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/chemical/util.cc
/// @brief Utilities for modifying and utilizing Residues and other core::chemical classes.


// Unit headers
#include <core/chemical/util.hh>

// Package Headers
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomTypeSet.hh>

// Project Headers
#include <core/types.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <utility/vector1.hh>
#include <utility/file/file_sys_util.hh>
#include <basic/Tracer.hh>

namespace core {
namespace chemical {

static basic::Tracer TR("core.chemical.util");

core::chemical::ResidueTypeSetCAP
rsd_set_from_cmd_line() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;

	std::string const type_set_name( option[ in::file::residue_type_set ]() );
	ResidueTypeSetCAP set = ChemicalManager::get_instance()->residue_type_set(
		type_set_name
	);

	return set;
}



void
add_atom_type_set_parameters_from_command_line(
																							 std::string const & atom_type_set_tag,
																							 AtomTypeSet & atom_type_set
																							 )
{
	if ( !( basic::options::option[ basic::options::OptionKeys::chemical::add_atom_type_set_parameters ].user() ) ) {
		/// do nothing if flag not present
		return;
	}
	utility::vector1< std::string > paramstring
		( basic::options::option[ basic::options::OptionKeys::chemical::add_atom_type_set_parameters ]() );
	if ( paramstring.size()%2 != 0 ) {
		utility_exit_with_message("bad format for -add_atom_type_set_parameters; should be: -add_atom_type_set_parameters <tag1> <filename1> <tag2> <filename2> ...");
	}
	Size const nparams( paramstring.size()/2 );
	for ( Size i=0; i< nparams; ++i ) {
		std::string const tag( paramstring[2*i+1] ), filename( paramstring[2*i+2] );
		TR.Trace << "add_atom_type_set_parameters_from_command_line: desired_tag= " << atom_type_set_tag <<
			" cmdline-tag= " << tag << " filename= " << filename << std::endl;
		if ( tag == atom_type_set_tag ) {
			if ( !utility::file::file_exists( filename ) ) {
				utility_exit_with_message("unable to locate/open file: "+filename );
			}
			TR.Trace << "add_atom_type_set_parameters_from_command_line: tag= " << tag << " filename= " << filename <<
				std::endl;
			atom_type_set.add_parameters_from_file( filename );
		}
	}
}

} // namespace chemical
} // namespace core
