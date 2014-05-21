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
#include <core/chemical/AtomType.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/Patch.hh>

// Project Headers
#include <core/types.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <utility/vector1.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/string_util.hh>
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


/// @brief  Add additional parameter files not present in <atom-set-name>/extras.txt. Called by ChemicalManager at time of AtomTypeSet creation.

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


/// @brief Modify atom_type properties from the command line. Called by ChemicalManager at time of AtomTypeSet creation.
void
modify_atom_properties_from_command_line(
	std::string const & atom_type_set_tag,
	AtomTypeSet & atom_type_set
)
{


	if ( basic::options::option[ basic::options::OptionKeys::chemical::set_atom_properties ].user() ) {
		utility::vector1< std::string > const & mods
			( basic::options::option[ basic::options::OptionKeys::chemical::set_atom_properties ]);

		std::string const errmsg( "-set_atom_properties format should be:: -set_atom_properties <set1>:<atom1>:<param1>:<setting1> <set2>:<atom2>:<param2>:<setting2> ...; for example: '-chemical:set_atom_properties fa_standard:OOC:LK_DGFREE:-5 fa_standard:ONH2:LJ_RADIUS:0.5' ");

		for ( Size i=1; i<= mods.size(); ++i ) {
			///
			/// mod should look like (for example):  "fa_standard:OOC:LK_RADIUS:4.5"
			///
			std::string const & mod( mods[i] );

			Size const pos1( mod.find(":") );
			if ( pos1 == std::string::npos ) utility_exit_with_message(errmsg);
			std::string const atomset_tag( mod.substr(0,pos1) );
			if ( atomset_tag != atom_type_set_tag ) continue;

			Size const pos2( mod.substr(pos1+1).find(":") );
			if ( pos2 == std::string::npos ) utility_exit_with_message(errmsg);
			std::string const atom_name( mod.substr(pos1+1,pos2) );
			if ( !atom_type_set.has_atom( atom_name ) ) utility_exit_with_message(errmsg+". Nonexistent atomname: "+atom_name);

			Size const pos3( mod.substr(pos1+1+pos2+1).find(":") );
			if ( pos3 == std::string::npos ) utility_exit_with_message(errmsg);
			std::string const param( mod.substr(pos1+1+pos2+1,pos3) );

			std::string const stringsetting( mod.substr(pos1+1+pos2+1+pos3+1) );
			if ( !ObjexxFCL::is_double( stringsetting ) ) utility_exit_with_message(errmsg);
			Real const setting( ObjexxFCL::double_of( stringsetting ) );

			TR.Trace << "modify_atom_properties_from_command_line: setting " << atomset_tag << ' ' << atom_name << ' ' <<
				param << ' ' << setting << std::endl;

			Size const atom_index( atom_type_set.atom_type_index( atom_name ) );

			/// I would like to uncomment the following if-check, but right now there is an extra parameter file
			/// that defines a parameter with the name LK_DGFREE (memb_fa_params.txt). That's kind of confusing...
			///
			// if ( atom_type_set.has_extra_parameter( param ) ) {
			// 	Size const param_index( atom_type_set.extra_parameter_index( param ) );
			// 	atom_type_set[ atom_index ].set_extra_parameter( param_index, setting );
			// } else {
			atom_type_set[ atom_index ].set_parameter( param, setting );
		}
	}
}

std::string
formatted_icoord_tree(core::chemical::ResidueType const & restype) {
	// General scheme - pick an atom, then move up the icoord tree until you hit the root
	// Then proceed down the tree depth first, outputting atom names as you go,
	// keeping a stack of atoms to go back to.
	// The complicated bit is keeping track of where the different attachment points are
	// such that you can accurately place the parenthesis for the tree.
	std::string output("Icoord tree: ");

	utility::vector1< bool > placed( restype.natoms(), false );
	utility::vector1< AtomICoor > icoords;
	utility::vector1< core::Size > deferred;
	for( core::Size ii(1); ii <= restype.natoms(); ++ii ) {
		icoords.push_back( restype.icoor(ii) );
	}

	core::Size curres(1); // The first atom is arbitrary (although it probably is the root.
	// The root atom has itself listed as the stub_atom1.
	while( icoords[curres].stub_atom1().atomno() != curres ) {
		curres = icoords[curres].stub_atom1().atomno();
	}

	do {
		output += restype.atom_name( curres );
		placed[curres] = true;

		utility::vector1< core::Size > possibles;
		for( core::Size ii(1); ii <= icoords.size(); ++ii ) {
			if( icoords[ii].stub_atom1().atomno() == curres && ! placed[ii] ) {
				possibles.push_back( ii );
			}
		}

		if( possibles.size() == 1 ) {
			curres = possibles[1];
		} else if ( possibles.size() ) {
			output += " (";
			curres = possibles.back();
			possibles.pop_back();
			if( deferred.size() ) {
				deferred.push_back( 0 ); // Sentinel, but not needed for first
			}
			deferred.insert( deferred.end() , possibles.begin(), possibles.end() );
		} else if( deferred.size() ) {
			output += ") ";
			while( deferred.size() && deferred.back() == 0 ) {
				output += ") ";
				deferred.pop_back();
			}
			if( ! deferred.size() ) {
				curres = 0;
				break;
			}
			output += " (";
			curres = deferred.back();
			deferred.pop_back();
		} else {
			output += ") ";
			curres = 0;
			break;
		}
	} while( curres != 0 );

	return output;
}


void print_chis(std::ostream & out, ResidueType const & res) {
	using namespace core::chemical;
	out << "Residue: " << res.name() << std::endl;
	out << formatted_icoord_tree(res) << std::endl;
	for( core::Size ii(1); ii <= res.nchi(); ++ii ) {
		AtomIndices const & indexes( res.chi_atoms(ii) );
		out << "Chi " << ii << ": " << res.atom_name( indexes[1] ) << " " << res.atom_name( indexes[2] ) << " "
				<< res.atom_name( indexes[3] ) << " " << res.atom_name( indexes[4] ) << " ";
		if( res.is_proton_chi( ii) ) { out << " PROTON";}
		out << std::endl;
	}
}

/////////////////////////////////////////////////////////
std::string
fixup_patches( std::string string_in ){
	std::string const string_out = utility::replace_in( string_in, "_p:", patch_linker );
	return string_out;
}

} // namespace chemical
} // namespace core
