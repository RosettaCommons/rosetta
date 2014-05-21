// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/chemical/sdf/MolFileIOReader.cc
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <core/chemical/sdf/MolFileIOReader.hh>
#include <core/chemical/sdf/MolFileIOData.hh>
#include <core/chemical/sdf/SDFParser.hh>

#include <core/chemical/ResidueType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ElementSet.hh>

#include <numeric/util.hh>
#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <utility/file/FileName.hh>
#include <utility/io/izstream.hh>

#include <ObjexxFCL/string.functions.hh>

#include <fstream>

namespace core {
namespace chemical {
namespace sdf {

static basic::Tracer TR("core.io.sdf.MolFileIOReader");

MolFileIOReader::MolFileIOReader()
{}

MolFileIOReader::~MolFileIOReader()
{}

/// @details This has to take a filename, as autodetection of file type may require reopening the file

utility::vector1< MolFileIOMoleculeOP > MolFileIOReader::parse_file( std::string const & filename, std::string type /*= "" */ ) {

	utility::io::izstream file( filename );
	if ( ! file ) {
		TR << "Error: cannot open file " << filename << std::endl;
		utility_exit_with_message( "Error opening " + filename );
	}
	if ( type == "" ) {
		utility::file::FileName fn( filename );
		std::string ext( fn.extension() ); // The fn.extension() machinery should seamlessly ignore the .gz ending.
		ObjexxFCL::lowercase( ext );
		if( ext == ".mol" || ext == ".sdf" || ext == ".pdb" || ext == ".mol2" || ext == ".params" ) {
			type = ext.substr(1); // Everything from the second character to the end of the string.
		} else {
			TR << "Extension of '" << ext << "' for file '" << filename << "' not recognized - attempting autodetection of type." << std::endl;
			utility::io::izstream test_type( filename );
			std::string line;
			bool M_END_line = false;
			bool dollars_line = false;
			bool ATOM_line = false;
			bool BOND_line = false;
			bool Tripos_line = false;
			for( getline(test_type, line);  test_type; getline(test_type, line) ) {
				if( utility::startswith(line, "M  END") ) { M_END_line = true; }
				if( utility::startswith(line, "$$$$") ) { dollars_line = true; }
				if( utility::startswith(line, "ATOM") ) { ATOM_line = true; }
				if( utility::startswith(line, "BOND") ) { BOND_line = true; }
				if( utility::startswith(line, "@<TRIPOS>") ) { Tripos_line = true; }
			}
			if( M_END_line || dollars_line ) {
				if( ATOM_line || BOND_line || Tripos_line ) {
					TR.Warning << "Warning: Ambiguous filetype for loading " << filename << ", assuming SDF." << std::endl;
				}
				type = "sdf";
			} else if( Tripos_line ) {
				if( ATOM_line || BOND_line || M_END_line || dollars_line ) {
					TR.Warning << "Warning: Ambiguous filetype for loading " << filename << ", assuming mol2." << std::endl;
				}
				type = "mol2";
			} else if( ATOM_line && BOND_line ) {
				if( Tripos_line || M_END_line || dollars_line ) {
					TR.Warning << "Warning: Ambiguous filetype for loading " << filename << ", assuming Rosetta params." << std::endl;
				}
				type = "params";
			} else if( ATOM_line && ! BOND_line ) {
				if( Tripos_line || M_END_line || dollars_line ) {
					TR.Warning << "Warning: Ambiguous filetype for loading " << filename << ", assuming PDB." << std::endl;
				}
				type = "pdb";
			} else {
				TR.Error << "ERROR: Unable to autodetermine filetype of molecule file " << filename << std::endl;
				utility_exit_with_message( "Can't determine filetype for file " + filename );
			}
		} // if/else extension recognized
	} // type = ""
	return parse_file( file, type ); // file will be closed on destruction.
}

utility::vector1< MolFileIOMoleculeOP > MolFileIOReader::parse_file( std::istream & file, std::string type ) {
	ObjexxFCL::lowercase( type );
	utility::vector1< MolFileIOMoleculeOP > molecules;

	if( type == "mol" || type == "sdf" ) {
		SDFParser parser;
		molecules = parser.parse( file );
	} else if( type == "mol2" ) {
		TR.Error << "Loading of mol2 files via this method not currently supported." << std::endl;
	} else if( type == "pdb" ) {
		TR.Error << "Loading of pdb files via this method not currently supported." << std::endl;
	} else if( type == "params" ) {
		TR.Error << "Loading of params files via this method not currently supported." << std::endl;
	} else {
		utility_exit_with_message( "Do not know how to handle molecule file of type " + type );
	}

	if( ! molecules.size() ) {
		TR.Error << "Error: Stream contained no recognized molecules!" << std::endl;
	}
	return molecules;
}


utility::vector1< ResidueTypeOP > convert_to_ResidueType( utility::vector1< MolFileIOMoleculeOP > molfile_data,
			std::string atom_type_tag,
			std::string elements_tag,
			std::string mm_atom_type_tag) {
	AtomTypeSetCAP atom_types = ChemicalManager::get_instance()->atom_type_set( atom_type_tag );
	ElementSetCAP elements = ChemicalManager::get_instance()->element_set( elements_tag );
	MMAtomTypeSetCAP mm_atom_types = ChemicalManager::get_instance()->mm_atom_type_set( mm_atom_type_tag );
	return convert_to_ResidueType(molfile_data, atom_types, elements, mm_atom_types);
}

utility::vector1< ResidueTypeOP > convert_to_ResidueType( utility::vector1< MolFileIOMoleculeOP > molfile_data,
			AtomTypeSetCAP atom_types,
			ElementSetCAP element_type_set,
			MMAtomTypeSetCAP mm_atom_types) {
	utility::vector1< ResidueTypeOP > restypes;
	utility::vector1< MolFileIOMoleculeOP >::iterator iter, end;
	for( iter = molfile_data.begin(), end = molfile_data.end(); iter != end; ++iter ) {
		restypes.push_back( (*iter)->convert_to_ResidueType(atom_types,
				element_type_set,
				mm_atom_types) );
	}

	return restypes;
}


}
}
}
