// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/chemical/sdf/MolFileIOReader.cc
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <core/chemical/sdf/MolFileIOReader.hh>
#include <core/chemical/sdf/MolFileIOData.hh>
#include <core/chemical/sdf/SDFParser.hh>

#include <core/chemical/ResidueType.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/chemical/rotamers/StoredRotamerLibrarySpecification.hh>

#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <utility/file/FileName.hh>
#include <utility/io/izstream.hh>

#include <ObjexxFCL/string.functions.hh>

namespace core {
namespace chemical {
namespace sdf {

static basic::Tracer TR( "core.io.sdf.MolFileIOReader" );

MolFileIOReader::MolFileIOReader() = default;

MolFileIOReader::~MolFileIOReader() = default;

/// @details This has to take a filename, as autodetection of file type may require reopening the file

utility::vector1< MolFileIOMoleculeOP > MolFileIOReader::parse_file( std::string const & filename, std::string type /*= "" */, core::Size n_entries /* = 0 */ ) {

	utility::io::izstream file( filename );
	if ( ! file ) {
		TR.Error << "cannot open file " << filename << std::endl;
		utility_exit_with_message( "Error opening " + filename );
	}
	if ( type == "" ) {
		utility::file::FileName fn( filename );
		std::string ext( fn.extension() ); // The fn.extension() machinery should seamlessly ignore the .gz ending.
		ObjexxFCL::lowercase( ext );
		if ( ext == ".mol" || ext == ".sdf" || ext == ".pdb" || ext == ".mol2" || ext == ".params" ) {
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
			for ( getline(test_type, line);  test_type; getline(test_type, line) ) {
				if ( utility::startswith(line, "M  END") ) { M_END_line = true; }
				if ( utility::startswith(line, "$$$$") ) { dollars_line = true; }
				if ( utility::startswith(line, "ATOM") ) { ATOM_line = true; }
				if ( utility::startswith(line, "BOND") ) { BOND_line = true; }
				if ( utility::startswith(line, "@<TRIPOS>") ) { Tripos_line = true; }
			}
			if ( M_END_line || dollars_line ) {
				if ( ATOM_line || BOND_line || Tripos_line ) {
					TR.Warning << "Ambiguous filetype for loading " << filename << ", assuming SDF." << std::endl;
				}
				type = "sdf";
			} else if ( Tripos_line ) {
				if ( ATOM_line || BOND_line || M_END_line || dollars_line ) {
					TR.Warning << "Ambiguous filetype for loading " << filename << ", assuming mol2." << std::endl;
				}
				type = "mol2";
			} else if ( ATOM_line && BOND_line ) {
				if ( Tripos_line || M_END_line || dollars_line ) {
					TR.Warning << "Ambiguous filetype for loading " << filename << ", assuming Rosetta params." << std::endl;
				}
				type = "params";
			} else if ( ATOM_line && ! BOND_line ) {
				if ( Tripos_line || M_END_line || dollars_line ) {
					TR.Warning << "Ambiguous filetype for loading " << filename << ", assuming PDB." << std::endl;
				}
				type = "pdb";
			} else {
				TR.Error << "Unable to autodetermine filetype of molecule file " << filename << std::endl;
				utility_exit_with_message( "Can't determine filetype for file " + filename );
			}
		} // if/else extension recognized
	} // type = ""
	return parse_file( file, type, n_entries ); // file will be closed on destruction.
}

utility::vector1< MolFileIOMoleculeOP > MolFileIOReader::parse_file( std::istream & file, std::string type, core::Size n_entries /*=0*/ ) {
	ObjexxFCL::lowercase( type );
	utility::vector1< MolFileIOMoleculeOP > molecules;

	if ( type == "mol" || type == "sdf" ) {
		SDFParser parser;
		molecules = parser.parse( file, n_entries );
	} else if ( type == "mol2" ) {
		TR.Error << "Loading of mol2 files via this method not currently supported." << std::endl;
	} else if ( type == "pdb" ) {
		TR.Error << "Loading of pdb files via this method not currently supported." << std::endl;
	} else if ( type == "params" ) {
		TR.Error << "Loading of params files via this method not currently supported." << std::endl;
	} else {
		utility_exit_with_message( "Do not know how to handle molecule file of type " + type );
	}

	if ( ! molecules.size() ) {
		TR.Error << "Stream contained no recognized molecules!" << std::endl;
	}
	return molecules;
}


ResidueTypeOP convert_to_ResidueType( utility::vector1< MolFileIOMoleculeOP > molfile_data,
	std::string atom_type_tag,
	std::string elements_tag,
	std::string mm_atom_type_tag) {

	AtomTypeSetCOP atom_types( ChemicalManager::get_instance()->atom_type_set( atom_type_tag ) );
	ElementSetCOP elements( ChemicalManager::get_instance()->element_set( elements_tag ) );
	MMAtomTypeSetCOP mm_atom_types( ChemicalManager::get_instance()->mm_atom_type_set( mm_atom_type_tag ) );
	return convert_to_ResidueType(molfile_data, atom_types, elements, mm_atom_types);
}

ResidueTypeOP convert_to_ResidueType( utility::vector1< MolFileIOMoleculeOP > molfile_data,
	AtomTypeSetCOP atom_types,
	ElementSetCOP element_type_set,
	MMAtomTypeSetCOP mm_atom_types) {

	if ( molfile_data.size() == 0 ) {
		// This indicates a bad call, rather than bad data - early error
		utility_exit_with_message("ERROR: Cannot convert an empty vector of molecules to a ResidueType.");
	}

	std::map< core::chemical::sdf::AtomIndex, std::string > index_name_map;
	ResidueTypeOP restype = molfile_data[1]->convert_to_ResidueType(index_name_map, atom_types, element_type_set, mm_atom_types);
	if ( ! restype ) {
		TR.Info << "Could not load molecule '" << molfile_data[1]->name() << "' as a residue type." << std::endl;
		return ResidueTypeOP(nullptr);
	}

	if ( molfile_data.size() > 1 ) {
		rotamers::StoredRotamerLibrarySpecificationOP rotlib( new rotamers::StoredRotamerLibrarySpecification );
		for ( core::Size ii(1); ii <= molfile_data.size(); ++ii ) {
			std::map< std::string, core::Vector > location_map;
			for ( auto const & map_elem : index_name_map ) {
				MolFileIOAtomOP atom = molfile_data[ ii ]->atom_index( map_elem.first );
				if ( atom ) {
					location_map[ map_elem.second ] = atom->position();
				}
			}
			rotlib->add_rotamer( location_map );
		}
		restype->rotamer_library_specification( rotlib );
	}

	return restype;
}

utility::vector1< ResidueTypeOP >
convert_to_ResidueTypes( utility::vector1< MolFileIOMoleculeOP > molfile_data,
	bool load_rotamers,
	std::string atom_type_tag,
	std::string elements_tag,
	std::string mm_atom_type_tag) {

	AtomTypeSetCOP atom_types( ChemicalManager::get_instance()->atom_type_set( atom_type_tag ) );
	ElementSetCOP elements( ChemicalManager::get_instance()->element_set( elements_tag ) );
	MMAtomTypeSetCOP mm_atom_types( ChemicalManager::get_instance()->mm_atom_type_set( mm_atom_type_tag ) );
	return convert_to_ResidueTypes(molfile_data, load_rotamers, atom_types, elements, mm_atom_types);
}

utility::vector1< ResidueTypeOP >
convert_to_ResidueTypes( utility::vector1< MolFileIOMoleculeOP > molfile_data,
	bool load_rotamers,
	AtomTypeSetCOP atom_types,
	ElementSetCOP element_types,
	MMAtomTypeSetCOP mm_atom_types) {

	std::map< std::string, core::Size > name_index_map;
	utility::vector1< utility::vector1< MolFileIOMoleculeOP > > separated_molecules;
	std::string previous_entry("");

	for ( core::Size ii(1); ii <= molfile_data.size(); ++ii ) {
		if ( ! load_rotamers ) {
			separated_molecules.resize( separated_molecules.size() + 1 );
			separated_molecules[ separated_molecules.size() ].push_back( molfile_data[ii] );
		} else if ( previous_entry != "" && molfile_data[ii]->name() == "" ) {
			// If we have an empty name, assume it's a rotamer of the previous item (but only if there is a previous one)
			separated_molecules[ name_index_map[ previous_entry ] ].push_back( molfile_data[ii] );
		} else if ( name_index_map.find( molfile_data[ii]->name() ) != name_index_map.end() ) {
			// Existing name
			separated_molecules[ name_index_map[ molfile_data[ii]->name() ] ].push_back( molfile_data[ii] );
			previous_entry =  molfile_data[ii]->name();
		} else {
			// Brand new name
			separated_molecules.resize( separated_molecules.size() + 1 );
			separated_molecules[ separated_molecules.size() ].push_back( molfile_data[ii] );
			name_index_map[ molfile_data[ii]->name() ] = separated_molecules.size();
			previous_entry =  molfile_data[ii]->name();
		}
	}

	utility::vector1< ResidueTypeOP > residue_types;
	for ( core::Size jj(1); jj <= separated_molecules.size(); ++jj ) {
		TR.Debug << "Converting " << separated_molecules[jj][1]->name() << std::endl;
		ResidueTypeOP restype;
		try {
			restype = convert_to_ResidueType(separated_molecules[jj], atom_types, element_types, mm_atom_types);
		} catch( utility::excn::Exception & e ) {
			TR << ">>>>>> Skipping " << separated_molecules[jj][1]->name() << " due to Error: " << e.msg() << std::endl;
			continue;
		} catch (...) {
			TR << ">>>>>> Skipping " << separated_molecules[jj][1]->name() << " due to unspecified error!" << std::endl;
			continue;
		}
		if ( restype ) {
			residue_types.push_back( restype );
		}
	}

	return residue_types;
}


}
}
}
