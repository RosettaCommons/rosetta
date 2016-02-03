// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// comments in header file!
/// @file   src/core/chemical/sdf/mmCIFParser.cc
/// @author Steven Combs


#include <core/chemical/mmCIF/mmCIFParser.hh>
#include <core/types.hh>
#include <string>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <basic/database/open.hh>
#include <core/chemical/sdf/MolFileIOData.hh>

#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <utility/file/FileName.hh>
#include <utility/io/izstream.hh>

//Utility functions
#include <utility/string_util.hh>


//external CIF includes
#include <cifparse/CifFile.h>
#include <cifparse/CifParserBase.h>

#include <ctime>


namespace core {
namespace chemical {
namespace mmCIF {

//Load up the tracer for this class
static THREAD_LOCAL basic::Tracer TR( "core.io.mmCIF.mmCIFParser" );


mmCIFParser::mmCIFParser() {
	bond_string_to_sdf_size_[ "SING"] = 1;
	bond_string_to_sdf_size_[ "DOUB"] = 2;
	bond_string_to_sdf_size_[ "TRIP"] = 3;
}

sdf::MolFileIOMoleculeOP
mmCIFParser::parse(  std::string const &lines, std::string const &pdb_id){
	std::string diagnostics; //output from the parser about errors, etc.
	sdf::MolFileIOMoleculeOP molecule( new sdf::MolFileIOMolecule());
	CifFileOP cifFile( new CifFile);
	{
		CifParserOP cifParser( new CifParser(cifFile.get()) );
		cifParser->ParseString( lines, diagnostics);
		if ( !diagnostics.empty() ) {
			TR.Error << diagnostics << std::endl;
		} else {
			std::vector< std::string> blocks;
			cifFile->GetBlockNames( blocks);
			//You have to add the "#" string to the pdb_id. This is because of how
			//the cifFile parser interpets the lines. For whatever reason, it adds
			//a # to the data_??? block, where ??? is the 3 letter code for the
			//file
			if ( !cifFile->IsBlockPresent( pdb_id+"#") ) {
				TR.Error << pdb_id << " not found in components.cif" << std::endl;
			} else {
				Block& block = cifFile->GetBlock( pdb_id+"#");
				molecule = get_molfile_molecule( block);
				molecule->name( pdb_id );
				return molecule;
			}
		}
	}
	return molecule;
}

utility::vector1< sdf::MolFileIOMoleculeOP>
mmCIFParser::parse(std::string const &filename){
	utility::vector1< sdf::MolFileIOMoleculeOP> molecules;
	sdf::MolFileIOMoleculeOP molecule( new sdf::MolFileIOMolecule());
	std::string diagnostics; //output from the parser about errors, etc.
	CifFileOP cifFile( new CifFile);
	{
		CifParserOP cifParser( new CifParser(cifFile.get()) );
		cifParser->Parse( filename, diagnostics);
		//this assumes that the very first block being passed is the one being used.
		std::vector< std::string> block_names;
		cifFile->GetBlockNames( block_names);
		for ( unsigned int ii=0; ii< block_names.size(); ++ii ) {
			Block& block = cifFile->GetBlock( block_names[ ii]);
			molecule = get_molfile_molecule( block);
			molecule->name( block_names[ ii]);
			molecules.push_back( molecule);
		}
	}
	return molecules;
}

sdf::MolFileIOMoleculeOP
mmCIFParser::get_molfile_molecule(Block& block){
	sdf::MolFileIOMoleculeOP molecule( new sdf::MolFileIOMolecule());
	//only proceed if the tables for bonds and atoms are present
	if ( block.IsTablePresent("chem_comp_atom") ) {
		//get the atom and bond composition table
		ISTable& atom_comp = block.GetTable("chem_comp_atom");

		//this is to map the atom names to a core::Size value (id). This is added when adding
		//bonds to the MolIO object.
		std::map< std::string, core::Size> atom_name_to_id;


		//prefer the ideal coordinates, but if not found, use cartesian coordinates
		std::string xyz_start_type, xyz_end_type;
		if ( atom_comp.IsColumnPresent( "pdbx_model_Cartn_x_ideal") ) {
			xyz_start_type = "pdbx_model_Cartn_x_ideal";
			xyz_end_type = "pdbx_model_Cartn_z_ideal";
			if ( atom_comp( 0, xyz_start_type ) == "?" ) {
				// for some entries they're present but undefined - try the model coordinates instead
				xyz_start_type = "model_Cartn_x";
				xyz_end_type = "model_Cartn_z";
			}
		} else {
			xyz_start_type = "model_Cartn_x";
			xyz_end_type = "model_Cartn_z";
		}
		if ( atom_comp( 0, xyz_end_type ) == "?" ) {
			utility_exit_with_message( "No usable coordinates for mmCIF file for " + block.GetName() );
		}

		//store the atom_id_type we will be using. atom id should be the atom name
		std::string atom_name_type( atom_comp.IsColumnPresent( "atom_id") ? "atom_id" : "pdbx_component_atom_id");

		//start atom block
		for ( unsigned int ii = 0; ii< atom_comp.GetNumRows(); ++ii ) {
			sdf::MolFileIOAtomOP atom( new sdf::MolFileIOAtom());
			//atom id is whatever the number we are on +1, because we dont do 0 based index
			atom->index( utility::string2int( atom_comp(ii, "pdbx_ordinal")));

			//set atom name
			std::string atom_name( atom_comp(ii, atom_name_type));
			atom->name( atom_name);

			//set map to index
			atom_name_to_id[ atom_name] = core::Size (ii+1); //pluse 1 because 0 index

			//set element name
			atom->element( atom_comp(ii, "type_symbol"));

			//get the xyz cordinates
			std::vector< std::string> atom_coords;
			atom_comp.GetRow(atom_coords, ii, xyz_start_type, xyz_end_type);
			core::Real x = utility::string2float( atom_coords[ 0]);
			core::Real y = utility::string2float( atom_coords[ 1]);
			core::Real z = utility::string2float( atom_coords[ 2]);
			//set xyz coordinates
			atom->position( core::Vector(x, y, z));

			std::string charge( atom_comp(ii, "charge"));
			if ( charge == "?" ) {
				atom->formal_charge( 0 );
			} else {
				atom->formal_charge( utility::string2int( charge ) );
			}

			molecule->add_atom( atom);
		}

		if ( block.IsTablePresent( "chem_comp_bond" ) ) {
			ISTable& bond_comp = block.GetTable("chem_comp_bond");

			//start bond block
			for ( unsigned int ii=0; ii < bond_comp.GetNumRows(); ++ii ) {
				sdf::MolFileIOBondOP bond( new sdf::MolFileIOBond());

				std::string source( bond_comp(ii, "atom_id_1")); //atom 1
				std::string target( bond_comp(ii, "atom_id_2")); //atom 2 - I guess thats self explanatory

				bond->atom1( atom_name_to_id[ source]);
				bond->atom2( atom_name_to_id[ target]);

				core::Size bond_type( bond_string_to_sdf_size_[ bond_comp(ii, "value_order")]);
				bond->sdf_type( bond_type); //bond order

				molecule->add_bond( bond);
			}
		} else if ( atom_comp.GetNumRows() > 1 ) {
			TR.Error << "Cannot parse CIF file. No bond block (chem_comp_bond) found for multi-atom entry " << block.GetName() << std::endl;
		} // else one atom entry without bond block
	} else {
		TR.Error << "Cannot parse CIF file. No atom block (chem_comp_atom) found for " << block.GetName() << std::endl;
	}
	return molecule;
}


}
}
}
