// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file core/chemical/sdf/sdf_parser.cc
///
/// @brief sdf file parser implementation
/// @author Sam DeLuca

#include <core/chemical/sdf/sdf_parser.hh>
#include <core/chemical/sdf/mol_parser.hh>

#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <basic/Tracer.hh>

#include <string>
#include <map>

// Boost Headers
#include <boost/foreach.hpp>

namespace core {
namespace chemical {
namespace sdf {

static basic::Tracer SDFParserTracer("core.io.sdf.sdf_parser");

SDFParser::SDFParser(utility::vector1<std::string> file_vector) : file_vector_(file_vector)
{ }

void SDFParser::GenerateResidues()
{
	SDFParserTracer << "generating residues " <<std::endl;

	core::chemical::AtomTypeSetCAP atom_types = core::chemical::ChemicalManager::get_instance()->atom_type_set("fa_standard");
	core::chemical::MMAtomTypeSetCAP mm_atom_types = core::chemical::ChemicalManager::get_instance()->mm_atom_type_set("fa_standard");

	SDFParserTracer << "mol_file_map_ contains " << mol_file_map_.size() << " entries" <<std::endl;

	std::map<std::string,utility::vector1<std::string> >::iterator mol_file_map_iterator;
	for(mol_file_map_iterator = mol_file_map_.begin(); mol_file_map_iterator != mol_file_map_.end(); ++mol_file_map_iterator)
	{
		utility::vector1<std::string> current_mol_file( mol_file_map_iterator->second);
		MolFileParser mol_file_converter(current_mol_file);
		mol_file_converter.parse_mol_file(atom_types, mm_atom_types);
		core::chemical::ResidueTypeOP current_molecule(mol_file_converter.GetResidueTypeOP());
		SDFParserTracer << "inserting " << mol_file_map_iterator->first <<std::endl;
		molecule_map_.insert(std::pair<std::string, core::chemical::ResidueTypeOP>(mol_file_map_iterator->first,current_molecule));
		SDFParserTracer << "size after insertion: " << molecule_map_.size() <<std::endl;
	}
}

core::Size SDFParser::GetNumberOfResidues()
{
	return mol_file_map_.size();
}

utility::vector1<std::string> SDFParser::GetResidueNameVector()
{
	return molecule_name_vector_;
}

void SDFParser::SplitSDF()
{
	std::string current_name;
	utility::vector1<std::string> current_mol_block;
	utility::vector1<std::string> current_data_block;
	core::Size line_counter = 1;
	bool data_block = false;

	BOOST_FOREACH(std::string current_line, file_vector_){

		if(line_counter == 1)
		{
			//if the line counter is 1, we have a name line
			current_name = current_line;
			molecule_name_vector_.push_back(current_name);
		}

		if(current_line == "$$$$")
		{
			//end of struct, clear add the block vectors to the map, reset the counter to 1
			line_counter = 1;
			mol_file_map_.insert(std::pair<std::string, utility::vector1<std::string> >(current_name,current_mol_block));
			extended_data_map_.insert(std::pair<std::string,utility::vector1<std::string> >(current_name,current_data_block));

			current_mol_block.clear();
			current_data_block.clear();
			continue;
		}

		if(data_block)
		{
			current_data_block.push_back(current_line);
		}else
		{
			current_mol_block.push_back(current_line);
		}
		line_counter++;

	}
}

core::chemical::ResidueTypeOP SDFParser::GetResidueByName(std::string name)
{
	return molecule_map_.find(name)->second;
}



}
}
}
