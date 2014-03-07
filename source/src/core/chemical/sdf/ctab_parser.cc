// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   core/io/sdf/ctab_parser.cc
///
/// @brief
/// @author Sam DeLuca

#include <basic/Tracer.hh>
#include <core/chemical/sdf/ctab_parser.hh>
#include <utility/string_util.hh>
//#include <core/chemical/ResidueType.hh>
//#include <protocols/ligand_docking/ColoredGraph.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomType.hh>
#include <utility/exit.hh>
#include <numeric/xyzVector.hh>
#include <vector>
#include <algorithm>
#include <map>
#include <sstream>
#include <cstring>
#include <string>
#include <stdlib.h>

// Boost Headers
#include <boost/foreach.hpp>

namespace core {
namespace chemical {
namespace sdf {


static basic::Tracer ctabParserTracer("core.io.sdf.ctab_parser");
/*
elementToType::elementToType()
{
	e_to_t["C"]="CH3"; //Mult
	e_to_t["N"]="Nbb"; //Mult
	e_to_t["O"]="OH"; //Mult
	e_to_t["S"]="S";
	e_to_t["P"]="P";
	e_to_t["H"]="Hpol"; //Mult
	e_to_t["F"]="F";
	e_to_t["Cl"]="Cl";
	e_to_t["Br"]="Br";
	e_to_t["I"]="I";
	e_to_t["Zn"]="Zn2p";
	e_to_t["Fe"]="Fe2p"; //Mult
	e_to_t["Mg"]="Mg2p";
	e_to_t["Ca"]="Ca2p";
	e_to_t["Na"]="Na1p";
	e_to_t["K"]="K1p";

}

std::string elementToType::get(std::string key)
{
	std::map<std::string, std::string>::iterator val = e_to_t.find(key);
	std::string type;
	if (val==e_to_t.end()) {type="VIRT";}
	else { type=val->second;}
	return type;
}

*/
ctabV2000Parser::ctabV2000Parser(const utility::vector1<std::string> connection_table_lines, core::chemical::ResidueTypeOP molecule_container, MolData mol_data) :
		connection_table_lines_(connection_table_lines), molecule_container_(molecule_container), mol_data_(mol_data)
{

}

void ctabV2000Parser::ParseTable()
{
	std::string counts_line = connection_table_lines_[1];
	core::Size atom_count = atoi(counts_line.substr(0,3).c_str());
	core::Size bond_count = atoi(counts_line.substr(3,3).c_str());

	ctabParserTracer.Debug << atom_count << " atoms and " << bond_count << " bonds" <<std::endl;

	std::string version_tag = counts_line.substr(34,5);
	if(version_tag != "V2000")
	{
		utility_exit_with_message("This doesnt look like a V2000 CTAB, bailing out");
	}

	core::Size last_atom = atom_count;
	core::Size last_bond = last_atom+bond_count;

	for(core::Size line_number = 1; line_number < connection_table_lines_.size(); ++line_number)
	{
		if(line_number <= last_atom)
		{
			ParseAtom(connection_table_lines_[line_number+1],line_number);
		}else if(line_number <= last_bond)
		{
			ParseBond(connection_table_lines_[line_number+1]);
		}
	}

	//Iterate across the atoms, adding Hs and finding their types
	std::map<core::Size, std::string>::iterator atom_iterator;
	for(atom_iterator = index_to_names_map_.begin();
			atom_iterator != index_to_names_map_.end();
			++atom_iterator)
	{
		std::string atomname=atom_iterator->second;
		core::Size atomno = molecule_container_->atom_index(atomname);
		set_atom_type(atomno, atomname);
	}

	BOOST_FOREACH(addedH added, added_H_){
		index_to_names_map_.insert(
				std::pair<core::Size,std::string>(added.atom_number,"H"+added.atom_number)
		);
	}


}

//Sets the atom type and adds H if necessary.
void ctabV2000Parser::set_atom_type(core::Size atomno, std::string atomname)
{
	atomTyper typer = atomTyper(atomno, molecule_container_);
	if(typer.get_element()=="C"||typer.get_element()=="N"||typer.get_element()=="O")
	{
		core::SSize total_bonds=0;
		char ele=typer.get_element().at(0);
		switch(ele)
		{
		case 'C':
			total_bonds++;
		case 'N':
			total_bonds++;
		case 'O':
			total_bonds+=2;
		}
		total_bonds+=molecule_container_->atom(atomno).charge();
		total_bonds-=typer.getNumBonds();
		while(total_bonds>0)
		{
			total_bonds--;
			addedH newH;
			newH.atom_number=++current_atom_;
			newH.atom_type=(ele=='C')?"Hapo":"Hpol";
			newH.bonded_atom_name=atomname;
			added_H_.push_back(newH);
			molecule_container_->add_atom("H"+newH.atom_number,
					newH.atom_type,DEFAULT_MM_ATOM_TYPE_,0);
			molecule_container_->add_bond(newH.bonded_atom_name,
					std::string("H"+newH.atom_number),core::chemical::SingleBond);
			typer = atomTyper(atomno, molecule_container_);
		}
	}
	molecule_container_->set_atom_type(atomname, typer.getType());

	//set default charge
	core::chemical::ChemicalManager* chemical_manager = core::chemical::ChemicalManager::get_instance();
	core::chemical::AtomTypeSetCAP atom_type_set = chemical_manager->atom_type_set("fa_standard");
	core::Size atom_type_index = atom_type_set->atom_type_index(typer.getType());
	core::Size parameter_index = atom_type_set->extra_parameter_index("CHARGE");
	core::Real charge = atom_type_set->operator[](atom_type_index).extra_parameter(parameter_index);
	molecule_container_->atom( atomname ).charge(charge));
}

void ctabV2000Parser::ParseAtom(std::string atom_line, core::Size atom_number)
{
	core::Real x_coord = atof(atom_line.substr(0,10).c_str());
	core::Real y_coord = atof(atom_line.substr(10,10).c_str());
	core::Real z_coord = atof(atom_line.substr(20,10).c_str());
	ctabParserTracer.Debug << "atom " <<atom_number << " has coordinates " << x_coord << ',' << y_coord << ',' << z_coord << std::endl;
	std::string element_name = atom_line.substr(31,3).c_str();
	if(element_name.at(1)==' ') {
		element_name=element_name.substr(0,1);
	} else if (element_name.at(2)==' ') {
		element_name=element_name.substr(0,2);
	}

	std::string atom_number_string;
	std::stringstream convert_stream;
	convert_stream << atom_number;
	atom_number_string = convert_stream.str();

	//Set the atom type to a default based on the element.
	std::string atom_type = element_to_default_type.get(element_name);

	current_atom_=atom_number;

	std::string element_id = element_name+ atom_number_string;
	utility::add_spaces_left_align(element_id,4);

	core::Size charge = atoi(atom_line.substr(36,3).c_str());

	numeric::xyzVector<core::Real> coordinates(x_coord,y_coord,z_coord);

	index_to_names_map_.insert(std::pair<core::Size,std::string>(atom_number,element_id));
	molecule_container_->add_atom(element_id,atom_type,DEFAULT_MM_ATOM_TYPE_,charge);
	molecule_container_->set_xyz(element_id,coordinates);
}

void ctabV2000Parser::ParseBond(std::string bond_line)
{
	core::Size atom1_index = atoi(bond_line.substr(0,3).c_str());
	core::Size atom2_index = atoi(bond_line.substr(3,3).c_str());
	core::chemical::BondName bond_type = static_cast<core::chemical::BondName>(atoi(bond_line.substr(6,3).c_str()));

	std::string atom1_id(index_to_names_map_.find(atom1_index)->second);
	std::string atom2_id(index_to_names_map_.find(atom2_index)->second);
	ctabParserTracer.Debug << "bond " << atom1_id << " to " << atom2_id << " of type " << bond_type <<  std::endl;
	molecule_container_->add_bond(atom1_id,atom2_id,bond_type);
}

 std::map<core::Size, std::string>  ctabV2000Parser::ParseAtomTypeData()
 {
	 std::string atom_type_data = "";
	 std::map<core::Size, std::string> data_map;
	 for(core::Size index = 1; index <= connection_table_lines_.size(); ++index )
	 {
		 std::string line = connection_table_lines_[index];
		 if(line.find("> <Rosetta AtomTypes>")!= std::string::npos)
		 {
			 atom_type_data = connection_table_lines_[index+1];
			 break;
		 }
	 }
	 if(atom_type_data == "")
	 {
		 return data_map;
	 }
	 else
	 {
		 utility::vector1<std::string> tokens=utility::split(atom_type_data);
		 for(core::Size index = 1; index <= tokens.size();++index)
		 {
			 std::string current_token = tokens[index];
			 utility::vector1<std::string> token_split=utility::split(current_token);
			 utility::trim(token_split[2],"(");
			 utility::trim(token_split[3],")");

			 core::Size atomno = atoi(token_split[2].c_str());
			 std::pair<core::Size, std::string> atom_type_point(atomno,token_split[3]);
			 data_map.insert(atom_type_point);
		 }

	 }
	 return data_map;
 }


core::chemical::ResidueTypeOP ctabV2000Parser::GetResidueType()
{
	return molecule_container_;
}


ctabV3000Parser::ctabV3000Parser(const utility::vector1<std::string> connection_table_lines, core::chemical::ResidueTypeOP molecule_container, MolData mol_data ) :
		connection_table_lines_(connection_table_lines), molecule_container_(molecule_container), mol_data_(mol_data)
{ }

core::chemical::ResidueTypeOP ctabV3000Parser::GetResidueType()
{
	return molecule_container_;

}

void ctabV3000Parser::set_atom_type(core::Size atomno, std::string atomname)
{
	atomTyper typer = atomTyper(atomno, molecule_container_);
	if(typer.get_element()=="C"||typer.get_element()=="N"||typer.get_element()=="O")
	{
		core::Size total_bonds=0;
		char ele=typer.get_element().at(0);
		switch(ele)
		{
		case 'C':
			total_bonds++;
		case 'N':
			total_bonds++;
		case 'O':
			total_bonds+=2;
		}
		total_bonds+=molecule_container_->atom(atomno).charge();
		total_bonds-=typer.getNumBonds();
		while(total_bonds>0)
		{
			total_bonds--;
			addedH newH;
			newH.atom_number=++current_atom_;
			newH.atom_type=(ele=='C')?"Hapo":"Hpol";
			newH.bonded_atom_name=atomname;
			added_H_.push_back(newH);
			molecule_container_->add_atom("H"+newH.atom_number,
					newH.atom_type,DEFAULT_MM_ATOM_TYPE_,0);
			molecule_container_->add_bond(newH.bonded_atom_name,
					std::string("H"+newH.atom_number),core::chemical::SingleBond);
			typer = atomTyper(atomno, molecule_container_);
		}
	}
	molecule_container_->set_atom_type(atomname, typer.getType());
}

core::Real ctabV3000Parser::FindExtraParameter(std::vector<std::string> extra_parameters,const std::string query )
{
	BOOST_FOREACH(std::string extra_parameter, extra_parameters){
		if(extra_parameter.find(query) == std::string::npos)
		{
			continue;
		}
		else
		{
			core::Size value_length = extra_parameter.size() - (query.size()+1);
			core::Size value_start = query.size()+1;

			std::string value_string = extra_parameter.substr(value_start,value_length);
			core::Real value = atof(value_string.c_str());
			return value;
		}
	}
	return 0.0;
}
void ctabV3000Parser::ParseTable()
{

	bool atom_block(false);
	bool bond_block(false);

	std::map<core::Size,std::string> atom_type_data = this->ParseAtomTypeData();

	BOOST_FOREACH(std::string current_line, connector_table_lines_){

		utility::vector1<std::string> line_vector(utility::split(current_line));
		//if(line_vector[1] != "M" || line_vector[2] != "V30")
		//{
		//	utility_exit_with_message("This doesn't look like a V3000 CTAB, bailing out");
		//}

		core::Size atom_count(0);
		core::Size bond_count(0);
		if(line_vector[2] == "END")
		{
			break;
		}
		if(line_vector[2] != "V30" || line_vector.size() <2)
		{
			continue;
		}
		if(line_vector[3] == "COUNTS")
		{
			atom_count = atoi(line_vector[4].c_str());
			bond_count = atoi(line_vector[5].c_str());
		}else if(line_vector[3] == "BEGIN")
		{
			if (line_vector[4] == "ATOM")
				atom_block = true;
			else if (line_vector[4] == "BOND")
				bond_block = true;
		}else if(line_vector[3] == "END")
		{
			if (line_vector[4] == "ATOM")
				atom_block = false;
			else if (line_vector[4] == "BOND")
				bond_block = false;
		}else
		{
			if (atom_block)
				ParseAtom(current_line);
			else if (bond_block)
				ParseBond(current_line);
		}
	}

	std::map<core::Size, std::string>::iterator atom_iterator;
	for(atom_iterator = index_to_names_map_.begin();
			atom_iterator != index_to_names_map_.end();
			++atom_iterator)
	{
		std::string atomname= atom_iterator->second;
		core::Size atomno = molecule_container_->atom_index(atomname);
		std::map<core::Size,std::string>::iterator atom_type_it = atom_type_data.find(atomno);
		if(atom_type_it != atom_type_data.end())
		{
			std::string atom_type = atom_type_it->second;
			molecule_container_->set_atom_type(atomname,atom_type);
		} else
		{
			set_atom_type(atomno,atomname);
		}
	}

	BOOST_FOREACH(addedH added, added_H_){
		index_to_names_map_.insert(std::pair<core::Size,std::string>(
				added.atom_number,"H"+added.atom_number));
	}

}

void ctabV3000Parser::ParseAtom(const std::string atom_line)
{
	utility::vector1<std::string> atom_vector(utility::split(atom_line));

	core::Size index = atoi(atom_vector[3].c_str());
	std::string element_name(atom_vector[4]);
	std::string element_id(element_name+atom_vector[3]);
	utility::add_spaces_left_align(element_id,4);
	core::Real x_coord = atof(atom_vector[5].c_str());
	core::Real y_coord = atof(atom_vector[6].c_str());
	core::Real z_coord = atof(atom_vector[7].c_str());

	numeric::xyzVector<core::Real> coordinates(x_coord,y_coord,z_coord);



	//std::vector<std::string> extra_parameters;
	//std::copy(atom_vector.begin()+8, atom_vector.end(),extra_parameters.begin());
	core::Real charge(FindExtraParameter(atom_vector,"CHG"));

	index_to_names_map_.insert(std::pair<core::Size,std::string>(index,element_id));
	//Set the atom type to a default based on the element.
	std::string atom_type = element_to_default_type.get(element_name);
	molecule_container_->add_atom(element_id,atom_type,"X",charge);
	molecule_container_->set_xyz(element_id,coordinates);

}

void ctabV3000Parser::ParseBond(const std::string bond_line)
{
	utility::vector1<std::string> bond_vector(utility::split(bond_line));

	core::Size index = atoi(bond_vector[3].c_str());
	core::chemical::BondName type = static_cast<core::chemical::BondName>(atoi(bond_vector[4].c_str()));

	core::Size atom1_index = atoi(bond_vector[5].c_str());
	core::Size atom2_index = atoi(bond_vector[6].c_str());

	std::string atom1_id(index_to_names_map_.find(atom1_index)->second);
	std::string atom2_id(index_to_names_map_.find(atom2_index)->second);

	std::vector<std::string> extra_parameters;
	std::copy(bond_vector.begin()+8, bond_vector.end(),extra_parameters.begin());

	molecule_container_->add_bond(atom1_id,atom2_id,type);
}

std::map<core::Size, std::string>  ctabV3000Parser::ParseAtomTypeData()
{
	 utility::vector1<std::string> tokens = mol_data_.get_mol_data_string_vector("Rosetta AtomTypes",' ');
	 std::map<core::Size, std::string> data_map;
	 /*
	 for(core::Size index = 1; index <= connection_table_lines_.size(); ++index )
	 {
		 std::string line = connection_table_lines_[index];
		 if(line.find("> <Rosetta AtomTypes>")!= std::string::npos)
		 {
			 atom_type_data = connection_table_lines_[index+1];
			 break;
		 }
	 }
	 */

	 if(tokens.size() == 0)
	 {
		 return data_map;
	 }
	 else
	 {
		 //utility::vector1<std::string> tokens=utility::split(atom_type_data);
		 for(core::Size index = 0; index < tokens.size();++index)
		 {
			 std::string current_token = tokens[index];
			 utility::vector1<std::string> token_split=utility::string_split(current_token,',');
			 utility::trim(token_split[1],"(");
			 utility::trim(token_split[2],")");

			 core::Size atomno = atoi(token_split[1].c_str());
			 std::pair<core::Size, std::string> atom_type_point(atomno,token_split[2]);
			 data_map.insert(atom_type_point);
		 }

	 }
	 return data_map;
}

} // sdf
} // io
} // core
