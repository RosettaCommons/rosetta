// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/chemical/sdf/v3_parser.cc
/// @author Sam DeLuca

#include <core/chemical/sdf/v3_parser.hh>

#include <numeric/xyzVector.hh>
#include <core/chemical/ResidueType.hh>

#include <utility/string_util.hh>

// Boost Headers
#include <boost/foreach.hpp>

#include <utility/vector1.hh>

namespace core {
namespace chemical {
namespace sdf {

V3Parser::V3Parser(utility::vector1<std::string> const & connection_table_lines,
		core::chemical::ResidueTypeOP molecule_container, MolData const & mol_data):
		CtabBase(connection_table_lines,molecule_container,mol_data)
{

}

void V3Parser::ParseTable()
{
	core::Size ctab_length = this->connection_table_length();

	bool atom_block(false);
	bool bond_block(false);
	//std::map<core::Size,std::string> atom_type_data = this->ParseAtomTypeData();
	//core::Size atom_count(0);
	//core::Size bond_count(0);
	for(core::Size line_number = 1; line_number <= ctab_length; ++line_number)
	{
		std::string current_line(this->connection_table_line(line_number));
		utility::vector1<std::string> line_vector(utility::split(current_line));

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
			//atom_count = atoi(line_vector[4].c_str());  // set but never used ~Labonte
			//bond_count = atoi(line_vector[5].c_str());  // set but never used ~Labonte
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
				ParseAtom(current_line,0 /*dummy arg*/);
			else if (bond_block)
				ParseBond(current_line);
		}
	}

	//Iterate across the atoms, finding their types
	this->fix_atom_types();
	/*
	core::chemical::ResidueTypeOP residue(this->GetResidueType());
	for(core::Size atom_index = 1; atom_index <=atom_count; ++atom_index)
	{
		std::string atom_name(this->atom_name_from_index(atom_index));
		core::Size atomno = residue->atom_index(atom_name);
		this->set_atom_type(atomno,atom_name);
	}
	*/
}

void V3Parser::ParseAtom(std::string const atom_line, core::Size const )
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

	core::Real charge(FindExtraParameter(atom_vector,"CHG"));

	this->add_index_name_pair(index,element_id);

	core::chemical::ResidueTypeOP molecule_container = this->GetResidueType();
	std::string atom_type(element_to_default_type().get(element_name));
	molecule_container->add_atom(element_id,atom_type,DEFAULT_MM_ATOM_TYPE_,charge);
	molecule_container->set_ideal_xyz(element_id,coordinates);
	//molecule_container->finalize();

}

void V3Parser::ParseBond(std::string const bond_line)
{
	utility::vector1<std::string> bond_vector(utility::split(bond_line));

	//core::Size index = atoi(bond_vector[2].c_str());
	core::chemical::BondName type = static_cast<core::chemical::BondName>(atoi(bond_vector[4].c_str()));

	core::Size atom1_index = atoi(bond_vector[5].c_str());
	core::Size atom2_index = atoi(bond_vector[6].c_str());
	if(this->check_for_aromatic(atom1_index,atom2_index))
	{
		type = core::chemical::AromaticBond;
	}

	std::string atom1_id(this->atom_name_from_index(atom1_index));
	std::string atom2_id(this->atom_name_from_index(atom2_index));

	core::chemical::ResidueTypeOP molecule_container = this->GetResidueType();
	molecule_container->add_bond(atom1_id,atom2_id,type);

}

core::Real V3Parser::FindExtraParameter(utility::vector1<std::string> const extra_parameters, std::string const query )
{
	BOOST_FOREACH(std::string current_parameter, extra_parameters){
		if(current_parameter.find(query) == std::string::npos)
		{
			continue;
		}
		else
		{
			core::Size value_length = current_parameter.size() - (query.size()+1);
			core::Size value_start = query.size()+1;

			std::string value_string = current_parameter.substr(value_start,value_length);
			core::Real value = atof(value_string.c_str());
			return value;
		}
	}
	return 0.0;
}

}
}
}
