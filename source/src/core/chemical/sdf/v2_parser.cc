// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   /src/core/chemical/sdf/v2_parser.cc
/// @author Sam DeLuca

#include <core/chemical/sdf/v2_parser.hh>
#include <core/chemical/ResidueType.hh>
#include <utility/string_util.hh>

#include <utility/vector1.hh>


namespace core {
namespace chemical {
namespace sdf {

V2Parser::V2Parser(utility::vector1<std::string> const & connection_table_lines,
		core::chemical::ResidueTypeOP molecule_container, MolData const & mol_data):
		CtabBase(connection_table_lines,molecule_container,mol_data)
{ }

void V2Parser::ParseTable()
{
	std::string counts_line = this->connection_table_line(1);
	core::Size atom_count = atoi(counts_line.substr(0,3).c_str());
	core::Size bond_count = atoi(counts_line.substr(3,3).c_str());

	//ctabParserTracer.Debug << atom_count << " atoms and " << bond_count << " bonds" <<std::endl;
	//std::map<core::Size, std::string> atom_type_data = this->ParseAtomTypeData();
	std::string version_tag = counts_line.substr(34,5);

	assert(version_tag == "V2000");

	core::Size last_atom = atom_count;
	core::Size last_bond = last_atom+bond_count;

	for(core::Size line_number = 1; line_number < this->connection_table_length(); ++line_number)
	{
		if(line_number <= last_atom)
		{
			ParseAtom(this->connection_table_line(line_number+1),line_number);
		}else if(line_number <= last_bond)
		{
			ParseBond(this->connection_table_line(line_number+1));
		}
	}

	//Iterate across the atoms, finding their types
	this->fix_atom_types();
	/*
	core::chemical::ResidueTypeOP residue(this->GetResidueType());
	for(core::Size atom_index = 1; atom_index <= atom_count; ++atom_index)
	{
		std::string atom_name(this->atom_name_from_index(atom_index));
		core::Size atomno = residue->atom_index(atom_name);
		this->set_atom_type(atomno,atom_name);
	}
	*/
}

void V2Parser::ParseAtom(std::string const atom_line, core::Size const atom_number)
{
	core::Real x_coord = atof(atom_line.substr(0,10).c_str());
	core::Real y_coord = atof(atom_line.substr(10,10).c_str());
	core::Real z_coord = atof(atom_line.substr(20,10).c_str());
	//ctabParserTracer.Debug << "atom " <<atom_number << " has coordinates " << x_coord << ',' << y_coord << ',' << z_coord << std::endl;
	std::string element_name = atom_line.substr(31,3).c_str();
	if(element_name.at(1)==' ') {
		element_name=element_name.substr(0,1);
	} else if (element_name.at(2)==' ') {
		element_name=element_name.substr(0,2);
	}

	std::string atom_number_string=	utility::to_string<core::Size>(atom_number);

	//atom_number_string = convert_stream.str();

	//Set the atom type to a default based on the element.
	std::string atom_type = element_to_default_type().get(element_name);

	std::string element_id = element_name+ atom_number_string;
	utility::add_spaces_left_align(element_id,4);

	core::Size charge = atoi(atom_line.substr(36,3).c_str());
	numeric::xyzVector<core::Real> coordinates(x_coord,y_coord,z_coord);
	this->add_index_name_pair(atom_number,element_id);

	core::chemical::ResidueTypeOP residue(this->GetResidueType());
	residue->add_atom(element_id,atom_type,DEFAULT_MM_ATOM_TYPE_,charge);
	residue->set_ideal_xyz(element_id,coordinates);
	//residue->finalize();
}

void V2Parser::ParseBond(std::string const bond_line)
{
	core::Size atom1_index = atoi(bond_line.substr(0,3).c_str());
	core::Size atom2_index = atoi(bond_line.substr(3,3).c_str());
	core::chemical::BondName bond_type = static_cast<core::chemical::BondName>(atoi(bond_line.substr(6,3).c_str()));
	if(this->check_for_aromatic(atom1_index,atom2_index))
	{
		bond_type = core::chemical::AromaticBond;
	}
	std::string atom1_id(this->atom_name_from_index(atom1_index));
	std::string atom2_id(this->atom_name_from_index(atom2_index));

	core::chemical::ResidueTypeOP residue(this->GetResidueType());
	residue->add_bond(atom1_id,atom2_id,bond_type);
	//residue->finalize();

}

}
}
}
