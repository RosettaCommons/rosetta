// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/chemical/sdf/ctab_base.cc
/// @author Sam DeLuca

#include <core/chemical/sdf/ctab_base.hh>
#include <core/chemical/sdf/ctab_typer.hh>
#include <core/chemical/Atom.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/AtomType.hh>
// AUTO-REMOVED #include <utility/string_util.hh>

#include <basic/Tracer.hh>

#include <core/chemical/sdf/mol_util.hh>
#include <utility/vector1.hh>


namespace core {
namespace chemical {
namespace sdf {

static basic::Tracer ctabParserTracer("core.io.sdf.ctab_parser");

elementToType & element_to_default_type()
{
	static elementToType element_to_default_type_;
	return element_to_default_type_;
}


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

CtabBase::CtabBase(utility::vector1<std::string> const & connection_table_lines,
		core::chemical::ResidueTypeOP molecule_container, MolData const & mol_data) :
		connection_table_lines_(connection_table_lines),
		molecule_container_(molecule_container),
		mol_data_(mol_data)
{
	atom_type_data_map_ = parse_atom_type_data(mol_data_.get_mol_data("Rosetta AtomTypes"));
	bond_type_data_set_ = parse_bond_type_data(mol_data_.get_mol_data("PUBCHEM_BONDANNOTATIONS"));
}

CtabBase::~CtabBase()
{

}

core::chemical::ResidueTypeOP CtabBase::GetResidueType()
{
	return molecule_container_;
}

core::Size CtabBase::connection_table_length() const
{
	return connection_table_lines_.size();
}

std::string CtabBase::connection_table_line(core::Size const line_number) const
{
	return connection_table_lines_[line_number];
}

void CtabBase::add_index_name_pair(core::Size const index, std::string const atomname)
{
	index_to_names_map_.insert(std::pair<core::Size,std::string>(index,atomname));
}

std::string CtabBase::atom_name_from_index(core::Size const index) const
{
	std::map<core::Size,std::string>::const_iterator atom_map_it = index_to_names_map_.find(index);
	if(atom_map_it == index_to_names_map_.end())
	{
		return "";

	}else
	{
		return atom_map_it->second;
	}
}

/*
std::map<core::Size, std::string> CtabBase::ParseAtomTypeData()
{
	 utility::vector1<std::string> tokens = mol_data_.get_mol_data_string_vector("Rosetta AtomTypes",' ');
	 std::map<core::Size, std::string> data_map;

	 if(tokens.size() == 0)
	 {
		 return data_map;
	 }
	 else
	 {
		 //utility::vector1<std::string> tokens=utility::split(atom_type_data);
		 for(core::Size index = 1; index <= tokens.size();++index)
		 {

			 std::string current_token = tokens[index];
			 if(current_token.size() <=1)
			 {
				 continue;
			 }
			 utility::vector1<std::string> token_split=utility::string_split(current_token,',');
			 utility::trim(token_split[0],"(");
			 utility::trim(token_split[1],")");
			 //std::cout << current_token<<std::endl;
			 core::Size atomno = atoi(token_split[0].c_str());
			 std::pair<core::Size, std::string> atom_type_point(atomno,token_split[1]);
			 data_map.insert(atom_type_point);
		 }

	 }
	 return data_map;
}

*/

bool CtabBase::check_for_aromatic(core::Size lower, core::Size upper)
{
	if(bond_type_data_set_.find(BondData(lower,upper,8)) != bond_type_data_set_.end())
	{
		return true;
	}else
	{
		return false;
	}
}

void CtabBase::set_atom_type(core::Size const atomno, std::string const atomname)
{
	atomTyper typer = atomTyper(atomno, molecule_container_);
	//%TODO fix hydrogen adding code eventually
	/*
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
	*/

	//set default charge
	core::chemical::ChemicalManager* chemical_manager = core::chemical::ChemicalManager::get_instance();
	core::chemical::AtomTypeSetCAP atom_type_set = chemical_manager->atom_type_set("fa_standard");
	core::Size atom_type_index = atom_type_set->atom_type_index(typer.getType());
	core::Size parameter_index = atom_type_set->extra_parameter_index("CHARGE");
	core::Real charge = atom_type_set->operator[](atom_type_index).extra_parameter(parameter_index);
	molecule_container_->atom( atomname ).charge(charge);
	molecule_container_->set_atom_type(atomname, typer.getType());
}

void CtabBase::fix_atom_types()
{

	std::map<core::Size,std::string>::iterator atom_it;
	for(atom_it = index_to_names_map_.begin(); atom_it != index_to_names_map_.end();++atom_it)
	{
		std::string atom_name = atom_it->second;
		core::Size atom_no = atom_it->first;
		std::map<core::Size, std::string>::iterator atom_data_it = atom_type_data_map_.find(atom_no);
		if(atom_data_it == atom_type_data_map_.end())
		{
			set_atom_type(atom_no,atom_name);
		}else
		{
			//std::cout << atom_name << " " << atom_data_it->second <<std::endl;
			molecule_container_->set_atom_type(atom_name,atom_data_it->second);
		}

	}
}

/*
void CtabBase::fix_bond_types()
{
	 utility::vector1<std::string> tokens = mol_data_.get_mol_data_string_vector("PUBCHEM_BONDANNOTATIONS",'\n');

	 utility::vector1<std::string>::iterator token_it;
	 for(token_it = tokens.begin();token_it != tokens.end();++token_it)
	 {
		 utility::vector1<std::string> current_token(utility::split(*token_it));
		 if(current_token.size() < 3)
		 {
			 continue;
		 }

		 core::Size lower_id = utility::from_string(current_token[1],core::Size(0));
		 core::Size upper_id = utility::from_string(current_token[2],core::Size(0));
		 core::Size type = utility::from_string(current_token[3],core::Size(0));

		 //Currently we only support pubchem bondtype 8 (aromatic)
		 if(type == 8)
		 {

		 }
	 }
}
*/
}
}
}
