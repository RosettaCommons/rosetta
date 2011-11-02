// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file core/chemical/sdf/mol_parser.cc
///
/// @brief Implementation of molfile parser
/// @author Sam DeLuca

#include <string>
#include <algorithm>
#include <utility/vector1.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/sdf/v2_parser.hh>
#include <core/chemical/sdf/v3_parser.hh>
#include <core/chemical/sdf/mol_parser.hh>
#include <basic/Tracer.hh>
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>
#include <fstream>

namespace core {
namespace chemical {
namespace sdf {

static basic::Tracer MolParserTracer("core.chemical.sdf.mol_parser");

MolFileParser::MolFileParser(const utility::vector1<std::string> mol_file_lines) : mol_file_lines_(mol_file_lines)
{
	MolParserTracer << "init mol file parser" << std::endl;
}

MolFileParser::MolFileParser(const std::string file_name)
{
	utility::io::izstream infile;
	utility::vector1<std::string> lines;
	infile.open(file_name.c_str());
	while(!infile.eof())
	{
		std::string line;
		getline(infile,line);
		lines.push_back(line);
	}
	infile.close();
	mol_file_lines_= lines;
}

core::chemical::ResidueType MolFileParser::GetResidueType()
{
	return *molecule_container_;
}

core::chemical::ResidueTypeOP MolFileParser::GetResidueTypeOP()
{
	MolParserTracer.Debug << "returning molecule name " << molecule_container_->name() <<std::endl;
	return molecule_container_;
}

void MolFileParser::parse_mol_file(core::chemical::AtomTypeSetCAP atom_types, core::chemical::ElementSetCAP elements, core::chemical::MMAtomTypeSetCAP mm_atom_types, core::chemical::orbitals::OrbitalTypeSetCAP orbital_type_set)
{

	MolParserTracer << "**WARNING**: this is extremely experimental, don't use for production purposes yet"<<std::endl;
	//molecule name
	molecule_name_ = mol_file_lines_[1];
	//molecule info
	molecule_info_ = mol_file_lines_[2];
	//molecule comments
	molecule_comments_ = mol_file_lines_[3];

	mol_data_.parse_mol_data(mol_file_lines_);
	//mol_data_.print();
	//make a new residue object, give it a name
	molecule_container_ = new core::chemical::ResidueType(atom_types, elements, mm_atom_types, orbital_type_set);
	molecule_container_->name(molecule_name_);
	molecule_container_->name3(molecule_comments_.substr(0,3));
	molecule_container_->name1(molecule_comments_[1]);

	molecule_container_->set_mol_data(mol_data_);
	molecule_container_->add_property("LIGAND");

	//pass the ctab and the new residue object to the ctab parser

	//first figure out if we have a V2000 or V3000 connection table
	MolParserTracer <<mol_file_lines_[4] <<std::endl;
	std::string ctab_version(mol_file_lines_[4].substr(34,5));
	MolParserTracer <<"ctab version is \"" << ctab_version <<"\"" << std::endl;
	if(ctab_version == "V2000")
	{
		MolParserTracer.Debug << "found a V2000 ctab" <<std::endl;
		//the rest is the ctab, copy it
		utility::vector1<std::string> connection_table;
		connection_table.resize(mol_file_lines_.size()-3);
		std::copy(mol_file_lines_.begin()+3,mol_file_lines_.end(),connection_table.begin());
		MolParserTracer.Debug << "connection table of size " << connection_table.size() << std::endl;
		V2Parser connection_table_parser(connection_table, molecule_container_,mol_data_);
		MolParserTracer.Debug << "Passing molecule " << molecule_container_->name() << " to ctab parser" <<std::endl;
		connection_table_parser.ParseTable();
		molecule_container_ = connection_table_parser.GetResidueType();
	}else if(ctab_version == "V3000")
	{
		//the rest is the ctab, copy it
		utility::vector1<std::string> connection_table;
		connection_table.resize(mol_file_lines_.size()-3);
		std::copy(mol_file_lines_.begin()+3,mol_file_lines_.end(),connection_table.begin());

		V3Parser connection_table_parser(connection_table, molecule_container_,mol_data_);
		connection_table_parser.ParseTable();
		molecule_container_ = connection_table_parser.GetResidueType();
	}

	//MolParserTracer <<"finished"<<std::endl;

	//fix the order and sort everything
	//molecule_container_->finalize();

	//set neighbor atom and internal coordinates



	molecule_container_->assign_internal_coordinates();

	//the ctab parser only sets default charges, charges need to be normalized to zero

	core::Real total_charge = 0.0;
	for(core::Size index =1; index <= molecule_container_->natoms(); ++index)
	{
		total_charge += molecule_container_->atomic_charge(index);
	}

	core::Real charge_offset = -total_charge/static_cast<core::Real>(molecule_container_->natoms());
	for(core::Size index = 1;index <= molecule_container_->natoms();++index)
	{
		std::string atom_name = molecule_container_->atom_name(index);
		core::Real starting_charge = molecule_container_->atomic_charge(index);
		molecule_container_->set_atomic_charge(atom_name,starting_charge+charge_offset);
	}
	molecule_container_->finalize();

	std::string preset_nbr = FindNbrAtom();
	if(preset_nbr != "")
	{

		//if you uncomment this tracer statement ChemicalManager segfaults and i have no idea why
		MolParserTracer.Debug << "assigning \"" <<preset_nbr<< "\" as nbr atom" <<std::endl;
		molecule_container_->nbr_atom(preset_nbr);
		//molecule_container_->assign_neighbor_atom();
	}else
	{
		molecule_container_->assign_neighbor_atom();
	}

}

std::string MolFileParser::GetMoleculeComments()
{
	return molecule_comments_;
}

std::string MolFileParser::GetMoleculeName()
{
	return molecule_name_;
}

std::string MolFileParser::GetMoleculeInfo()
{
	return molecule_info_;
}

std::string MolFileParser::FindNbrAtom()
{
	std::string nbr_atom_data = mol_data_.get_mol_data("Rosetta nbr_atom");
	if(nbr_atom_data == "")
	{
		return nbr_atom_data;
	}
	core::Size nbr_atom_index(utility::from_string<core::Size>(nbr_atom_data,core::Size(0) ) );
	std::string nbr_atom = molecule_container_->atom_name(nbr_atom_index);
	/*
	for(core::Size index=1;index <= mol_file_lines_.size();++index)
	{
		std::string line = mol_file_lines_[index];

		if(line.find("> <Rosetta nbr_atom>")!=std::string::npos)
		{
			return mol_file_lines_[index+1];
		}

	}
	*/
	return nbr_atom;
}



} // sdf
} // io
} // core
