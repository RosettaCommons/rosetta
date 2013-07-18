// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/chemical/sdf/mol_writer.cc
/// @author Sam DeLuca
/// @details this class outputs a residue in the form of a V3000 molfile, for details of the spec, see http://www.symyx.com/downloads/public/ctfile/ctfile.pdf
#include <core/chemical/sdf/mol_writer.hh>
#include <core/chemical/sdf/mol_util.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ResidueType.hh>
#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
#include <set>

// Boost Headers
#include <boost/foreach.hpp>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>

#define foreach BOOST_FOREACH

namespace core {
namespace chemical {
namespace sdf {



MolWriter::MolWriter() : line_header_("M  V30 ")
{}


void MolWriter::output_residue(utility::io::ozstream & output_stream, core::conformation::ResidueCOP residue)
{
	std::list<std::string> prepared_lines;
	std::list<std::string> metadata = this->compose_metadata(residue);
	std::list<std::string> ctab = this->compose_ctab(residue);
	std::list<std::string> typeinfo = this->compose_typeinfo(residue);
	std::list<std::string> nbr_atom = this->compose_nbr_atom(residue);
	std::list<std::string> job_data = this->compose_job_info();

	prepared_lines.insert(prepared_lines.end(),metadata.begin(),metadata.end());
	prepared_lines.insert(prepared_lines.end(),ctab.begin(),ctab.end());
	prepared_lines.insert(prepared_lines.end(),typeinfo.begin(),typeinfo.end());
	prepared_lines.insert(prepared_lines.end(),nbr_atom.begin(),nbr_atom.end());
	prepared_lines.insert(prepared_lines.end(),job_data.begin(),job_data.end());


	foreach(std::string line, prepared_lines){
		output_stream << line;
	}

}





void MolWriter::output_residue(utility::io::ozstream & output_stream, core::chemical::ResidueTypeOP residue_type)
{
	core::conformation::ResidueCOP residue_ptr(new  core::conformation::Residue(*residue_type,false) );
	//core::conformation::ResidueCOP residue_ptr= &residue;
	//std::cout <<residue_ptr->name3() <<std::endl;
	output_residue(output_stream,residue_ptr);
}


void MolWriter::output_residue(std::string file_name,core::conformation::ResidueCOP residue)
{
	utility::io::ozstream outfile;

	outfile.open(file_name.c_str(),std::ios::out | std::ios::binary);
	if(!outfile)
	{
		throw utility::excn::EXCN_FileNotFound("Cannot open file"+file_name);
	}
	output_residue(outfile,residue);
	outfile.close();
}


void MolWriter::output_residue(std::string file_name, core::chemical::ResidueTypeOP residue_type)
{
	utility::io::ozstream outfile;
	outfile.open(file_name.c_str(), std::ios::out | std::ios::binary);
	if(!outfile)
	{
		throw utility::excn::EXCN_FileNotFound("Cannot open file"+file_name);
	}
	output_residue(outfile,residue_type);
	outfile.close();
}

std::list<std::string> MolWriter::compose_metadata(core::conformation::ResidueCOP residue)
{
	std::list<std::string> lines;

	std::string const name_line = residue->name()+"\n";
	//for specification see page 35 of the sybl molfile specification document
	std::string const info_line = "  Rosetta           3D                              \n";
	std::string const name3_line = residue->name3()+"\n";
	std::string const counts_line = "  0  0  0     0  0            999 V3000\n";

	lines.push_back(name_line);
	lines.push_back(info_line);
	lines.push_back(name3_line);
	lines.push_back(counts_line);
	return lines;
}

std::list<std::string> MolWriter::compose_ctab(core::conformation::ResidueCOP residue)
{
	std::list<std::string> lines;
	std::string begin_header = line_header_+"BEGIN CTAB\n";
	std::string end_header = line_header_+"END CTAB\n";

	core::Size n_atoms = residue->natoms();
	core::Size n_bonds = 0;
	for(core::Size index = 1; index <= n_atoms;++index)
	{
		n_bonds += residue->bonded_neighbor(index).size();
	}
	//bonds are seen twice in the bond matrix (a-b and b-a)
	n_bonds = n_bonds/2;

	std::string counts = line_header_+"COUNTS "+ utility::to_string<core::Size>(n_atoms)+" "+
						utility::to_string<core::Size>(n_bonds)+" " + "0" + " " + "0" + " "+ "0"+"\n";

	//compose bonds and atoms, append all this to the ctab

	std::list<std::string> atom_lines = this->compose_atoms(residue);
	std::list<std::string> bond_lines = this->compose_bonds(residue);

	lines.push_back(begin_header);
	lines.push_back(counts);
	lines.insert(lines.end(),atom_lines.begin(),atom_lines.end());
	lines.insert(lines.end(),bond_lines.begin(),bond_lines.end());
	lines.push_back(end_header);
	lines.push_back("M  END\n");

	return lines;
}

std::list<std::string> MolWriter::compose_atoms(core::conformation::ResidueCOP residue)
{
	std::list<std::string> lines;
	std::string begin_header = line_header_+"BEGIN ATOM"+"\n";
	std::string end_header = line_header_+"END ATOM"+"\n";

	lines.push_back(begin_header);

	for(core::Size index = 1; index <= residue->natoms(); ++index)
	{
		core::Vector xyz_coords(residue->xyz(index));
		core::chemical::AtomType const atom_type = residue->atom_type(index);
		core::chemical::ResidueTypeCOP residue_type = & residue->type();
		std::string element = atom_type.element();
		core::Real charge = residue_type->atom(index).charge();

		std::string atom_string = line_header_ + " " + utility::to_string<core::Size>(index)+" "+
									element+" "+ utility::to_string<core::Real>(xyz_coords.x())+ " "+
									utility::to_string<core::Real>(xyz_coords.y())+ " "	+
									utility::to_string<core::Real>(xyz_coords.z())+ " " +
									"0"+" "+"CHG="+utility::to_string<core::Real>(charge)+"\n";
		lines.push_back(atom_string);
	}
	lines.push_back(end_header);
	return lines;
}

std::list<std::string> MolWriter::compose_bonds(core::conformation::ResidueCOP residue)
{
	std::list<std::string> lines;
	std::string begin_header = line_header_+"BEGIN BOND"+"\n";
	std::string end_header = line_header_+"END BOND"+"\n";

	lines.push_back(begin_header);

	//bonds are non-directional so we want to make sure they only get counted once.
	//(bond from atom 1 to atom 5 is the same as the one from 5 to 1)
	//stick all the bonds in a set and then iterate through the set to output
	//the smaller atom index is always designated lower in the struct


	std::set<BondData> bond_data_set;
	bond_data_set.clear();
	for(core::Size index = 1; index <= residue->natoms();++index)
	{
		core::chemical::AtomIndices const bonded_neighbors = residue->bonded_neighbor(index);
		utility::vector1<core::chemical::BondName> const bonded_neighbor_types = residue->type().bonded_neighbor_types(index);
		assert(bonded_neighbors.size()== bonded_neighbor_types.size());

		for(core::Size neighbor_index = 1; neighbor_index <= bonded_neighbors.size();++neighbor_index)
		{
			core::Size type = bonded_neighbor_types[neighbor_index];
			core::Size neighbor = bonded_neighbors[neighbor_index];
			BondData bond(index,neighbor,type);
			if(bond_data_set.find(bond)!= bond_data_set.end())
			{
				//This is a terrible hack. If it's not here, and there
				//is a bond from index 1 to index 4, it gets inserted twice
				//I don't know why, I'll fix it later, I promise.

				//Sorry :(
				continue;
			}
			bond_data_set.insert(bond);
		}
	}

	core::Size bond_index = 1;
	foreach(BondData current_bond, bond_data_set){
		std::string bond_line = line_header_ + utility::to_string<core::Size>(bond_index) + " " +
								utility::to_string<core::Size>(current_bond.bondType) + " " +
								utility::to_string<core::Size>(current_bond.lower) + " " +
								utility::to_string<core::Size>(current_bond.upper) + " " + "\n";
		lines.push_back(bond_line);
		++bond_index;
	}

	lines.push_back(end_header);

	return lines;
}

std::list<std::string> MolWriter::compose_typeinfo(core::conformation::ResidueCOP residue)
{
	std::list<std::string> lines;

	std::string header = "> <Rosetta AtomTypes>\n";
	std::string type_data = "";
	core::chemical::ResidueTypeCOP residue_type = & residue->type();
	for(core::Size index =1; index <= residue->natoms();++index)
	{
		std::string atom_type_name = residue_type->atom_type(index).name();
		std::string data_string = "("+utility::to_string<core::Size>(index)+","+atom_type_name+") ";
		type_data.append(data_string);
	}
	type_data.append("\n");

	lines.push_back(header);
	lines.push_back(type_data);
	lines.push_back("\n");

	return lines;
}

std::list<std::string> MolWriter::compose_nbr_atom(core::conformation::ResidueCOP residue)
{
	std::list<std::string>	lines;

	std::string header = "> <Rosetta nbr_atom>\n";
	std::string nbr_atom = utility::to_string<core::Size>(residue->nbr_atom()) + "\n";

	lines.push_back(header);
	lines.push_back(nbr_atom);
	lines.push_back("\n");

	return lines;
}

std::list<std::string> MolWriter::compose_job_info()
{
	std::list<std::string>	lines;
	for(std::map<std::string,std::string>::const_iterator data_it = job_data_.begin(); data_it != job_data_.end();++data_it)
	{
		std::string header_name(data_it->first);
		std::string data(data_it->second+"\n");
		std::string header("> <"+header_name+">\n");
		lines.push_back(header);
		lines.push_back(data);
		lines.push_back("\n");
	}
	return lines;
}


}
}
}
