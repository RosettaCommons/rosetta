// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   core/chemical/sdf/ctab_parsher.hh
///
/// @brief parse the CTAB table into a residue object
/// @author Sam DeLuca


#ifndef INCLUDED_CORE_chemical_sdf_ctab_parser_HH
#define INCLUDED_CORE_chemical_sdf_ctab_parser_HH

#include <core/chemical/sdf/ctab_parser.fwd.hh>
#include <core/chemical/sdf/ctab_typer.hh>
#include <core/chemical/sdf/ctab_base.hh>
#include <core/chemical/sdf/MolData.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <utility/vector1.hh>

#include <string>
#include <map>
#include <vector>


namespace core {
namespace chemical {
namespace sdf {

static std::string const DEFAULT_ATOM_TYPE_="CH3";
static std::string const DEFAULT_MM_ATOM_TYPE_="X";

/*
struct addedH {
public:
	core::Size atom_number;
	std::string bonded_atom_name;
	std::string atom_type;
};


class elementToType {
public:
	elementToType();
	std::string get(std::string key);
private:
	std::map<std::string, std::string> e_to_t;
};
*/
//static elementToType element_to_default_type;

class ctabV2000Parser {

public:
	ctabV2000Parser(const utility::vector1<std::string> connection_table_lines, core::chemical::ResidueTypeOP molecule_container,MolData mol_data);
	void ParseTable();

	core::chemical::ResidueTypeOP GetResidueType();
	//protocols::ligand_docking::ColoredGraphOP GetResidueGraph();


private:
	void set_atom_type(core::Size atomno, std::string atomname);
	void ParseAtom(const std::string atom_line, core::Size atom_number);
	void ParseBond(const std::string bond_line);
	std::map<core::Size,std::string> ParseAtomTypeData();

private:
	utility::vector1<std::string> connection_table_lines_;
	MolData mol_data_;
	std::map<core::Size, std::string> index_to_names_map_;
	core::chemical::ResidueTypeOP molecule_container_;
	utility::vector1<addedH> added_H_;
	core::Size current_atom_;

};

class ctabV3000Parser {
public:
	ctabV3000Parser(const utility::vector1<std::string> connection_table_lines, core::chemical::ResidueTypeOP molecule_container, MolData mol_data );

	void ParseTable();

	core::chemical::ResidueTypeOP GetResidueType();
	//protocols::ligand_docking::ColoredGraphOP GetResidueGraph();


private:
	void set_atom_type(core::Size atomno, std::string atomname);
	void ParseAtom(const std::string atom_line);
	void ParseBond(const std::string bond_line );
	core::Real FindExtraParameter( std::vector<std::string> extra_parameters, const std::string query);
	std::map<core::Size,std::string> ParseAtomTypeData();

private:
	utility::vector1<std::string> connection_table_lines_;
	MolData mol_data_;
	std::map<core::Size, std::string> index_to_names_map_;
	core::chemical::ResidueTypeOP molecule_container_;
	utility::vector1<addedH> added_H_;
	core::Size current_atom_;

};



}
}
}


#endif /* CTAB_PARSER_HH_ */
