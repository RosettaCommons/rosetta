// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file core/chemical/sdf/mol_parser.hh
///
/// @brief header file for molfile parser
/// @author Sam DeLuca



#ifndef INCLUDED_core_chemical_sdf_mol_parser_HH_
#define INCLUDED_core_chemical_sdf_mol_parser_HH_

#include <core/chemical/ResidueType.hh>
#include <utility/vector1.hh>
#include <core/chemical/sdf/mol_parser.fwd.hh>
#include <core/chemical/sdf/MolData.hh>
//#include <protocols/ligand_docking/ColoredGraph.fwd.hh>


#include <string>


namespace core {
namespace chemical {
namespace sdf {

class MolFileParser
{
public:

	MolFileParser(const utility::vector1<std::string> mol_file_lines);
	MolFileParser(const std::string file_name);
	void parse_mol_file(core::chemical::AtomTypeSetCAP atom_types, ElementSetCAP elements, core::chemical::MMAtomTypeSetCAP mm_atom_types, core::chemical::orbitals::OrbitalTypeSetCAP orbital_type_set);
	core::chemical::ResidueTypeOP GetResidueTypeOP();
	core::chemical::ResidueType GetResidueType();

	std::string GetMoleculeName();
	std::string GetMoleculeInfo();
	std::string GetMoleculeComments();

private:

	std::string FindNbrAtom();

	utility::vector1<std::string> mol_file_lines_;
	core::chemical::ResidueTypeOP molecule_container_;
	MolData mol_data_;
	//protocols::ligand_docking::ColoredGraphOP molecule_container_;
	std::string molecule_name_;
	std::string molecule_info_;
	std::string molecule_comments_;
};

} // sdf
} // io
} // core

#endif /* MOL_PARSER_HH_ */
