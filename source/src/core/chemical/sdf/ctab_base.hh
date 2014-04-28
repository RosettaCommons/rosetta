// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/chemical/sdf/ctab_base.hh
/// @author Sam DeLuca



#ifndef INCLUDED_core_chemical_sdf_ctab_base_hh
#define INCLUDED_core_chemical_sdf_ctab_base_hh

#include <core/chemical/sdf/ctab_base.fwd.hh>

#include <core/chemical/ResidueType.fwd.hh>
#include <utility/vector1.hh>
#include <core/chemical/sdf/MolData.hh>
// AUTO-REMOVED #include <core/chemical/sdf/mol_util.hh>

#include <core/chemical/sdf/mol_util.fwd.hh>
#include <set>


namespace core {
namespace chemical {
namespace sdf {

static std::string const DEFAULT_ATOM_TYPE_="CH3";
static std::string const DEFAULT_MM_ATOM_TYPE_="X";

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

elementToType & element_to_default_type();

class CtabBase {
public:
	CtabBase(utility::vector1<std::string> const & connection_table_lines, core::chemical::ResidueTypeOP molecule_container, MolData const & mol_data);

	virtual ~CtabBase();

	virtual void ParseTable() = 0;
	core::chemical::ResidueTypeOP GetResidueType();

	core::Size connection_table_length() const ;
	std::string connection_table_line(core::Size const line_number) const;

	void add_index_name_pair(core::Size const index, std::string const atomname);
	std::string atom_name_from_index(core::Size const index) const;

	bool check_for_aromatic(core::Size lower, core::Size upper);

	void set_atom_type(core::Size const atomno, std::string const atomname);
	void fix_atom_types();


private:
	virtual void ParseAtom(std::string const atom_line, core::Size const atom_number) = 0;
	virtual void ParseBond(std::string const bond_line) = 0;


private:

	utility::vector1<std::string> connection_table_lines_;
	core::chemical::ResidueTypeOP molecule_container_;
	MolData mol_data_;
	std::map<core::Size, std::string> index_to_names_map_;

	std::map<core::Size,std::string> atom_type_data_map_;
	std::set<BondData> bond_type_data_set_;

};

}
}
}


#endif /* CTAB_BASE_HH_ */
