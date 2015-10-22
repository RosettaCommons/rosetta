// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/chemical/sdf/mol_util.hh
/// @author Sam DeLuca

#ifndef INCLUDED_core_chemical_sdf_mol_util_hh
#define INCLUDED_core_chemical_sdf_mol_util_hh

#include <core/chemical/sdf/mol_util.fwd.hh>
#include <core/types.hh>
#include <map>
#include <set>
#include <iosfwd>

namespace core {
namespace chemical {
namespace sdf {

struct BondData
{
	BondData(core::Size index1, core::Size index2, core::Size type);
	bool operator<(const BondData& other) const;
	bool operator ==(const core::chemical::sdf::BondData & other) const;
	core::Size lower;
	core::Size upper;
	core::Size bondType;

};

std::set<BondData> parse_bond_type_data(std::string raw_data);
std::map<core::Size, std::string> parse_atom_type_data(std::string raw_data);
}
}
}


#endif /* MOL_UTIL_HH_ */
