// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   /grid/src/core/chemical/sdf/MolData.hh
/// @author Sam DeLuca

#ifndef INCLUDED_core_chemical_sdf_MolData_hh
#define INCLUDED_core_chemical_sdf_MolData_hh

#include <core/chemical/sdf/MolData.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <map>
#include <string>

namespace core {
namespace chemical {
namespace sdf {

class MolData
{
public:
	MolData();
	virtual ~MolData();
	void parse_mol_data(utility::vector1<std::string> const &  file_lines);
	core::Size size();
	void clear();
	std::string get_mol_data(std::string const & key) const;
	utility::vector1<std::string> get_mol_data_string_vector(std::string const & key,char const & splitter) const;
	void print() const;
private:
	std::map<std::string,std::string> mol_data_map_;

};

}
}
}

#endif /* MOLDATA_HH_ */
