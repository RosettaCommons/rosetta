// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file core/chemical/sdf/sdf_parser.hh
///
/// @brief sdf file parser header
/// @author Sam DeLuca


#ifndef INCLUDED_core_chemical_sdf_sdf_parser_HH_
#define INCLUDED_core_chemical_sdf_sdf_parser_HH_

#include <core/chemical/sdf/sdf_parser.fwd.hh>
#include <core/chemical/ResidueType.hh>


#include <map>
#include <string>
#include <utility/vector1.hh>

namespace core {
namespace chemical {
namespace sdf {

class SDFParser {
public:
	SDFParser(utility::vector1<std::string> file_vector);
	void SplitSDF();
	core::Size GetNumberOfResidues();
	void GenerateResidues();
	core::chemical::ResidueTypeOP GetResidueByName(std::string name);
	utility::vector1<std::string> GetResidueNameVector();

private:
	std::map<std::string,utility::vector1<std::string> > extended_data_map_;
	std::map<std::string,utility::vector1<std::string> > mol_file_map_;
	utility::vector1<std::string> file_vector_;

	std::map<std::string,core::chemical::ResidueTypeOP> molecule_map_;
	utility::vector1<std::string> molecule_name_vector_;
};

}
}
}

#endif /* SDF_PARSER_HH_ */
