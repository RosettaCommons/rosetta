// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/chemical/sdf/mol_util.cc
/// @author Sam DeLuca

#include <core/chemical/sdf/mol_util.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <map>
#include <set>
#include <core/types.hh>

namespace core {
namespace chemical {
namespace sdf {

BondData::BondData(core::Size index1, core::Size index2, core::Size type)
{
	bondType = type;
	if ( index1 <= index2 ) {
		lower=index1;
		upper=index2;
	} else {
		lower=index2;
		upper=index1;
	}
}

bool BondData::operator <(const core::chemical::sdf::BondData & other) const
{
	return (this->lower < other.lower) || (this->upper < other.upper);
}

bool BondData::operator ==(const core::chemical::sdf::BondData& other) const
{
	return (this->lower == other.lower) && (this->upper == other.upper) && (this->bondType == other.bondType);
}

std::set<BondData> parse_bond_type_data(std::string raw_data)
{
	std::set<BondData> bond_set;
	if ( raw_data == "" ) {
		return bond_set;
	}
	utility::vector1<std::string> tokens(utility::string_split(raw_data,'\n'));
	if ( tokens.size() == 0 ) {
		return bond_set;
	}
	for ( core::Size index = 1; index <= tokens.size(); ++index ) {
		utility::vector1<std::string> current_token(utility::split(tokens[index]));
		if ( current_token.size() < 3 ) {
			continue;
		}

		core::Size lower_id = utility::from_string(current_token[1],core::Size(0));
		core::Size upper_id = utility::from_string(current_token[2],core::Size(0));
		core::Size type = utility::from_string(current_token[3],core::Size(0));

		bond_set.insert(BondData(lower_id,upper_id,type));
	}
	return bond_set;
}

std::map<core::Size,std::string> parse_atom_type_data(std::string raw_data)
{
	std::map<core::Size, std::string> data_map;
	if ( raw_data == "" ) {
		return data_map;
	}
	utility::vector1<std::string> tokens(utility::string_split(raw_data,' '));
	if ( tokens.size() == 0 ) {
		return data_map;
	} else {
		//utility::vector1<std::string> tokens=utility::split(atom_type_data);
		for ( core::Size index = 1; index <= tokens.size(); ++index ) {

			std::string current_token = tokens[index];
			if ( current_token.size() <=1 ) {
				continue;
			}
			utility::vector1<std::string> token_split=utility::string_split(current_token,',');
			utility::trim(token_split[1],"(");
			utility::trim(token_split[2],")");
			//std::cout << current_token<<std::endl;
			core::Size atomno = atoi(token_split[1].c_str());
			std::pair<core::Size, std::string> atom_type_point(atomno,token_split[2]);
			data_map.insert(atom_type_point);
		}

	}
	return data_map;
}


}
}
}
