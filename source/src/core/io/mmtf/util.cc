// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/io/mmtf/util.cc
/// @brief Functions for MMTF writing/reading
/// @author Daniel Farrell (danpf@uw.edu)

#include <core/io/mmtf/util.hh>

#include <utility/string_util.hh>

#include <mmtf.hpp>

namespace core {
namespace io {
namespace mmtf {


bool
sd_index::operator==(sd_index const & rhs) const {
	return std::tie(model_index, chain_index, group_index, group_atom_index) ==
		std::tie(rhs.model_index, rhs.chain_index, rhs.group_index, rhs.group_atom_index);
}


std::string
sd_index::to_string() const {
	std::stringstream ss;
	ss << "mmtf::sd_index mi: " << model_index << " ci: " << chain_index << " gi: " << group_index << " gai: " << group_atom_index;
	return ss.str();
}

std::ostream & operator<<( std::ostream & os, sd_index const & sd_i ) {
	os << sd_i.to_string();
	return os;
}


std::map<core::Size, sd_index>
make_atom_num_to_sd_map(::mmtf::StructureData const & sd) {
	int32_t modelIndex = 0;
	int32_t chainIndex = 0;
	int32_t groupIndex = 0;
	core::Size atomSerial = 0;

	std::map<core::Size, sd_index> ret_map;

	for ( int32_t i = 0; i < sd.numModels; i++, modelIndex++ ) {
		// traverse chains
		for ( int32_t j = 0; j < sd.chainsPerModel[modelIndex]; j++, chainIndex++ ) {
			// traverse groups
			for ( int32_t k = 0; k < sd.groupsPerChain[chainIndex]; k++, groupIndex++ ) {
				const ::mmtf::GroupType& group =
					sd.groupList[sd.groupTypeList[groupIndex]];
				int const groupAtomCount = group.atomNameList.size();
				for ( int32_t l = 0; l < groupAtomCount; l++, atomSerial++ ) {
					ret_map[atomSerial] = sd_index(modelIndex, chainIndex, groupIndex, l);
				}
			}
		}
	}
	return ret_map;
}

std::string
resid_to_tag(ResID const & resid) {
	// Not the historical one, but something which is more compatible with multi-chain designations
	return std::to_string(resid.seq()) + "|" + std::string{resid.icode()} + "|" + resid.chain();
}

/// @brief Convert an mmtf-read tag into the internal ResID structure
ResID
resid_from_tag(std::string const & input) {
	if ( input.find('|') != std::string::npos ) {
		utility::vector1< std::string > split = utility::string_split(input, '|');
		if ( split.size() != 3 ) {
			utility_exit_with_message("Issue reading MMTF file - can't understand `"+input+"` as resid");
		}
		return ResID( atof(split[1].c_str()), split[2][0], split[3] );
	} else {
		// Older, smashed together version
		if ( input.size() != 6 ) {
			utility_exit_with_message("Issue reading MMTF file - can't understand `"+input+"` as resid");
		}
		return ResID( atof( input.substr(1,4).c_str() ), input[4], input.substr(5,1) );
	}
}

} // mmtf
} // io
} // core

