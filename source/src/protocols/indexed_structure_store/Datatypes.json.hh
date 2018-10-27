// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @file protocols/indexed_structure_store/Datatypes.json.hh
/// @brief json-support functions for indexed_structure_store primitive datatypes.
/// @author Alex Ford <fordas@uw.edu>

#pragma once

#include <algorithm>
#include <string>
#include <json.hpp>
#include <protocols/indexed_structure_store/Datatypes.hh>


namespace protocols
{
namespace indexed_structure_store
{

inline void to_json(nlohmann::json& j, const ResidueBackboneEntry& i) {
	j = nlohmann::json{
		{"phi" , i.phi},
		{"psi" , i.psi},
		{"omega" , i.omega},
		};
}

inline void from_json(const nlohmann::json& j, ResidueBackboneEntry& i) {
	i.phi = j["phi"].get<float>();
	i.psi = j["psi"].get<float>();
	i.omega = j["omega"].get<float>();
}

inline void to_json(nlohmann::json& j, const ResidueSidechainEntry& i) {
	j = nlohmann::json{
		{"chi1" , i.chi1},
		{"chi2" , i.chi2},
		{"chi3" , i.chi3},
		{"chi4" , i.chi4},
		{"aa" , std::string(1, i.aa)},
		};
}

inline void from_json(const nlohmann::json& j, ResidueSidechainEntry& i) {
	i.chi1 = j["chi1"].get<float>();
	i.chi2 = j["chi2"].get<float>();
	i.chi3 = j["chi3"].get<float>();
	i.chi4 = j["chi4"].get<float>();
	i.aa = j["aa"].get<std::string>()[0];
}

inline void to_json(nlohmann::json& j, const ResidueOrientEntry& i) {
	j = nlohmann::json{
		{"N" , i.N},
		{"C" , i.C},
		{"CA" , i.CA},
		{"O" , i.O},
		};
}

inline void from_json(const nlohmann::json& j, ResidueOrientEntry& i) {
	//TODO alexford add array accessors
	i.N[0] = j["N"][0].get<float>();
	i.N[1] = j["N"][1].get<float>();
	i.N[2] = j["N"][2].get<float>();
	i.C[0] = j["C"][0].get<float>();
	i.C[1] = j["C"][1].get<float>();
	i.C[2] = j["C"][2].get<float>();
	i.CA[0] = j["CA"][0].get<float>();
	i.CA[1] = j["CA"][1].get<float>();
	i.CA[2] = j["CA"][2].get<float>();
	i.O[0] = j["O"][0].get<float>();
	i.O[1] = j["O"][1].get<float>();
	i.O[2] = j["O"][2].get<float>();
}

inline void to_json(nlohmann::json& j, const ResidueEntry& i) {
	j = nlohmann::json{
		{"structure_id" , i.structure_id},
		{"residue_id" , i.residue_id},
		{"bb" , i.bb},
		{"sc" , i.sc},
		{"orient" , i.orient},
		{"chain_ending" , i.chain_ending},
		};
}

inline void from_json(const nlohmann::json& j, ResidueEntry& i) {
	i.structure_id = j["structure_id"].get<uint32_t>();
	i.residue_id = j["residue_id"].get<uint32_t>();
	i.bb = j["bb"].get<ResidueBackboneEntry>();
	i.sc = j["sc"].get<ResidueSidechainEntry>();
	i.orient = j["orient"].get<ResidueOrientEntry>();
	i.chain_ending = j["chain_ending"].get<bool>();
}

inline void to_json(nlohmann::json& j, const StructureEntry& i) {
	std::string name(i.name, 32);
	name.erase(name.find_first_of('\0'));

	j = nlohmann::json{
		{"name" , name},
		{"id" , i.id}
		};
}

inline void from_json(const nlohmann::json& j, StructureEntry& i) {
	i.id = j["id"].get<uint32_t>();
	std::string name = j["name"].get<std::string>();
	name.resize(32);
	name.copy(i.name, 32, 0);
}
}
}
