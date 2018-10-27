// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @file
/// @brief
/// @author Alex Ford <fordas@uw.edu>

#pragma once

#include <json.hpp>
#include <protocols/indexed_structure_store/search/QueryDatabase.hh>


namespace protocols { namespace indexed_structure_store { namespace search {
inline void to_json(nlohmann::json& j, const SingleQuerySummaryStatistics& s) {
	j = nlohmann::json{
		{"structures_considered", s.structures_considered},
		{"fragments_considered", s.fragments_considered},
		{"result_count", s.result_count}
		};
}

inline void from_json(const nlohmann::json& j, SingleQuerySummaryStatistics& s) {
	s.structures_considered = j["structures_considered"].get<Index>();
	s.fragments_considered = j["fragments_considered"].get<Index>();
	s.result_count = j["result_count"].get<Index>();
}

inline void to_json(nlohmann::json& j, const StructureSingleQueryResult& s) {
	j = nlohmann::json{
		{"fragment_start", s.fragment_start},
		{"structure_index", s.structure_index},
		{"fragment_structure_start", s.fragment_structure_start},
		{"result_rmsd", s.result_rmsd}
		};
}

inline void from_json(const nlohmann::json& j, StructureSingleQueryResult& s) {
	s.fragment_start = j["fragment_start"].get<Index>();
	s.structure_index = j["structure_index"].get<Index>();
	s.fragment_structure_start = j["fragment_structure_start"].get<Index>();
	s.result_rmsd = j["result_rmsd"].get<SearchReal>();
}

inline void to_json(nlohmann::json& j, const PairQuerySummaryStatistics& s) {
	j = nlohmann::json{
		{"structures_considered", s.structures_considered},
		{"fragments_considered", s.fragments_considered},
		{"fragments_expanded", s.fragments_expanded},
		{"pairs_considered", s.pairs_considered},
		{"pairs_aligned", s.pairs_aligned},
		{"result_count", s.result_count}
		};
}

inline void from_json(const nlohmann::json& j, PairQuerySummaryStatistics& s) {
	s.structures_considered = j["structures_considered"].get<Index>();
	s.fragments_considered = j["fragments_considered"].get<Index>();
	s.fragments_expanded = j["fragments_expanded"].get<Index>();
	s.pairs_considered = j["pairs_considered"].get<Index>();
	s.pairs_aligned = j["pairs_aligned"].get<Index>();
	s.result_count = j["result_count"].get<Index>();
}

inline void to_json(nlohmann::json& j, const StructurePairQueryResult& s) {
	j = nlohmann::json{
		{"fragment_a_start", s.fragment_a_start},
		{"fragment_b_start", s.fragment_b_start},
		{"structure_index", s.structure_index},
		{"fragment_a_structure_start", s.fragment_a_structure_start},
		{"fragment_b_structure_start", s.fragment_b_structure_start},
		{"result_rmsd", s.result_rmsd}
		};
}

inline void from_json(const nlohmann::json& j, StructurePairQueryResult& s) {
	s.fragment_a_start = j["fragment_a_start"].get<Index>();
	s.fragment_b_start = j["fragment_b_start"].get<Index>();
	s.structure_index = j["structure_index"].get<Index>();
	s.fragment_a_structure_start = j["fragment_a_structure_start"].get<Index>();
	s.fragment_b_structure_start = j["fragment_b_structure_start"].get<Index>();
	s.result_rmsd = j["result_rmsd"].get<SearchReal>();
}

} } }
