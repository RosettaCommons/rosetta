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
#include <protocols/indexed_structure_store/DirectSegmentLookup.hh>


namespace protocols
{
namespace indexed_structure_store
{

inline void to_json(nlohmann::json& j, const DirectSegmentLookupConfig& c) {
	j = nlohmann::json{
		{"rmsd_tolerance" , c.rmsd_tolerance},
		{"segment_cluster_tolerance" , c.segment_cluster_tolerance},
		{"max_insertion_length" , c.max_insertion_length}
		};
}

inline void from_json(const nlohmann::json& j, DirectSegmentLookupConfig& c) {
	c.rmsd_tolerance = j["rmsd_tolerance"].get<search::SearchReal>();
	c.segment_cluster_tolerance = j["segment_cluster_tolerance"].get<search::SearchReal>();
	c.max_insertion_length = j["max_insertion_length"].get<search::Index>();
}

}
}
