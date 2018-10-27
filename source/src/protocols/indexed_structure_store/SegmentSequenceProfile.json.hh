// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ./src/protocols/indexed_structure_store/SegmentSequenceProfile.json.hh
/// @brief Json-support for segment sequence profile configuration.
/// @author Alex Ford (fordas@uw.edu)
//

#pragma once

#include <json.hpp>
#include <protocols/indexed_structure_store/SegmentSequenceProfile.hh>


namespace protocols
{
namespace indexed_structure_store
{

inline void to_json(nlohmann::json& j, const SegmentSequenceProfileConfig& c) {
	j = nlohmann::json{
		{"rmsd_tolerance" , c.rmsd_tolerance},
		{"pseudocount" , c.pseudocount},
		};
}

inline void from_json(const nlohmann::json& j, SegmentSequenceProfileConfig& c) {
	c.rmsd_tolerance = j["rmsd_tolerance"].get<search::SearchReal>();
	c.pseudocount = j["pseudocount"].get<search::SearchReal>();
}

}
}
