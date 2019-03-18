// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ./src/protocols/indexed_structure_store/SegmentSequenceProfile.hh
/// @brief Generate lookup-based sequence profiles for contiguous structureal elements.
/// @author Alex Ford (fordas@uw.edu)
//

#pragma once

#include <ndarray.h>

#include <core/chemical/AA.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

#include <protocols/indexed_structure_store/StructureStore.fwd.hh>
#include <protocols/indexed_structure_store/search/QueryDatabase.hh>

namespace protocols { namespace indexed_structure_store {


struct SegmentSequenceProfileConfig {
	search::SearchReal rmsd_tolerance = .5;
	search::SearchReal pseudocount = 1;
};

struct SegmentSequenceProfileResult {
	typedef Eigen::Array<core::Real, Eigen::Dynamic, core::chemical::num_canonical_aas> ArrayXaa;
	typedef Eigen::Array<core::Real, 1, core::chemical::num_canonical_aas> Array1aa;

	std::vector<search::StructureSingleQueryResult> query_results;
	ArrayXaa counts;
	ArrayXaa frequencies;
	ArrayXaa log_odds;
};

struct SegmentSequenceProfile {

	typedef SegmentSequenceProfileConfig Config;

	static std::map<core::chemical::AA, core::Real> aa_background_distribution;

	SegmentSequenceProfile() = default;
	SegmentSequenceProfile(Config conf) : config(conf) { }

	SegmentSequenceProfileResult segment_profile(
		StructureStore & structure_store,
		search::StructureDatabase & structure_db,
		core::pose::Pose & context,
		core::Size segment_start, core::Size segment_end
	);

	search::PairQuerySummaryStatistics query_stats;
	SegmentSequenceProfileConfig config;
};
} }
