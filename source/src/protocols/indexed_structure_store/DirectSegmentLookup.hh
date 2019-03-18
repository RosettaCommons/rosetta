// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ./src/protocols/indexed_structure_store/DirectSegmentLookup.hh
/// @brief Support class for direct rmsd-based segment lookup.
/// @details
/// @author Alex Ford (fordas@uw.edu)

#pragma once

#include <ndarray.h>
#include <numeric/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/indexed_structure_store/StructureStore.fwd.hh>
#include <protocols/indexed_structure_store/search/QueryDatabase.hh>

namespace protocols { namespace indexed_structure_store {

struct DirectSegmentLookupConfig {
	search::SearchReal rmsd_tolerance;
	search::SearchReal segment_cluster_tolerance;
	search::Index max_insertion_length;
};

struct DirectSegmentLookupResult {
	std::vector<ResidueEntry> result_residues;
	std::vector<search::StructurePairQueryResult> query_results;
};

struct DirectSegmentLookup {
	typedef DirectSegmentLookupConfig Config;

	DirectSegmentLookup(Config conf) : config(conf) { }

	std::vector<DirectSegmentLookupResult> segment_lookup(
		ndarray::Array<ResidueEntry, 1> source_residues,
		search::StructureDatabase & structure_db,
		core::pose::Pose & context,
		platform::Size n_start_res, platform::Size n_end_res,
		platform::Size c_start_res, platform::Size c_end_res
	);

	search::PairQuerySummaryStatistics query_stats;
	DirectSegmentLookupConfig config;
};

} }
