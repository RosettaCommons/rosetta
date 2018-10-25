// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ./src/protocols/indexed_structure_store/SegmentSequenceProfile.cc
/// @brief Generate lookup-based sequence profiles for contiguous structureal elements.
/// @author Alex Ford (fordas@uw.edu)
//
#include <vector>
#include <numeric>

#include <boost/range/combine.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range.hpp>

#include <utility/backtrace.hh>
#include <utility/exit.hh>
#undef eigen_assert
#define eigen_assert(x) \
	runtime_assert(x)

#include <json.hpp>

// Project headers
#include <basic/Tracer.hh>

#include <numeric/alignment/rmsd_calc.hh>

#include <core/pose/xyzStripeHashPose.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/id/types.hh>
#include <core/conformation/Residue.hh>

#include <protocols/indexed_structure_store/Datatypes.hh>
#include <protocols/indexed_structure_store/Datatypes.json.hh>
#include <protocols/indexed_structure_store/StructureStore.hh>
#include <protocols/indexed_structure_store/search/QueryDatabase.hh>
#include <protocols/indexed_structure_store/search/QueryDatabase.json.hh>
#include <protocols/indexed_structure_store/utility.hh>
#include <protocols/indexed_structure_store/orient_array.hh>
#include <protocols/indexed_structure_store/pose_utility.hh>
#include <protocols/indexed_structure_store/vector_tools.hh>

#include <protocols/indexed_structure_store/SegmentSequenceProfile.hh>
#include <protocols/indexed_structure_store/SegmentSequenceProfile.json.hh>

#include "ndarray.h"
#include "ndarray/eigen.h"

static basic::Tracer TR("protocols.indexed_structure_store.SegmentSequenceProfile");

using nlohmann::json;

namespace protocols { namespace indexed_structure_store {

using namespace boost::adaptors;
using namespace boost::range;

using namespace protocols::indexed_structure_store;
using namespace protocols::indexed_structure_store::search;

std::map<core::chemical::AA, core::Real> SegmentSequenceProfile::aa_background_distribution = {
{ core::chemical::aa_ala, 0.082949},
{ core::chemical::aa_arg, 0.046585},
{ core::chemical::aa_asn, 0.046900},
{ core::chemical::aa_asp, 0.059111},
{ core::chemical::aa_cys, 0.017223},
{ core::chemical::aa_gln, 0.037391},
{ core::chemical::aa_glu, 0.061225},
{ core::chemical::aa_gly, 0.079431},
{ core::chemical::aa_his, 0.022011},
{ core::chemical::aa_ile, 0.055373},
{ core::chemical::aa_leu, 0.082300},
{ core::chemical::aa_lys, 0.059822},
{ core::chemical::aa_met, 0.021037},
{ core::chemical::aa_phe, 0.040169},
{ core::chemical::aa_pro, 0.046761},
{ core::chemical::aa_ser, 0.061387},
{ core::chemical::aa_thr, 0.058805},
{ core::chemical::aa_trp, 0.014722},
{ core::chemical::aa_tyr, 0.036899},
{ core::chemical::aa_val, 0.069899}
};

SegmentSequenceProfileResult SegmentSequenceProfile::segment_profile(
	StructureStore & structure_store,
	StructureDatabase & structure_db,
	core::pose::Pose & context,
	core::Size segment_start_resi, core::Size segment_end_resi
) {
	using ArrayXaa = SegmentSequenceProfileResult::ArrayXaa;
	using Array1aa = SegmentSequenceProfileResult::Array1aa;

	TR << "config: " << json(config) << " "
		<< "segment_start_resi: " << segment_start_resi << " "
		<< "segment_end_resi: " << segment_end_resi << " " << std::endl;

	core::Size segment_len = segment_end_resi - segment_start_resi;

	ndarray::Array<ResidueEntry, 1, 1> segment_res(segment_len);
	for ( core::Size i = 0; i < segment_len; ++i ) {
		segment_res[i] = extract_residue_entry(context.residue(segment_start_resi + i));
	}

	ndarray::Array<float, 3, 3> segment_orient = ndarray::copy(orient_array(segment_res));

	StructureSingleQuery query(segment_orient, config.rmsd_tolerance);

	SingleQueryExecutor executor(query);
	executor.execute(structure_db);

	TR << "query_stats:" << json(executor.query_stats) << std::endl;

	SegmentSequenceProfileResult result;
	result.query_results = executor.query_results;
	result.counts = ArrayXaa::Zero(segment_len, core::chemical::num_canonical_aas);

	for ( StructureSingleQueryResult qresult : result.query_results ) {
		for ( core::Size i = 0; i < segment_len; ++i ) {
			char aa = structure_store.residue_entries[qresult.fragment_start + i].sc.aa;
			result.counts(i, core::chemical::aa_from_oneletter_code(aa) - core::chemical::first_l_aa) += 1;
		}
	}

	Array1aa aa_bg_dist = Array1aa::Zero();
	for ( auto aa_freq : aa_background_distribution ) {
		aa_bg_dist[aa_freq.first - core::chemical::first_l_aa] = aa_freq.second;
	}
	Array1aa pseudocounts = (aa_bg_dist * config.pseudocount * core::chemical::num_canonical_aas);
	ArrayXaa mod_counts = result.counts;
	mod_counts.rowwise() += pseudocounts;

	result.frequencies = mod_counts.colwise() / mod_counts.rowwise().sum();
	result.log_odds = result.frequencies;
	result.log_odds.rowwise() /= aa_bg_dist;
	result.log_odds = result.log_odds.log() / std::log(2);

	// transformed range adaptor is broken on icc 17
#ifndef __INTEL_COMPILER
	TR.Debug << "counts:\n";
	boost::copy(
		boost::irange((int)core::chemical::first_l_aa, (int)core::chemical::num_canonical_aas + 1) |
		transformed([](int i){ return core::chemical::oneletter_code_from_aa((core::chemical::AA)i);}),
		std::ostream_iterator<char>(TR.Debug, " ")
	);
	TR.Debug << "\n" << result.counts << std::endl;

	TR.Debug << "frequencies:\n";
	boost::copy(
		boost::irange((int)core::chemical::first_l_aa, (int)core::chemical::num_canonical_aas + 1) |
		transformed([](int i){ return core::chemical::oneletter_code_from_aa((core::chemical::AA)i);}),
		std::ostream_iterator<char>(TR.Debug, " ")
	);
	TR.Debug << "\n" << result.frequencies << std::endl;

	TR.Debug << "log_odds:\n";
	boost::copy(
		boost::irange((int)core::chemical::first_l_aa, (int)core::chemical::num_canonical_aas + 1) |
		transformed([](int i){ return core::chemical::oneletter_code_from_aa((core::chemical::AA)i);}),
		std::ostream_iterator<char>(TR.Debug, " ")
	);
	TR.Debug << "\n" << result.log_odds << std::endl;
#endif

	return result;
}

} }
