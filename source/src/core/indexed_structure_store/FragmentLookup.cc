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

#include <numeric/xyzVector.hh>

#include <core/indexed_structure_store/FragmentLookup.hh>
#include <core/indexed_structure_store/FragmentStore.hh>

#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>

#include <vector>

namespace core {
namespace indexed_structure_store {

FragmentLookupResult::FragmentLookupResult() :
	found_match(false),
	match_score(0),
	match_index(0),
	match_rmsd(std::numeric_limits<numeric::Real>::infinity()),
	match_rmsd_threshold(std::numeric_limits<numeric::Real>::infinity())
{

}

FragmentLookup::FragmentLookup(FragmentStoreOP store) :
	store_(std::move(store)),
	lookup_(
	&(store_->fragment_coordinates[0].x()),
	&(store_->fragment_threshold_distances[0]),
	store_->fragment_threshold_distances.size(),
	store_->fragment_specification.coordinates_per_fragment())
{

}

std::vector< std::pair< Size, Size > > FragmentLookup::get_fragment_residue_spans(core::pose::Pose const & pose)
{
	// Search for all contiguous spans of protein residues in the pose structure

	// Fragment discontinuities fall into two categories:
	//   breaks - Chain endings where fragments spanning r and r+1 are invalid,
	//        but r and r+1 may be within a fragment.
	//   residues - Invalid residues where any fragment containing r is invalid.
	//
	//   For each break insert the break point, and the starting residue after
	//   the break point. For breaks, this will be the break residue, for invalid
	//   residues this will be the residue after the break.
	//
	//   std::map stores the results in increasing order.
	std::map<Size, Size> span_breaks_to_next_start;
	typedef std::map<Size, Size>::value_type SpanBreakValue;

	for ( Size i = 1; i <= pose.conformation().num_chains(); i++ ) {
		span_breaks_to_next_start[pose.conformation().chain_end(i) + 1] = pose.conformation().chain_end(i) + 1;
	}

	for ( Size i = 1; i <= pose.size(); i++ ) {
		if ( !pose.residue(i).is_protein() ) {
			span_breaks_to_next_start[i] = i + 1;
		}
	}

	if ( core::pose::symmetry::is_symmetric(pose) ) {
		conformation::symmetry::SymmetryInfoCOP symm_info = core::pose::symmetry::symmetry_info( pose );
		span_breaks_to_next_start[symm_info->last_independent_residue() + 1] = pose.size() + 1;
	}

	std::vector< ResidueSpan > valid_residue_spans;
	// Unpack endings into list of valid spans within the pose
	// where a span runs along [start, end).
	numeric::Size span_start = 1;
	for ( SpanBreakValue const & break_and_start : span_breaks_to_next_start ) {
		Size next_break = break_and_start.first;
		Size next_start = break_and_start.second;

		if ( next_break >= span_start ) {
			valid_residue_spans.emplace_back(span_start, next_break);
			span_start = next_start;
		}
	}

	return valid_residue_spans;
}

FragmentLookupResult FragmentLookup::lookup_fragment_v(std::vector< numeric::xyzVector<numeric::Real> > & query_coordinates)
{
	FragmentLookupResult result;

	numeric::Real* coordinate_start = &(query_coordinates[0].x());
	result.found_match = lookup_.first_match(coordinate_start, result.match_index, result.match_rmsd);
	if ( result.found_match ) {
		result.match_score = 1;
		result.match_rmsd_threshold = store_->fragment_threshold_distances[result.match_index];
	}
	return result;
}

FragmentLookupResult FragmentLookup::lookup_closest_fragment_v(std::vector< numeric::xyzVector<numeric::Real> > & query_coordinates)
{
	FragmentLookupResult result;
	numeric::Real* coordinate_start = &(query_coordinates[0].x());
	result.found_match = lookup_.closest_match(coordinate_start, result.match_index, result.match_rmsd);

	if ( result.found_match ) {
		result.match_score = 1;
		result.match_rmsd_threshold = store_->fragment_threshold_distances[result.match_index];
	}
	return result;
}

FragmentLookupResult FragmentLookup::lookup_closest_fragment_subset_v(std::vector< numeric::xyzVector<numeric::Real> > & query_coordinates,std::vector<bool> subset)
{
	FragmentLookupResult result;
	numeric::Real* coordinate_start = &(query_coordinates[0].x());
	result.found_match = lookup_.closest_match_subset(coordinate_start, result.match_index, result.match_rmsd,subset);

	if ( result.found_match ) {
		result.match_score = 1;
		result.match_rmsd_threshold = store_->fragment_threshold_distances[result.match_index];
	}
	return result;
}


std::vector<FragmentLookupResult> FragmentLookup::lookup_close_fragments_v(std::vector< numeric::xyzVector<numeric::Real> > & query_coordinates,  Real rms_threshold, Size max_matches)
{
	std::vector<FragmentLookupResult> results;
	std::vector<Size>match_indexes;
	std::vector<Real>distances;

	numeric::Real* coordinate_start = &(query_coordinates[0].x());

	lookup_.all_matches_below_threshold(coordinate_start, match_indexes, distances,rms_threshold);
	for ( Size ii=0; (ii<match_indexes.size() && results.size()<max_matches); ++ii ) {
		FragmentLookupResult result;
		result.match_index =match_indexes[ii];
		result.match_rmsd=distances[ii];
		result.found_match = true;
		results.push_back(result);
	}
	return(results);
}


}
}
