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

#ifndef INCLUDED_core_indexed_structure_store_FragmentLookup_hh
#define INCLUDED_core_indexed_structure_store_FragmentLookup_hh

#include <iterator>
#include <vector>
#include <map>
#include <limits>

#include <utility/pointer/ReferenceCount.hh>

#include <numeric/types.hh>
#include <numeric/xyzVector.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/indexed_structure_store/FragmentStore.hh>
#include <core/indexed_structure_store/FragmentLookup.fwd.hh>

#include <numeric/coordinate_fitting/CoordinateArray_RMSD_FlatLookup.hh>
#include <vector>

namespace core
{
namespace indexed_structure_store
{

class FragmentLookupResult
{
public:
	FragmentLookupResult();

	// @brief True if lookup resulted in match.
	//
	// All other values are undefined if false.
	bool found_match;

	// @brief Match score on range [0,1], higher is more favorable.
	numeric::Real match_score;

	// @brief Matching fragment index in lookup.
	numeric::Size match_index;

	// @brief Match fragment RMSD.
	numeric::Real match_rmsd;
	// @brief Matching fragment threshold.
	numeric::Real match_rmsd_threshold;
};

class FragmentLookup : public utility::pointer::ReferenceCount
{
public:
	// @brief Basic structure store, holds a collection of structure and associated residue entries.
	FragmentLookup(FragmentStoreOP store);

	// @brief Copy fragment_specification().fragment_length coordinates from input vector and lookup.
	template<class xyzVectorIterator>
	FragmentLookupResult lookup_fragment(xyzVectorIterator input)
	{
		// Lookup modifies coordinates in-place during lookup, copy input vector for query.
		std::vector< numeric::xyzVector<numeric::Real> > query_coordinates(store_->fragment_specification.coordinates_per_fragment());

		for ( Size i = 0; i < store_->fragment_specification.coordinates_per_fragment(); i++, input++ ) {
			query_coordinates[i] = *input;
		}

		return lookup_fragment(query_coordinates);
	}

	// @brief Copy fragment_specification().fragment_length coordinates from input vector and lookup.
	template<class xyzVectorIterator>
	FragmentLookupResult lookup_closest_fragment(xyzVectorIterator input)
	{
		// Lookup modifies coordinates in-place during lookup, copy input vector for query.
		std::vector< numeric::xyzVector<numeric::Real> > query_coordinates(store_->fragment_specification.coordinates_per_fragment());

		for ( Size i = 0; i < store_->fragment_specification.coordinates_per_fragment(); i++, input++ ) {
			query_coordinates[i] = *input;
		}

		return lookup_closest_fragment(query_coordinates);
	}

	// @brief Copy fragment_specification().fragment_length coordinates from input vector and lookup.
	template<class xyzVectorIterator>
	FragmentLookupResult lookup_closest_fragment_subset(xyzVectorIterator input, std::vector<bool> subset)
	{
		// Lookup modifies coordinates in-place during lookup, copy input vector for query.
		std::vector< numeric::xyzVector<numeric::Real> > query_coordinates(store_->fragment_specification.coordinates_per_fragment());

		for ( Size i = 0; i < store_->fragment_specification.coordinates_per_fragment(); i++, input++ ) {
			query_coordinates[i] = *input;
		}

		return lookup_closest_fragment_subset(query_coordinates,subset);
	}

	// @brief Copy fragment_specification().fragment_length coordinates from input vector and lookup.
	template<class xyzVectorIterator> std::vector<FragmentLookupResult> lookup_close_fragments(xyzVectorIterator input,Real rms_threshold)
	{
		// Lookup modifies coordinates in-place during lookup, copy input vector for query.
		std::vector< numeric::xyzVector<numeric::Real> > query_coordinates(store_->fragment_specification.coordinates_per_fragment());

		for ( Size i = 0; i < store_->fragment_specification.coordinates_per_fragment(); i++, input++ ) {
			query_coordinates[i] = *input;
		}

		return lookup_close_fragments(query_coordinates,rms_threshold);
	}

	template<class FragmentLookupOutputIterator, class ResidueNumberOutputIterator>
	void lookup_pose_fragments(core::pose::Pose const & pose, FragmentLookupOutputIterator result_out, ResidueNumberOutputIterator fragment_start_out)
	{
		std::vector< ResidueSpan > valid_residue_spans = get_fragment_residue_spans(pose);

		// Traverse each residue span, extract fragment atom coordinates and score.
		for ( ResidueSpan const & residue_span : valid_residue_spans ) {
			// Short-circuit if the span is shorter than the minimum fragment length.
			if ( residue_span.second - residue_span.first < store_->fragment_specification.fragment_length ) {
				continue;
			}

			std::vector< numeric::xyzVector<numeric::Real> > query_coordinates;
			for ( Size i = residue_span.first; i < residue_span.second; i++ ) {
				for ( std::string const & atom_name : store_->fragment_specification.fragment_atoms ) {
					query_coordinates.push_back(pose.residue(i).xyz(atom_name));
				}
			}

			for ( Size i = 0; residue_span.first + i + (store_->fragment_specification.fragment_length - 1) < residue_span.second; i++ ) {
				*result_out = lookup_fragment(
					&query_coordinates[i * store_->fragment_specification.fragment_atoms.size()]);
				result_out++;

				*fragment_start_out = residue_span.first + i;
				fragment_start_out++;
			}
		}
	}

	template<class FragmentLookupOutputIterator, class ResidueNumberOutputIterator>
	void lookup_closest_pose_fragments(core::pose::Pose const & pose, FragmentLookupOutputIterator result_out, ResidueNumberOutputIterator fragment_start_out)
	{
		std::vector< ResidueSpan > valid_residue_spans = get_fragment_residue_spans(pose);

		// Traverse each residue span, extract fragment atom coordinates and score.
		for ( ResidueSpan const & residue_span : valid_residue_spans ) {
			// Short-circuit if the span is shorter than the minimum fragment length.
			if ( residue_span.second - residue_span.first < store_->fragment_specification.fragment_length ) {
				continue;
			}

			std::vector< numeric::xyzVector<numeric::Real> > query_coordinates;
			for ( Size i = residue_span.first; i < residue_span.second; i++ ) {
				for ( std::string const & atom_name : store_->fragment_specification.fragment_atoms ) {
					query_coordinates.push_back(pose.residue(i).xyz(atom_name));
				}
			}

			for ( Size i = 0; residue_span.first + i + (store_->fragment_specification.fragment_length - 1) < residue_span.second; i++ ) {
				*result_out = lookup_closest_fragment(
					&query_coordinates[i * store_->fragment_specification.fragment_atoms.size()]);
				result_out++;

				*fragment_start_out = residue_span.first + i;
				fragment_start_out++;
			}
		}
	}


	FragmentStoreCOP store() { return store_; }
	FragmentSpecification const & fragment_specification() { return store_->fragment_specification; }

protected:

	typedef std::pair<Size, Size> ResidueSpan;
	// @brief Extract valid residue spans for the fragment lookup from the given source pose.
	//
	// Returns [start, end) residue number spans on the target pose for lookup. This handles
	// chain-breaks, non-protein residues, and symmetry.
	//
	// Ex: A pose with a virtual root and two chains ([VRT] AAAA AAAA) will produce:
	//         [(2,6), (6, 10)]
	std::vector< ResidueSpan > get_fragment_residue_spans(core::pose::Pose const & target_pose);

	FragmentLookupResult lookup_fragment(std::vector< numeric::xyzVector<numeric::Real> > & query_coordinates);

	FragmentLookupResult lookup_closest_fragment(std::vector< numeric::xyzVector<numeric::Real> > & query_coordinates);

	FragmentLookupResult lookup_closest_fragment_subset(std::vector< numeric::xyzVector<numeric::Real> > & query_coordinates, std::vector<bool> subset);

	std::vector<FragmentLookupResult> lookup_close_fragments(std::vector< numeric::xyzVector<numeric::Real> > & query_coordinates,  Real rms_threshold, Size max_matches=99999999);

private:
	FragmentStoreOP store_;
	numeric::coordinate_fitting::CoordinateArray_RMSD_FlatLookup<numeric::Real> lookup_;
};

}
}
#endif
