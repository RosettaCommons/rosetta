// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ./src/protocols/indexed_structure_store/DirectSegmentLookup.cc
/// @brief
/// @details
/// @author Alex Ford (fordas@uw.edu)

#include <vector>
#include <numeric>

#include <boost/range/combine.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range.hpp>
#include <boost/format.hpp>

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
#include <protocols/indexed_structure_store/orient_array.hh>
#include <protocols/indexed_structure_store/vector_tools.hh>
#include <protocols/indexed_structure_store/pose_utility.hh>

#include <core/kinematics/MoveMap.hh>

#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/indexed_structure_store/DirectSegmentLookup.hh>
#include <protocols/indexed_structure_store/DirectSegmentLookup.json.hh>

#include "ndarray.h"
#include "ndarray/eigen.h"

static basic::Tracer TR("protocols.indexed_structure_store.DirectSegmentLookup");

namespace ndarray {

/**
*  @brief Create a view into an array with interior contiguous dimensions merged.
*
*  The first template parameter sets the dimension of the output array and must
*  be specified directly. Only row-major contiguous dimensions can be flattened.
*
*  The second template parameter specified the number of trailing dimensions to preserve.
*/
template <int Nf, int Np, typename T, int N, int C>
inline typename boost::enable_if_c< ((C+Nf-N-Np)>=1), ArrayRef<T,Nf,(C+Nf-N)> >::type
flatten_center(Array<T,N,C> const & input) {
	typedef detail::ArrayAccess< ArrayRef<T,Nf,(C+Nf-N)> > Access;
	typedef typename Access::Core Core;
	BOOST_STATIC_ASSERT(C+Nf-N-Np >= 1);

	Vector<Size,N> oldShape = input.getShape();
	Vector<Offset,N> oldStrides = input.getStrides();
	Vector<Size,Nf> newShape;
	Vector<Offset,Nf> newStrides;

	for ( int i = 0; i < Nf - Np - 1; ++i ) {
		newShape[i] = oldShape[i];
		newStrides[i] = oldStrides[i];
	}

	int nfi = Nf - Np - 1;
	newShape[nfi] = 1;
	newStrides[nfi] = oldStrides[N - Np - 1];

	for ( int fi = Nf - Np - 1; fi < N - Np; ++fi ) {
		newShape[nfi] *= oldShape[fi];
	}

	for ( int k = 0; k < Np; ++k ) {
		newShape[Nf - k - 1] = oldShape[N - k - 1];
		newStrides[Nf - k - 1] = oldStrides[N - k - 1];
	}

	return Access::construct(input.getData(), Core::create(newShape, newStrides, input.getManager()));
}
}

namespace protocols { namespace indexed_structure_store {

using namespace protocols::indexed_structure_store;
using namespace protocols::indexed_structure_store::search;

ndarray::Array<ResidueEntry, 1, 1>
_pose_res(
	core::pose::Pose & pose,
	core::Size start_res, core::Size end_res
) {
	runtime_assert(start_res < end_res);
	runtime_assert(end_res <= pose.total_residue() + 1);

	ndarray::Array<ResidueEntry, 1, 1> res(end_res - start_res);
	for ( core::Size i = 0; i < end_res - start_res; ++i ) {
		res[i] = extract_residue_entry(pose.residue(start_res + i));
	}

	return res;
}


std::vector<DirectSegmentLookupResult> DirectSegmentLookup::segment_lookup(
	ndarray::Array<ResidueEntry, 1> source_residues,
	StructureDatabase & structure_db,
	core::pose::Pose & context,
	core::Size n_start_res, core::Size n_end_res,
	core::Size c_start_res, core::Size c_end_res
) {
	TR << json(config) << std::endl;
	TR << boost::str(boost::format(
		"n_start_res: %i n_end_res: %i c_start_res %i c_end_res %i")
		% n_start_res % n_end_res % c_start_res % c_end_res
		) << std::endl;

	ndarray::Array<float, 3, 1> source_orient = orient_array(source_residues);

	ndarray::Array<ResidueEntry, 1, 1> n_res = _pose_res(context, n_start_res, n_end_res);
	ndarray::Array<float, 3, 3> n_orient = ndarray::copy(orient_array(n_res));
	core::Size n_context = n_end_res - n_start_res;

	ndarray::Array<ResidueEntry, 1, 1> c_res = _pose_res(context, c_start_res, c_end_res);
	ndarray::Array<float, 3, 3> c_orient = ndarray::copy(orient_array(c_res));
	core::Size c_context = c_end_res - c_start_res;

	TR <<"n_context: " << n_context;
	TR << std::endl;
	TR.Debug <<"n_context res:";
	for ( auto r : n_res ) {
		TR.Debug << "\n" << json(r).dump(2);
	}
	TR.Debug << std::endl;

	TR <<"c_context: " << c_context;
	TR << std::endl;
	TR.Debug <<"c_context res:";
	for ( auto r : c_res ) {
		TR.Debug << "\n" << json(r).dump(2);
	}
	TR.Debug << std::endl;

	utility::vector1<int> context_map(context.total_residue(), 0);
	std::iota(context_map.begin(), context_map.end(), 1);
	for ( core::Size i = n_start_res; i < n_end_res; ++i ) {
		context_map[i] = 0;
	}
	for ( core::Size i = c_start_res; i < c_end_res; ++i ) {
		context_map[i] = 0;
	}
	core::pose::xyzStripeHashPose context_hash(context, context_map, core::pose::PoseCoordPickMode_N_CA_C_CB);

	StructurePairQuery query(
		n_orient, c_orient,
		config.rmsd_tolerance,
		0, config.max_insertion_length + n_context
	);

	PairQueryExecutor executor(query);
	executor.execute(structure_db);

	query_stats = executor.query_stats;
	TR << json(query_stats) << std::endl;


	typedef PairQueryExecutor::QueryResult QResult;
	boost::sort(
		executor.query_results,
		[](const QResult & a, const QResult & b) {
			return a.result_rmsd < b.result_rmsd;
		});
	if ( executor.query_results.size() > 0 ) {
		TR << "min_rmsd:" << executor.query_results[0].result_rmsd << std::endl;
	} else {
		TR << "no results" << std::endl;
	}

	std::map<int, std::vector<QResult>> query_results_by_length;
	for ( auto & r : executor.query_results ) {
		query_results_by_length[r.fragment_b_start - r.fragment_a_start + c_context].push_back(r);
	}

	std::vector<DirectSegmentLookupResult> results;
	for ( auto & r : query_results_by_length ) {
		core::Size segment_length = r.first;
		auto & len_results = r.second;

		TR << boost::str(boost::format("grouping segment length: %s count: %s") % segment_length % len_results.size()) << std::endl;
		// Load query, lookup result endpoint and lookup result segment coordinates into buffer.
		ndarray::Array<float, 4, 4> query_coords(1, n_context + c_context, 4, 3);
		query_coords[0][ndarray::view(0, n_context)] = n_orient;
		query_coords[0][ndarray::view(n_context, n_context + c_context)] = c_orient;

		ndarray::Array<ResidueEntry, 2, 2> segment_residues(len_results.size(), segment_length);
		for ( core::Size i = 0; i < len_results.size(); ++i ) {
			segment_residues[i] = source_residues[ndarray::view(
				len_results[i].fragment_a_start, len_results[i].fragment_b_start + c_context)];
		}

		ndarray::Array<float, 4, 4> segment_coords(len_results.size(), segment_length, 4, 3);
		segment_coords.deep() = orient_array(segment_residues);

		ndarray::Array<float, 4, 4> endpoint_coords(len_results.size(), n_context + c_context, 4, 3);
		endpoint_coords[ndarray::view()(0, n_context)] =
			segment_coords[ndarray::view()(0, n_context)];
		endpoint_coords[ndarray::view()(n_context, n_context + c_context)] =
			segment_coords[ndarray::view()(segment_length - c_context, segment_length)];

		// Setup flattened views of coordinate buffers w/o per-atom dimension
		ndarray::Array<float, 3, 1> f_query_coords = ndarray::flatten_center<3, 1>(query_coords);
		ndarray::Array<float, 3, 1> f_endpoint_coords = ndarray::flatten_center<3, 1>(endpoint_coords);
		ndarray::Array<float, 3, 1> f_segment_coords = ndarray::flatten_center<3, 1>(segment_coords);
		ndarray::Array<float, 1> rmsd_out(len_results.size());

		// Align lookup results onto query endpoints
		numeric::alignment::coordinate_array_superimpose(
			f_endpoint_coords,
			f_query_coords,
			f_segment_coords,
			rmsd_out);

		orient_array(segment_residues) = segment_coords;

		// Iterate in rmsd-order though result set, pruning to result clusters
		std::vector<std::vector<core::Size>> cluster_indicies;
		float cluster_tolerance = pow(config.segment_cluster_tolerance, 2);
		int chainbreak_prune_count = 0;
		int clash_prune_count = 0;
		int cluster_prune_count = 0;

		for ( core::Size ri = 0; ri < f_segment_coords.getSize<0>(); ++ri ) {
			// Perform chainbreak check
			int chain_ending = 0;
			for ( core::Size residue_i = 0; residue_i < segment_residues.getSize<1>(); ++ residue_i ) {
				if ( segment_residues(ri, residue_i).chain_ending ) {
					chain_ending += 1;
				}
			}

			if ( chain_ending > 0 ) {
				TR.Trace << "Pruned segment for chainbreak." << std::endl;
				chainbreak_prune_count += 1;
				continue;
			}

			// Perform context clash check
			int clash = 0;
			for ( core::Size ai = 1; ai < f_segment_coords.getSize<1>() - 1; ++ai ) {
				numeric::geometry::hashing::Ball ball;
				ball.x() = f_segment_coords(ri, ai, 0);
				ball.y() = f_segment_coords(ri, ai, 1);
				ball.z() = f_segment_coords(ri, ai, 2);
				ball.radius(1.5);

				if ( context_hash.clash_check_ball(ball) > 0 ) {
					clash += 1;
				}
			}

			if ( clash > 3 ) {
				TR.Trace << "Pruned segment for clashes: " << clash << std::endl;
				clash_prune_count += 1;
				continue;
			}

			// Perform segment cluster identity check
			bool unique = true;
			for ( auto & cluster : cluster_indicies ) {

				auto coord_msd = (f_segment_coords[cluster.front()].asEigen() - f_segment_coords[ri].asEigen()).rowwise().squaredNorm().mean();
				if ( coord_msd < cluster_tolerance ) {
					unique = false;
					cluster.push_back(ri);
					break;
				}
			}

			if ( !unique ) {
				TR.Trace << "Pruned segment to cluster." << std::endl;
				cluster_prune_count += 1;
				continue;
			}
			cluster_indicies.push_back({ri});
		}

		TR << boost::format("pruned to unique segments: %i chainbreaks: %i clashed: %i clustered: %i")
			% cluster_indicies.size() % chainbreak_prune_count % clash_prune_count % cluster_prune_count << std::endl;

		for ( auto & cluster : cluster_indicies ) {
			DirectSegmentLookupResult result;
			result.result_residues.assign(
				segment_residues[cluster.front()].begin(),
				segment_residues[cluster.front()].end()
			);

			for ( auto & member_index : cluster ) {
				result.query_results.push_back( len_results.at(member_index) );
			}

			results.push_back(result);
		}
	}

	return results;
}

} }
