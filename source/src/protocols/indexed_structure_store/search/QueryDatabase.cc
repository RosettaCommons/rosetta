// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ./src/protocols/indexed_structure_store/search/QueryDatabase.cc
/// @brief Database for generic alignment-based queries over structure stores.
/// @author Alex Ford (fordas@uw.edu)
//
#include <iostream>

#include <cstdlib>
#include <cmath>
#include <vector>
#include <map>

#include <boost/foreach.hpp>

#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "ndarray.h"
#include "ndarray/eigen.h"

#include <numeric/alignment/QCPKernel.hh>

#include <protocols/indexed_structure_store/search/QueryDatabase.hh>
#include <protocols/indexed_structure_store/Datatypes.hh>
#include <protocols/indexed_structure_store/orient_array.hh>
#include <protocols/indexed_structure_store/vector_tools.hh>

namespace protocols { namespace indexed_structure_store { namespace search {

using namespace protocols::indexed_structure_store;


void StructureData::initialize(
	Index src_structure_offset,
	ndarray::Array<SearchReal, 3, 3> src_coordinates,
	IndexArray src_chain_endpoints
)
{
	structure_offset = src_structure_offset;

	if ( src_coordinates.template getSize<2>() != 3 ) {
		throw std::invalid_argument("src_coordinates shape[2] != 3.");
	}
	n_entry = src_coordinates.template getSize<0>();
	c_per_entry = src_coordinates.template getSize<1>();

	coordinate_buffer.conservativeResize(3, n_entry * c_per_entry);
	for ( Index i = 0; i < n_entry; ++i ) {
		coordinate_buffer.block(0, i * c_per_entry, 3, c_per_entry) =
			src_coordinates[i].template asEigen<Eigen::MatrixXpr>().transpose();
	}

	chain_endpoints = src_chain_endpoints;
	// Insert additional chain endpoint at structure endpoint if needed.
	if ( (chain_endpoints.size() == 0) || (chain_endpoints[chain_endpoints.size() - 1] != n_entry - 1) ) {
		chain_endpoints.conservativeResize(chain_endpoints.size() + 1);
		chain_endpoints[chain_endpoints.size() - 1] = n_entry - 1;
	}
}

void StructureData::prepare_for_query(Index fragment_length)
{
	if ( fragment_index_cache.count(fragment_length) == 0 ) {
		IndexArray & fragment_indicies = fragment_index_cache[fragment_length];

		Index chain_index = 0;
		Index chain_start = 0;
		Index chain_end = 0;

		while ( chain_index < chain_endpoints.size() ) {
			chain_end = chain_endpoints[chain_index];

			int num_fragments_in_chain = ((chain_end + 1) - chain_start) - fragment_length + 1;

			if ( num_fragments_in_chain > 0 ) {
				fragment_indicies.conservativeResize(fragment_indicies.size() + num_fragments_in_chain);
				fragment_indicies.bottomRightCorner(num_fragments_in_chain, 1) = IndexArray::LinSpaced(num_fragments_in_chain, 0, num_fragments_in_chain - 1) + chain_start;
			}

			chain_index += 1;
			chain_start = chain_end + 1;
		}
	}

	if ( fragment_center_of_mass_cache.count(fragment_length) == 0 ) {
		IndexArray & fragment_indicies = fragment_index_cache[fragment_length];
		CoordinateMatrix & fragment_center_of_mass = fragment_center_of_mass_cache[fragment_length];
		fragment_center_of_mass.resize(3, fragment_indicies.size());

		for ( Index fragment_index = 0; fragment_index < fragment_indicies.size(); fragment_index++ ) {
			fragment_center_of_mass.col(fragment_index) =
				coordinate_buffer.block(0, fragment_indicies[fragment_index] * c_per_entry, 3, fragment_length * c_per_entry )
				.rowwise().sum() / (fragment_length * c_per_entry);
		}
	}
}

//void StructureData::prepare_spatial_index( Index fragment_length )
//{
//prepare_for_query(fragment_length);

//if(fragment_center_of_mass_index_cache.count(fragment_length) == 0) {
//CoordinateMatrix & fragment_centers_of_mass = fragment_center_of_mass_cache[fragment_length];

//std::vector<SpatialIndexEntry> index_entries;
//for( Index i = 0; i < fragment_centers_of_mass.cols(); i++) {
//index_entries.push_back( SpatialIndexEntry(fragment_centers_of_mass.col(i), i));
//}

//fragment_center_of_mass_index_cache[fragment_length] = SpatialIndex(
//index_entries.begin(),
//index_entries.end());
//}
//};

void StructureDatabase::initialize(ndarray::Array<ResidueEntry, 1> residue_entries) {
	ndarray::Array<SearchReal, 3, 3> coordinates(residue_entries.getSize<0>(), 4, 3);
	coordinates.deep() = orient_array(residue_entries);

	std::vector<Index> structure_endpoints;
	std::vector<Index> chain_endpoints;

	for ( ndarray::Size i=0; i < residue_entries.getSize<0>() - 1; ++i ) {
		if ( residue_entries[i].structure_id != residue_entries[i+1].structure_id ) {
			structure_endpoints.push_back(i);
			chain_endpoints.push_back(i);
		} else if ( residue_entries[i].chain_ending ) {
			chain_endpoints.push_back(i);
		}
	}

	initialize( coordinates, v_to_a(structure_endpoints), v_to_a(chain_endpoints) );
}

void StructureDatabase::initialize(
	ndarray::Array<SearchReal, 3, 3> coordinate_buffer,
	ndarray::Array<Index, 1> structure_endpoints,
	ndarray::Array<Index, 1> chain_endpoints
)
{
	if ( !structure_data.empty() ) {
		throw std::invalid_argument("StructureDatabase already initialized.");
	}

	if ( coordinate_buffer.template getSize<2>() != 3 ) {
		throw std::invalid_argument("coordinate_buffer shape[2] != 3.");
	}

	Index current_structure_id = 0;
	Index current_structure_start = 0;
	Index current_structure_endpoint_index = 0;
	Index current_chain_endpoint_index = 0;

	while ( current_structure_endpoint_index < (Index)structure_endpoints.template getSize<0>() )
			{
		initialize_structure_data_from_buffer(
			structure_data[current_structure_id],
			coordinate_buffer,
			structure_endpoints,
			chain_endpoints,
			current_structure_endpoint_index,
			current_chain_endpoint_index,
			current_structure_start
		);

		current_structure_id++;

	}

	if ( current_structure_start < (Index)coordinate_buffer.template getSize<0>() ) {
		initialize_structure_data_from_buffer(
			structure_data[current_structure_id],
			coordinate_buffer,
			structure_endpoints,
			chain_endpoints,
			current_structure_endpoint_index,
			current_chain_endpoint_index,
			current_structure_start
		);
	}
}

void StructureDatabase::initialize_structure_data_from_buffer(
	StructureData & struct_data,
	ndarray::Array<SearchReal, 3, 3> coordinate_buffer,
	ndarray::Array<Index, 1> structure_endpoints,
	ndarray::Array<Index, 1> chain_endpoints,
	Index & structure_endpoint_index,
	Index & chain_endpoint_index,
	Index & current_structure_start
)
{
	Index coordinate_span_start = current_structure_start;
	Index coordinate_span_end;
	if ( structure_endpoint_index < (Index)structure_endpoints.getSize<0>() ) {
		coordinate_span_end = structure_endpoints[structure_endpoint_index] + 1;

		if ( coordinate_span_end > (Index)coordinate_buffer.getSize<0>() ) {
			throw std::invalid_argument("structure_endpoints value exceeded coordinate buffer size.");
		}

	} else {
		coordinate_span_end = coordinate_buffer.template getSize<0>();
	}

	Index chain_endpoints_start = chain_endpoint_index;
	Index chain_endpoints_end = chain_endpoint_index;
	while ( (chain_endpoints_end < (Index)chain_endpoints.getSize<0>()) && (chain_endpoints[chain_endpoints_end] < coordinate_span_end) ) {
		chain_endpoints_end++;
	}

	struct_data.initialize(
		coordinate_span_start,
		coordinate_buffer[ndarray::view(coordinate_span_start, coordinate_span_end)],
		chain_endpoints[ndarray::view(chain_endpoints_start, chain_endpoints_end)].template asEigen<Eigen::ArrayXpr>() - coordinate_span_start
	);

	structure_endpoint_index += 1;
	chain_endpoint_index = chain_endpoints_end;
	current_structure_start = coordinate_span_end;
}

StructurePairQuery::StructurePairQuery(
	ndarray::Array<SearchReal, 3, 3> query_coordinates_a,
	ndarray::Array<SearchReal, 3, 3> query_coordinates_b,
	SearchReal rmsd_tolerance_
) :
	n_entry_a(),
	n_entry_b(),
	q_buffer_a(),
	q_buffer_b(),
	rmsd_tolerance( rmsd_tolerance_ ),
	limit_primary_distance(false),
	min_primary_distance(),
	max_primary_distance()
{
	init(query_coordinates_a, query_coordinates_b);
}

StructurePairQuery::StructurePairQuery(
	ndarray::Array<SearchReal, 3, 3> query_coordinates_a,
	ndarray::Array<SearchReal, 3, 3> query_coordinates_b,
	SearchReal rmsd_tolerance_,
	Index min_primary_distance_,
	Index max_primary_distance_
) :
	n_entry_a(),
	n_entry_b(),
	q_buffer_a(),
	q_buffer_b(),
	rmsd_tolerance( rmsd_tolerance_ ),
	limit_primary_distance(true),
	min_primary_distance(min_primary_distance_),
	max_primary_distance(max_primary_distance_)
{
	init(query_coordinates_a, query_coordinates_b);
}

void StructurePairQuery::init(
	ndarray::Array<SearchReal, 3, 3> query_coordinates_a,
	ndarray::Array<SearchReal, 3, 3> query_coordinates_b
) {
	// Validate query shapes
	if ( query_coordinates_a.template getSize<2>() != 3 ) {
		throw std::invalid_argument("query_coordinates_a shape[2] != 3.");
	}

	if ( query_coordinates_b.template getSize<2>() != 3 ) {
		throw std::invalid_argument("query_coordinates_b shape[2] != 3.");
	}

	if ( query_coordinates_a.template getSize<1>() != query_coordinates_b.template getSize<1>() ) {
		throw std::invalid_argument("query_coordinates_a shape[1] != query_coordinates_b shape[1]");

	}

	// Unpack input arrays into coordinate buffers
	c_per_entry = query_coordinates_a.template getSize<1>();

	n_entry_a = query_coordinates_a.template getSize<0>();
	q_buffer_a.resize(3, n_entry_a * c_per_entry);
	for ( Index i = 0; i < n_entry_a; ++i ) {
		q_buffer_a.block(0, i * c_per_entry, 3, c_per_entry) =
			query_coordinates_a[i].template asEigen<Eigen::MatrixXpr>().transpose();
	}

	n_entry_b = query_coordinates_b.template getSize<0>();
	q_buffer_b.resize(3, n_entry_b * c_per_entry);
	for ( Index i = 0; i < n_entry_b; ++i ) {
		q_buffer_b.block(0, i * c_per_entry, 3, c_per_entry) =
			query_coordinates_b[i].template asEigen<Eigen::MatrixXpr>().transpose();
	}
}

StructureSingleQuery::StructureSingleQuery(
	ndarray::Array<SearchReal, 3, 3> query_coordinates,
	SearchReal rmsd_tolerance_
) :
	n_entry(),
	c_per_entry(),
	q_buffer(),
	rmsd_tolerance( rmsd_tolerance_ )
{
	// Validate query shapes
	if ( query_coordinates.template getSize<2>() != 3 ) {
		throw std::invalid_argument("query_coordinates shape[2] != 3.");
	}

	// Unpack input arrays into coordinate buffers
	c_per_entry = query_coordinates.template getSize<1>();

	n_entry = query_coordinates.template getSize<0>();
	q_buffer.resize(3, n_entry * c_per_entry);
	for ( Index i = 0; i < n_entry; ++i ) {
		q_buffer.block(0, i * c_per_entry, 3, c_per_entry) =
			query_coordinates[i].template asEigen<Eigen::MatrixXpr>().transpose();
	}
}

PairQueryExecutor::PairQueryExecutor(Query const & q) : query(q)
{
	query_coordinate_buffer.resize(3, q.q_buffer_a.cols() + q.q_buffer_b.cols());

	query_coordinate_buffer.leftCols(q.q_buffer_a.cols()) = q.q_buffer_a;
	query_coordinate_buffer.rightCols(q.q_buffer_b.cols()) = q.q_buffer_b;
	query_coordinate_com = query_coordinate_buffer.rowwise().sum() / query_coordinate_buffer.cols();

	structure_coordinate_buffer.resize(3, query_coordinate_buffer.cols());
}

void PairQueryExecutor::execute(StructureDatabase & database)
{
	for ( auto & structure_data_pair : database.structure_data ) {
		execute_structure( structure_data_pair.first, structure_data_pair.second);
	}

	query_stats.result_count = query_results.size();
}

void PairQueryExecutor::execute_structure(Index structure_index, StructureData & target_structure )
{
	StructureData::IndexArray & fragment_indicies_a = target_structure.get_fragment_indicies(query.n_entry_a);
	StructureData::CoordinateMatrix & fragment_centers_of_mass_a = target_structure.get_fragment_centers_of_mass(query.n_entry_a);

	StructureData::IndexArray & fragment_indicies_b = target_structure.get_fragment_indicies(query.n_entry_b);
	StructureData::CoordinateMatrix & fragment_centers_of_mass_b = target_structure.get_fragment_centers_of_mass(query.n_entry_b);

	query_stats.structures_considered++;

	for ( Index fragment_index_a = 0; fragment_index_a < fragment_indicies_a.size(); fragment_index_a++ ) {
		query_stats.fragments_considered++;
		query_stats.fragments_expanded++;

		for ( Index fragment_index_b = 0; fragment_index_b < fragment_indicies_b.size(); fragment_index_b++ ) {
			Index a = fragment_indicies_a[fragment_index_a];
			Index b = fragment_indicies_b[fragment_index_b];

			// Skip b if fragment overlaps with fragment a
			if ( (b >= a) & (b < a + query.n_entry_a) ) {
				continue;
			} else if ( ((b + query.n_entry_b >= a) & (b + query.n_entry_b < a + query.n_entry_a)) ) {
				continue;
			} else if ( query.limit_primary_distance ) {
				Index primary_distance = b - a;
				if ( primary_distance < query.min_primary_distance || primary_distance > query.max_primary_distance ) {
					continue;
				}
			}

			query_stats.pairs_considered++;

			SearchReal result_rmsd = this->perform_structure_rmsd(
				target_structure,
				fragment_indicies_a,
				fragment_indicies_b,
				fragment_centers_of_mass_a,
				fragment_centers_of_mass_b,
				fragment_index_a,
				fragment_index_b);

			if ( result_rmsd < query.rmsd_tolerance ) {
				query_results.push_back(
					QueryResult{
					a + target_structure.structure_offset,
					b + target_structure.structure_offset,
					structure_index, a, b,
					result_rmsd,
					});
			}
		}
	}
}

SearchReal PairQueryExecutor::perform_structure_rmsd(
	StructureData & target_structure,
	StructureData::IndexArray & fragment_indicies_a,
	StructureData::IndexArray & fragment_indicies_b,
	StructureData::CoordinateMatrix & fragment_centers_of_mass_a,
	StructureData::CoordinateMatrix & fragment_centers_of_mass_b,
	Index fragment_index_a,
	Index fragment_index_b)
{
	query_stats.pairs_aligned++;
	Index a = fragment_indicies_a[fragment_index_a];
	Index b = fragment_indicies_b[fragment_index_b];

	structure_coordinate_buffer.leftCols(query.q_buffer_a.cols()) =
		target_structure.coordinate_buffer.block(0, a * query.c_per_entry, 3, query.q_buffer_a.cols());
	structure_coordinate_buffer.rightCols(query.q_buffer_b.cols()) =
		target_structure.coordinate_buffer.block(0, b * query.c_per_entry, 3, query.q_buffer_b.cols());

	structure_coordinate_com = (
		(fragment_centers_of_mass_a.col(fragment_index_a) * query.q_buffer_a.cols()) +
		(fragment_centers_of_mass_b.col(fragment_index_b) * query.q_buffer_b.cols()))
		/ (query.q_buffer_a.cols() + query.q_buffer_b.cols());

	return numeric::alignment::QCPKernel<SearchReal>::calc_coordinate_rmsd( query_coordinate_buffer, query_coordinate_com, structure_coordinate_buffer, structure_coordinate_com);
}

SingleQueryExecutor::SingleQueryExecutor(Query const & q) : query(q)
{
	query_coordinate_buffer = q.q_buffer;
	query_coordinate_com = query_coordinate_buffer.rowwise().sum() / query_coordinate_buffer.cols();

	structure_coordinate_buffer.resize(3, query_coordinate_buffer.cols());
}

void SingleQueryExecutor::execute(StructureDatabase & database)
{
	for ( auto & structure_data_pair : database.structure_data ) {
		execute_structure( structure_data_pair.first, structure_data_pair.second) ;
	}

	query_stats.result_count = query_results.size();
}

void SingleQueryExecutor::execute_structure(Index structure_index, StructureData & target_structure )
{
	StructureData::IndexArray & fragment_indicies = target_structure.get_fragment_indicies(query.n_entry);
	StructureData::CoordinateMatrix & fragment_centers_of_mass = target_structure.get_fragment_centers_of_mass(query.n_entry);

	query_stats.structures_considered++;

	for ( Index fragment_index = 0; fragment_index < fragment_indicies.size(); fragment_index++ ) {
		query_stats.fragments_considered++;
		Index a = fragment_indicies[fragment_index];
		SearchReal result_rmsd = this->perform_structure_rmsd(
			target_structure,
			fragment_indicies,
			fragment_centers_of_mass,
			fragment_index);

		if ( result_rmsd < query.rmsd_tolerance ) {
			query_results.push_back(
				QueryResult{
				a + target_structure.structure_offset,
				structure_index, a,
				result_rmsd
				});
		}
	}
}

SearchReal SingleQueryExecutor::perform_structure_rmsd(
	StructureData & target_structure,
	StructureData::IndexArray & fragment_indicies,
	StructureData::CoordinateMatrix & fragment_centers_of_mass,
	Index fragment_index)
{
	Index a = fragment_indicies[fragment_index];

	structure_coordinate_buffer = target_structure.coordinate_buffer.block(0, a * query.c_per_entry, 3, query.q_buffer.cols());
	structure_coordinate_com = fragment_centers_of_mass.col(fragment_index);

	return numeric::alignment::QCPKernel<SearchReal>::calc_coordinate_rmsd(
		query_coordinate_buffer,
		query_coordinate_com,
		structure_coordinate_buffer,
		structure_coordinate_com
	);
}

} } }
