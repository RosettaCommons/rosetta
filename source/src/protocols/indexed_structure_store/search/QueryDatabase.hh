// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ./src/protocols/indexed_structure_store/search/QueryDatabase.hh
/// @brief Database for generic alignment-based queries over structure stores.
/// @author Alex Ford (fordas@uw.edu)
//
#pragma once

#include <iostream>

#include <cstdlib>
#include <cmath>
#include <vector>
#include <map>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "ndarray.h"

#include <protocols/indexed_structure_store/search/QueryDatabase.fwd.hh>
#include <protocols/indexed_structure_store/Datatypes.fwd.hh>

namespace protocols { namespace indexed_structure_store { namespace search {

typedef int64_t Index;
typedef float SearchReal;

class StructureData
{
public:

	typedef Eigen::Matrix<SearchReal, 3, 1> Coordinate;
	typedef Eigen::Matrix<SearchReal, 3, Eigen::Dynamic> CoordinateMatrix;
	typedef Eigen::Array<Index, Eigen::Dynamic, 1> IndexArray;

	Index structure_offset;
	CoordinateMatrix coordinate_buffer;
	Index n_entry;
	Index c_per_entry;
	IndexArray chain_endpoints;

	std::map<Index, IndexArray > fragment_index_cache;
	std::map<Index, CoordinateMatrix > fragment_center_of_mass_cache;

	//std::map<Index, SpatialIndex > fragment_center_of_mass_index_cache;

	void initialize(
		Index src_structure_offset,
		ndarray::Array<SearchReal, 3, 3> src_coordinates,
		IndexArray src_chain_endpoints
	);

	void prepare_for_query(Index fragment_length);

	//void prepare_spatial_index( Index fragment_length );

	IndexArray & get_fragment_indicies(Index fragment_length)
	{
		prepare_for_query(fragment_length);
		return fragment_index_cache[fragment_length];
	}

	CoordinateMatrix & get_fragment_centers_of_mass(Index fragment_length)
	{
		prepare_for_query(fragment_length);
		return fragment_center_of_mass_cache[fragment_length];
	}

	//SpatialIndex & get_spatial_index(Index fragment_length)
	//{
	//prepare_spatial_index( fragment_length );

	//return fragment_center_of_mass_index_cache[fragment_length];
	//}

	IndexArray get_fragment_indicies_copy(Index fragment_length)
	{
		prepare_for_query(fragment_length);
		return fragment_index_cache[fragment_length];
	}

	CoordinateMatrix get_fragment_centers_of_mass_copy(Index fragment_length)
	{
		prepare_for_query(fragment_length);
		return fragment_center_of_mass_cache[fragment_length];
	}
};

class StructureDatabase
{
public:
	typedef std::map<Index, StructureData> StructureDataMap;

	void initialize(ndarray::Array<ResidueEntry, 1> residue_entries);

	void initialize(
		ndarray::Array<SearchReal, 3, 3> coordinate_buffer,
		ndarray::Array<Index, 1> structure_endpoints,
		ndarray::Array<Index, 1> chain_endpoints
	);

	StructureData & get_structure_data(Index structure_id){
		return structure_data.at(structure_id);
	}

	StructureDataMap structure_data;

private:
	void initialize_structure_data_from_buffer(
		StructureData & structure_data,
		ndarray::Array<SearchReal, 3, 3> coordinate_buffer,
		ndarray::Array<Index, 1> structure_endpoints,
		ndarray::Array<Index, 1> chain_endpoints,
		Index & structure_endpoint_index,
		Index & chain_endpoint_index,
		Index & current_structure_start
	);
};

class StructurePairQuery
{
	typedef Eigen::Matrix<SearchReal, 3, 1> Coordinate;
	typedef Eigen::Matrix<SearchReal, 3, Eigen::Dynamic> CoordinateMatrix;

public:

	StructurePairQuery(
		ndarray::Array<SearchReal, 3, 3> query_coordinates_a,
		ndarray::Array<SearchReal, 3, 3> query_coordinates_b,
		SearchReal rmsd_tolerance_
	);

	StructurePairQuery(
		ndarray::Array<SearchReal, 3, 3> query_coordinates_a,
		ndarray::Array<SearchReal, 3, 3> query_coordinates_b,
		SearchReal rmsd_tolerance_,
		Index min_primary_distance_,
		Index max_primary_distance_
	);

	void init(
		ndarray::Array<SearchReal, 3, 3> query_coordinates_a,
		ndarray::Array<SearchReal, 3, 3> query_coordinates_b
	);

	Index n_entry_a;
	Index n_entry_b;

	CoordinateMatrix q_buffer_a;
	CoordinateMatrix q_buffer_b;

	Index c_per_entry;

	SearchReal rmsd_tolerance;

	bool limit_primary_distance;
	Index min_primary_distance;
	Index max_primary_distance;
};

class StructureSingleQuery
{
	typedef Eigen::Matrix<SearchReal, 3, 1> Coordinate;
	typedef Eigen::Matrix<SearchReal, 3, Eigen::Dynamic> CoordinateMatrix;

public:

	StructureSingleQuery(
		ndarray::Array<SearchReal, 3, 3> query_coordinates,
		SearchReal rmsd_tolerance_
	);

	Index n_entry;
	Index c_per_entry;

	CoordinateMatrix q_buffer;

	SearchReal rmsd_tolerance;
};

struct StructurePairQueryResult
{
	Index fragment_a_start;
	Index fragment_b_start;
	Index structure_index;
	Index fragment_a_structure_start;
	Index fragment_b_structure_start;
	SearchReal result_rmsd;
};

struct StructureSingleQueryResult
{
	Index fragment_start;
	Index structure_index;
	Index fragment_structure_start;
	SearchReal result_rmsd;
};

struct PairQuerySummaryStatistics
{
	Index structures_considered = 0;
	Index fragments_considered = 0;
	Index fragments_expanded = 0;
	Index pairs_considered = 0;
	Index pairs_aligned = 0;
	Index result_count = 0;
};

struct SingleQuerySummaryStatistics
{
	Index structures_considered = 0;
	Index fragments_considered = 0;
	Index result_count = 0;
};

class PairQueryExecutor
{
public:
	typedef PairQuerySummaryStatistics SummaryStatistics;
	typedef StructurePairQueryResult QueryResult;
	typedef StructurePairQuery Query;

	PairQueryExecutor(Query const & query);

	void execute(StructureDatabase & database);
	void execute_structure( Index structure_index, StructureData & target_structure );

	SearchReal perform_structure_rmsd(
		StructureData & target_structure,
		StructureData::IndexArray & fragment_indicies_a,
		StructureData::IndexArray & fragment_indicies_b,
		StructureData::CoordinateMatrix & fragment_centers_of_mass_a,
		StructureData::CoordinateMatrix & fragment_centers_of_mass_b,
		Index fragment_index_a,
		Index fragment_index_b
	);

	Query query;

	Eigen::Matrix<SearchReal, 3, Eigen::Dynamic> query_coordinate_buffer;
	Eigen::Matrix<SearchReal, 3, Eigen::Dynamic> structure_coordinate_buffer;

	Eigen::Matrix<SearchReal, 3, 1> query_coordinate_com;
	Eigen::Matrix<SearchReal, 3, 1> structure_coordinate_com;

	std::vector<QueryResult> query_results;

	SummaryStatistics query_stats;
};

class SingleQueryExecutor
{
public:
	typedef SingleQuerySummaryStatistics SummaryStatistics;
	typedef StructureSingleQueryResult QueryResult;
	typedef StructureSingleQuery Query;

	SingleQueryExecutor(Query const & query);

	void execute(StructureDatabase & database);

	void execute_structure( Index structure_index, StructureData & target_structure );

	SearchReal perform_structure_rmsd(
		StructureData & target_structure,
		StructureData::IndexArray & fragment_indicies,
		StructureData::CoordinateMatrix & fragment_centers_of_mass,
		Index fragment_index
	);

	Query query;

	Eigen::Matrix<SearchReal, 3, Eigen::Dynamic> query_coordinate_buffer;
	Eigen::Matrix<SearchReal, 3, Eigen::Dynamic> structure_coordinate_buffer;

	Eigen::Matrix<SearchReal, 3, 1> query_coordinate_com;
	Eigen::Matrix<SearchReal, 3, 1> structure_coordinate_com;

	std::vector<QueryResult> query_results;

	SummaryStatistics query_stats;
};

}}}
