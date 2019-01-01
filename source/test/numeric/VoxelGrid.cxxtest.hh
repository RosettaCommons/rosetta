// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/VoxelGrid.cxxtest.hh
/// @brief  test suite for numeric::VoxelGrid
/// @author Georg Kuenze (georg.kuenze@vanderbilt.edu)


// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <numeric/VoxelGrid.impl.hh>

// Package headers
#include <numeric/xyzVector.hh>
#include <utility/vector1.hh>

// C++ headers
#include <iostream>
#include <algorithm>

class SimplePoint {

public: // Methods

	inline
	SimplePoint() :
		id_(),
		coords_()
	{}

	inline
	SimplePoint(
		size_t const & id, numeric::xyzVector<double> & coords
	) :
		id_(id),
		coords_(coords)
	{}

	inline
	SimplePoint(
		size_t const && id, numeric::xyzVector<double> && coords
	) :
		id_(id),
		coords_(coords)
	{}

	inline
	~SimplePoint()
	{}

	inline
	size_t const &
	id() const
	{
		return id_;
	}

	inline
	numeric::xyzVector<double> const &
	coords() const
	{
		return coords_;
	}

	inline
	bool
	is_relevant() const
	{
		return true;
	}

private: // Data
	size_t id_;
	numeric::xyzVector<double> coords_;

};

class VoxelGridSimplePoints : public numeric::VoxelGrid< SimplePoint > {

public:

	inline
	VoxelGridSimplePoints(
		double const & resolution,
		utility::vector1< SimplePoint > const & points,
		bool const & cache_edges = false
	) :
		numeric::VoxelGrid< SimplePoint >(resolution, cache_edges)
	{
		utility::vector1< SimplePoint const * > points_ptr;
		points_ptr.reserve(points.size());
		for ( auto & p : points ) {
			points_ptr.push_back(&p);
		}
		SetObjects(points_ptr);
	}

	inline
	VoxelGridSimplePoints(
		double const & resolution,
		utility::vector1< SimplePoint const * > const & points,
		bool const & cache_edges = false
	) :
		numeric::VoxelGrid< SimplePoint >(resolution, cache_edges)
	{
		SetObjects(points);
	}

	inline
	~VoxelGridSimplePoints() override
	{}

	inline
	numeric::xyzVector<double> const *
	ExtractPosition(
		SimplePoint const & point
	) const override
	{
		return &(point.coords());
	}

	inline
	bool
	IsSameItem(
		SimplePoint const & point1,
		SimplePoint const & point2
	) const override
	{
		return point1.id() == point2.id();
	}

	inline
	bool
	IsRelevantItem(
		SimplePoint const & point
	) const override
	{
		return point.is_relevant();
	}

};

class VoxelGridTests : public CxxTest::TestSuite {

public: // Typedefs

	typedef VoxelGridSimplePoints VGSPs;
	typedef numeric::xyzVector<double> Vector3D;
	typedef std::pair< SimplePoint const *, double > PointDistPair;
	typedef std::pair< std::pair< SimplePoint const *, SimplePoint const * >, double > PointPointDistTriple;

public:

	// Shared initialization goes here.
	void setUp() {
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	// Test the creation of a 3D grid
	void test_grid_creation() {
		using namespace numeric;
		using namespace utility;

		vector1< SimplePoint > points;
		points.reserve(8000);

		size_t count(0);
		for ( size_t i(0); i < size_t(10); ++i ) {
			for ( size_t j(0); j < size_t(10); ++j ) {
				for ( size_t k(0); k < size_t(10); ++k ) {
					points.emplace_back(SimplePoint(count++, Vector3D(-4.6+i, -4.6+j, -4.4+k)));
					points.emplace_back(SimplePoint(count++, Vector3D(-4.4+i, -4.6+j, -4.4+k)));
					points.emplace_back(SimplePoint(count++, Vector3D(-4.6+i, -4.4+j, -4.4+k)));
					points.emplace_back(SimplePoint(count++, Vector3D(-4.4+i, -4.4+j, -4.4+k)));
					points.emplace_back(SimplePoint(count++, Vector3D(-4.6+i, -4.6+j, -4.6+k)));
					points.emplace_back(SimplePoint(count++, Vector3D(-4.4+i, -4.6+j, -4.6+k)));
					points.emplace_back(SimplePoint(count++, Vector3D(-4.6+i, -4.4+j, -4.6+k)));
					points.emplace_back(SimplePoint(count++, Vector3D(-4.4+i, -4.4+j, -4.6+k)));
				}
			}
		}

		TS_TRACE("Creating 3D VoxelGrid with 8000 members.");
		VGSPs grid(1.0, points);
		TS_ASSERT_EQUALS(grid.GetNumberItems(), 8000);
		TS_ASSERT_EQUALS(grid.GetDimension(), 3);
		TS_ASSERT_DELTA(grid.GetBoundingBox().lower().x(), -4.6, 1.0e-6);
		TS_ASSERT_DELTA(grid.GetBoundingBox().lower().y(), -4.6, 1.0e-6);
		TS_ASSERT_DELTA(grid.GetBoundingBox().lower().z(), -4.6, 1.0e-6);
		TS_ASSERT_DELTA(grid.GetBoundingBox().upper().x(),  4.6, 1.0e-6);
		TS_ASSERT_DELTA(grid.GetBoundingBox().upper().y(),  4.6, 1.0e-6);
		TS_ASSERT_DELTA(grid.GetBoundingBox().upper().z(),  4.6, 1.0e-6);
	}

	// Test the retrieval of neighbors from a multidimensional grid
	void test_neighbor_detection_multidimensional() {
		using namespace numeric;
		using namespace utility;

		vector1< SimplePoint > points1;
		points1.reserve(8000);

		size_t count(0);
		for ( size_t i(0); i < size_t(10); ++i ) {
			for ( size_t j(0); j < size_t(10); ++j ) {
				for ( size_t k(0); k < size_t(10); ++k ) {
					points1.emplace_back(SimplePoint(count++, Vector3D(-4.6+i, -4.6+j, -4.4+k)));
					points1.emplace_back(SimplePoint(count++, Vector3D(-4.4+i, -4.6+j, -4.4+k)));
					points1.emplace_back(SimplePoint(count++, Vector3D(-4.6+i, -4.4+j, -4.4+k)));
					points1.emplace_back(SimplePoint(count++, Vector3D(-4.4+i, -4.4+j, -4.4+k)));
					points1.emplace_back(SimplePoint(count++, Vector3D(-4.6+i, -4.6+j, -4.6+k)));
					points1.emplace_back(SimplePoint(count++, Vector3D(-4.4+i, -4.6+j, -4.6+k)));
					points1.emplace_back(SimplePoint(count++, Vector3D(-4.6+i, -4.4+j, -4.6+k)));
					points1.emplace_back(SimplePoint(count++, Vector3D(-4.4+i, -4.4+j, -4.6+k)));
				}
			}
		}

		// Setup the grid
		VGSPs grid1(1.0, points1);

		// Test for a point in the center of the grid
		SimplePoint p1(8000, Vector3D(0.5, 0.5, 0.5));
		TS_TRACE("Retrieve neighbors of point (0.5, 0.5, 0.5) within 1.0 neighborhood radius.");
		vector1< PointDistPair > neighbor_vec1 = grid1.GetNeighbors(p1, 1.0);
		TS_ASSERT_EQUALS(neighbor_vec1.size(), 32);
		auto comparefunc1 = [](PointDistPair const & a, PointDistPair const & b) { return a.second < b.second; };
		double max_d1 = (*std::max_element(neighbor_vec1.begin(), neighbor_vec1.end(), comparefunc1)).second;
		TS_ASSERT_LESS_THAN_EQUALS(max_d1, 1.0);

		// Test for a point at one of the corners of the grid
		SimplePoint p2(8001, Vector3D(4.5, 4.5, 4.5));
		TS_TRACE("Retrieve neighbors of point (4.5, 4.5, 4.5) within 1.0 neighborhood radius.");
		vector1< PointDistPair > neighbor_vec2 = grid1.GetNeighbors(p2, 1.0);
		TS_ASSERT_EQUALS(neighbor_vec2.size(), 20);
		double max_d2 = (*std::max_element(neighbor_vec2.begin(), neighbor_vec2.end(), comparefunc1)).second;
		TS_ASSERT_LESS_THAN_EQUALS(max_d2, 1.0);

		// Test for a point at one of the faces of the grid
		SimplePoint p3(8003, Vector3D(4.5, 0.5, 0.5));
		TS_TRACE("Retrieve neighbors of point (4.5, 0.5, 0.5) within 1.0 neighborhood radius.");
		vector1< PointDistPair > neighbor_vec3 = grid1.GetNeighbors(p3, 1.0);
		TS_ASSERT_EQUALS(neighbor_vec3.size(), 28);
		double max_d3 = (*std::max_element(neighbor_vec3.begin(), neighbor_vec3.end(), comparefunc1)).second;
		TS_ASSERT_LESS_THAN_EQUALS(max_d3, 1.0);

		// Create another grid and search for neighbors between the two grids
		vector1< SimplePoint > points2;
		points2.reserve(125);
		count=0;
		for ( size_t i(0); i < size_t(5); ++i ) {
			for ( size_t j(0); j < size_t(5); ++j ) {
				for ( size_t k(0); k < size_t(5); ++k ) {
					points2.emplace_back(SimplePoint(count++, Vector3D(5.0+i, 5.0+j, 5.0+k)));
				}
			}
		}
		VGSPs grid2(1.0, points2);
		TS_TRACE("Retrieve neighbors between two voxel grids within 1.0 neighborhood radius.");
		vector1< PointPointDistTriple > neighbor_vec4 = grid1.GetNeighborsIn(grid2, 1.0);
		TS_ASSERT_EQUALS(neighbor_vec4.size(), 7);
		auto comparefunc2 = [](PointPointDistTriple const & a, PointPointDistTriple const & b ) { return a.second < b.second; };
		double max_d4 = (*std::max_element(neighbor_vec4.begin(), neighbor_vec4.end(), comparefunc2)).second;
		TS_ASSERT_LESS_THAN_EQUALS(max_d4, 1.0);
		Vector3D v(4.4, 4.4, 4.4);
		auto identityfunc = [v](PointPointDistTriple const & triple) { return triple.first.second->coords().distance_squared(v) <= 1.0e-6; };
		auto point_in_vector(std::find_if(neighbor_vec4.begin(), neighbor_vec4.end(), identityfunc) != neighbor_vec4.end());
		TS_ASSERT(!point_in_vector);

		// Retrieve all neighbor pairs of this grid
		// We use a small neighborhood radius of 0.30 and 0.20.
		// This should only return neighbors within the same voxel.
		TS_TRACE("Retrieve all neighbors of this grid within 0.30 neighborhood radius.");
		vector1< PointPointDistTriple > neighbor_vec5 = grid1.GetNeighbors(0.30);
		TS_ASSERT_EQUALS(neighbor_vec5.size(), 24000);
		double max_d5 = (*std::max_element(neighbor_vec5.begin(), neighbor_vec5.end(), comparefunc2)).second;
		TS_ASSERT_LESS_THAN_EQUALS(max_d5, 0.30);

		TS_TRACE("Retrieve all neighbors of this grid within 0.20 neighborhood radius.");
		vector1< PointPointDistTriple > neighbor_vec6 = grid1.GetNeighbors(0.2);
		TS_ASSERT_EQUALS(neighbor_vec6.size(), 12000);
		double max_d6 = (*std::max_element(neighbor_vec6.begin(), neighbor_vec6.end(), comparefunc2)).second;
		TS_ASSERT_LESS_THAN_EQUALS(max_d6, 0.2);
	}

	// Test the retrieval of neighbors from a one-dimensional grid
	void test_neighbor_detection_onedimensional() {
		using namespace numeric;
		using namespace utility;

		// Create a grid with a checkerboard pattern and 6 elements per cell
		vector1< SimplePoint > points1;
		points1.reserve(84);
		size_t count(0);
		for ( size_t i(1); i < size_t(4); ++i ) {
			for ( size_t j(1); j < size_t(4); ++j ) {
				for ( size_t k(1); k < size_t(4); ++k ) {
					if ( ( ((i+j)%2 == 0) && (k%2 != 0) ) ||
							( ((i+j)%2 != 0) && (k%2 == 0) ) ) {
						points1.emplace_back(SimplePoint(count++, Vector3D(0.4+i-1, 0.4+j-1, 0.5+k-1)));
						points1.emplace_back(SimplePoint(count++, Vector3D(0.4+i-1, 0.6+j-1, 0.5+k-1)));
						points1.emplace_back(SimplePoint(count++, Vector3D(0.6+i-1, 0.4+j-1, 0.5+k-1)));
						points1.emplace_back(SimplePoint(count++, Vector3D(0.6+i-1, 0.6+j-1, 0.5+k-1)));
						points1.emplace_back(SimplePoint(count++, Vector3D(0.5+i-1, 0.5+j-1, 0.6+k-1)));
						points1.emplace_back(SimplePoint(count++, Vector3D(0.5+i-1, 0.5+j-1, 0.4+k-1)));
					}
				}
			}
		}
		TS_TRACE("Creating 1D VoxelGrid with 84 members.");
		VGSPs grid1(1.0, points1);
		TS_ASSERT_EQUALS(grid1.GetNumberItems(), 84);
		TS_ASSERT_EQUALS(grid1.GetDimension(), 0);

		// Test for a point of the grid
		SimplePoint p1(84, Vector3D(1.5, 1.5, 1.5));
		TS_TRACE("Retrieve neighbors of point (1.5, 1.5, 1.5) within 1.0 neighborhood radius.");
		vector1< PointDistPair > neighbor_vec1 = grid1.GetNeighbors(p1, 1.0);
		TS_ASSERT_EQUALS(neighbor_vec1.size(), 10);
		auto comparefunc1 = [](PointDistPair const & a, PointDistPair const & b) { return a.second < b.second; };
		double max_d1 = (*std::max_element(neighbor_vec1.begin(), neighbor_vec1.end(), comparefunc1)).second;
		TS_ASSERT_LESS_THAN_EQUALS(max_d1, 1.0);

		// Create another grid and search for neighbors between the two grids
		vector1< SimplePoint > points2;
		points2.reserve(27);
		count=0;
		for ( size_t i(0); i < size_t(5); ++i ) {
			for ( size_t j(0); j < size_t(5); ++j ) {
				for ( size_t k(0); k < size_t(5); ++k ) {
					points2.emplace_back(SimplePoint(count++, Vector3D(1.5+i, 1.5+j, 1.5+k)));
				}
			}
		}
		VGSPs grid2(1.0, points2);
		// With a neighborhood radius of 0.2 we search only within each voxel
		TS_TRACE("Retrieve neighbors between two voxel grids within 0.2 neighborhood radius.");
		vector1< PointPointDistTriple > neighbor_vec2 = grid1.GetNeighborsIn(grid2, 0.2);
		TS_ASSERT_EQUALS(neighbor_vec2.size(), 24);
		auto comparefunc2 = [](PointPointDistTriple const & a, PointPointDistTriple const & b ) { return a.second < b.second; };
		double max_d2 = (*std::max_element(neighbor_vec2.begin(), neighbor_vec2.end(), comparefunc2)).second;
		TS_ASSERT_LESS_THAN_EQUALS(max_d2, 0.2);
		Vector3D v(2.5, 1.5, 2.5);
		auto identityfunc = [v](PointPointDistTriple const & triple) { return triple.first.second->coords().distance_squared(v) <= 1.0e-6; };
		auto point_in_vector(std::find_if(neighbor_vec2.begin(), neighbor_vec2.end(), identityfunc) != neighbor_vec2.end());
		TS_ASSERT(!point_in_vector);

		// Retrieve all neighbor pairs of this grid
		// We use a small neighborhood radius of 0.30 and 0.18.
		// This should only return neighbors within the same voxel.
		TS_TRACE("Retrieve all neighbors of this grid within 0.30 neighborhood radius.");
		vector1< PointPointDistTriple > neighbor_vec3 = grid1.GetNeighbors(0.30);
		TS_ASSERT_EQUALS(neighbor_vec3.size(), 210);
		double max_d3 = (*std::max_element(neighbor_vec3.begin(), neighbor_vec3.end(), comparefunc2)).second;
		TS_ASSERT_LESS_THAN_EQUALS(max_d3, 0.3);

		TS_TRACE("Retrieve all neighbors of this grid within 0.18 neighborhood radius.");
		vector1< PointPointDistTriple > neighbor_vec4 = grid1.GetNeighbors(0.18);
		TS_ASSERT_EQUALS(neighbor_vec4.size(), 112);
		double max_d4 = (*std::max_element(neighbor_vec4.begin(), neighbor_vec4.end(), comparefunc2)).second;
		TS_ASSERT_LESS_THAN_EQUALS(max_d4, 0.18);
	}
};
