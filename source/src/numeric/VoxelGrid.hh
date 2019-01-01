// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    numeric/VoxelGrid.hh
/// @brief   Class that provides optimal space hashing for rapid identification of neighbors
/// @details Provides retrieval of objects in 3D space by space hashing (1,2, or 3D Voxel grid
///          depending on object density. Object references and coordinates are stored on the
///          grid and allow for very rapid retrieval of neighbors within a distance <= the
///          resolution of the grid.
/// @author  Jeff Mendenhall (jeffrey.l.mendenhall@vanderbilt.edu)
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)
///

#ifndef INCLUDED_numeric_VoxelGrid_hh
#define INCLUDED_numeric_VoxelGrid_hh

// Unit headers
#include <numeric/VoxelGrid.fwd.hh>

// Core headers
#include <numeric/geometry/BoundingBox.hh>
#include <platform/types.hh>

// Package headers
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>

namespace numeric {

template< typename T >
class VoxelGrid
{

public: // Typedefs

	typedef T t_DataType;
	typedef platform::Size Size;
	typedef platform::Real Real;
	typedef xyzVector<Real> Vector;
	typedef std::pair< T const *, T const * > t_DataTypeCOPtrPair;
	typedef std::pair< T const *, Vector const * > t_DataTypeVectorCOPtrPair;

private:

	enum BorderFlag
	{
		e_None = 0,
		e_XMin = 1,
		e_XMax = 2,
		e_YMin = 4,
		e_YMax = 8,
		e_ZMin = 16,
		e_ZMax = 32
	};

public: // Methods

	// Construction

	/// @brief default constructor
	VoxelGrid(
		Real const & resolution = 4.0,
		bool const & cache_edges = false
	);

	/// @brief destructor
	virtual
	~VoxelGrid();

	// Operations

	/// @brief extract the 3D coordinates of a given t_DataType input. TO BE IMPLEMENTED BY DERIVED CLASSES.
	/// @param input: reference to object T
	virtual
	Vector const *
	ExtractPosition(t_DataType const & input) const = 0;

	/// @brief check if two grid items are the same. TO BE IMPLEMENTED BY DERIVED CLASSES.
	/// @param item1: first comparison item
	/// @param item2: second comparison item
	virtual
	bool
	IsSameItem(t_DataType const & item1, t_DataType const & item2) const = 0;

	/// @brief is item relevant to consider as neighbor. TO BE IMPLEMENTED BY DERIVED CLASSES.
	/// @param item: item to check
	virtual
	bool
	IsRelevantItem(t_DataType const & item) const = 0;

	/// @brief Get the number of items in the voxel grid
	Size
	GetNumberItems() const;

	/// @brief Get the actual internally used dimension of the voxel grid
	Size
	GetDimension() const;

	/// @brief Get the requested resolution by the user.
	Real
	GetResolution() const;

	/// @brief Get the bounding box of this grid
	geometry::BoundingBox< Vector > const &
	GetBoundingBox() const;

	/// @brief Get the items of this grid
	utility::vector1< utility::vector1< t_DataTypeVectorCOPtrPair > > const &
	GetGridItems() const;

	/// @brief Find neighbors within a specified distance of a given DATA
	/// @param input: datapoint of interest
	/// @param neighborhood: range in which to search for neighbors; must be <= resolution
	/// @param check_relevance: check and only return neighboring objects which are relevant
	utility::vector1< std::pair< t_DataType const *, Real > >
	GetNeighbors(
		t_DataType const & input,
		Real const & neighborhood,
		bool check_relevance = false
	) const;

	/// @brief Has given DATA a neighbor in the grid within a specified distance
	/// @param input: datapoint of interest
	/// @param neighborhood: range in which to search for neighbors; must be <= resolution
	/// @param check_relevance: check and only return true if neighboring objects are relevant
	bool
	HasNeighbors(
		t_DataType const & input,
		Real const & neighborhood,
		bool check_relevance = false
	) const;

	/// @brief Find all neighbor pairs within a specified distance
	/// @param grid: other grid to search for neighbors with this grid
	/// @param neighborhood: range in which to search for neighbors; must be <= resolution
	/// @param check_relevance: check and only return neighboring objects which are relevant
	utility::vector1< std::pair< t_DataTypeCOPtrPair, Real > >
	GetNeighborsIn(
		VoxelGrid const & grid,
		Real const & neighborhood,
		bool check_relevance = false
	) const;

	/// @brief Find all neighbor pairs within a specified distance
	/// @param neighborhood: range in which to search for neighbors; must be <= resolution
	/// @param check_relevance: check and only return neighboring objects which are relevant
	utility::vector1< std::pair< t_DataTypeCOPtrPair, Real> >
	GetNeighbors(
		Real const & neighborhood,
		bool check_relevance = false
	) const;

	/// @brief Update the VoxelGrid with new data
	/// @param new_data: Vector of the new data
	/// @param check_relevance: check if objects are relevant for neighbor calculation and only insert those
	void
	SetObjects(
		utility::vector1< t_DataType const * > const & new_data,
		bool check_relevance = false
	);

	/// @brief insert a new t_DataType object into the manager
	/// @param new_item: reference to item to insert
	/// @param check_relevance: check and only insert if object is relevant for neighbor calculation
	void
	InsertObject(
		t_DataType const & new_item,
		bool check_relevance = false
	);

	/// @brief inserts multiple t_DataType objects into the manager
	/// @param new_items: vector of references to items to insert
	/// @param check_relevance: check and only insert if objects are relevant for neighbor calculation
	void
	InsertObjects(
		utility::vector1< t_DataType const * > const & new_items,
		bool check_relevance = false
	);

	/// @brief remove a t_DataType object from the manager
	/// @param item_to_remove: reference to the item to remove
	void
	RemoveObject(t_DataType const & item_to_remove);

	/// @brief remove multiple t_DataType objects from the manager
	/// @param items_to_remove: vector of references to items to remove
	void
	RemoveObjects(utility::vector1< t_DataType const * > const & items_to_remove);

	/// @brief Remove all objects from the voxel grid while keeping the grid itself intact
	void
	Clear();

	/// @brief Translate the grid
	/// @param translation: amount to translate by
	void
	Translate(Vector const & translation);

private: // Methods

	/// @brief Get the position for the given coordinate vector.
	///        Places out-of-bounds vector in the bin closest to where they belong.
	/// @param X, Y, Z the position of interest
	Size
	GetIndex(
		Real const & X,
		Real const & Y,
		Real const & Z
	) const;

	/// @brief Get the position for the given coordinate vector.
	///        Places out-of-bounds vector in the bin closest to where they belong.
	/// @param pos: the position of interest
	Size
	GetIndex(Vector const & pos) const;

	/// @brief Try to allocate a grid of a given size. If this size would
	///        cause the object to definitely exceed the max memory allotment,
	///        the grid will be truncated to include only the core. It is assumed
	///        that Resolution_ has been setup prior to this.
	/// @param min_vals: dimension of bounding box
	/// @param max_vals: dimension of bounding box
	void
	SetupGrid(
		Vector const & min_vals,
		Vector const & max_vals,
		Size const & expected_n_elements
	);

	/// @brief Compute distance from this point to the nearest edge of the bounding box
	///        of all given points, provided the point is outside the box.
	/// @param point: point of interest
	Real
	GetSqDistanceOutsideBoundingBox(Vector const & point) const;

	/// @brief Compute distance from another grid's bounding box to the nearest edge of this
	///        bounding box provided that the MIN and MAX is outside of this bounding box.
	/// @param min_box: min of the grid of interest
	/// @param min_box: max of the grid of interest
	Real
	GetSqDistanceOutsideBoundingBox(
		Vector const & min_box,
		Vector const & max_box
	) const;

	/// @brief Find all neighbor pairs of this grid within a specified distance
	/// @param neighborhood: range in which to search for neighbors; must be <= resolution
	/// @param check_relevance: check and only return neighboring objects which are relevant
	utility::vector1< std::pair< t_DataTypeCOPtrPair, Real > >
	GetNeighborsMultiDimensional(
		Real const & neighborhood,
		bool check_relevance = false
	) const;

	/// @brief Find all neighbors of this object within a specified distance
	/// @param obj: object to look for neighbors
	/// @param neighborhood: range in which to search for neighbors; must be <= resolution
	/// @param check_relevance: check and only return neighboring objects which are relevant
	utility::vector1< std::pair< t_DataType const *, Real > >
	GetNeighborsMultiDimensional(
		t_DataType const & obj,
		Real const & neighborhood,
		bool check_relevance = false
	) const;

	/// @brief Find all neighbors between the input and this grid within a specified distance
	/// @param grid: input grid in which to look for objects that are neighbors of this grid
	/// @param neighborhood: range in which to search for neighbors; must be <= resolution
	/// @param check_relevance: check and only return neighboring objects which are relevant
	utility::vector1< std::pair< t_DataTypeCOPtrPair, Real > >
	GetNeighborsMultiDimensional(
		VoxelGrid< t_DataType > const & grid,
		Real const & neighborhood,
		bool check_relevance = false
	) const;

	/// @brief Has given object a neighbor in the grid within a specified distance
	/// @param obj: object to look for neighbors
	/// @param neighborhood: range in which to search for neighbors; must be <= resolution
	/// @param check_relevance: check and only return true if neighboring objects are relevant
	bool
	HasNeighborsMultiDimensional(
		t_DataType const & obj,
		Real const & neighborhood,
		bool check_relevance = false
	) const;

	/// @brief Find all neighbor pairs of this 1D grid within a specified distance
	/// @param neighborhood: range in which to search for neighbors; must be <= resolution
	/// @param check_relevance: check and only return neighboring objects which are relevant
	utility::vector1< std::pair< t_DataTypeCOPtrPair, Real > >
	GetNeighbors1D(
		Real const & neighborhood,
		bool check_relevance = false
	) const;

	/// @brief Find all neighbors of this object within a specified distance
	/// @param obj: object to look for neighbors
	/// @param neighborhood: range in which to search for neighbors; must be <= resolution
	/// @param check_relevance: check and only return neighboring objects which are relevant
	utility::vector1< std::pair< t_DataType const *, Real > >
	GetNeighbors1D(
		t_DataType const & obj,
		Real const & neighborhood,
		bool check_relevance = false
	) const;

	/// @brief Find all neighbors between the input and this grid within a specified distance
	/// @param grid: input grid in which to look for objects that are neighbors of this grid
	/// @param neighborhood: range in which to search for neighbors; must be <= resolution
	/// @param check_relevance: check and only return neighboring objects which are relevant
	utility::vector1< std::pair< t_DataTypeCOPtrPair, Real > >
	GetNeighbors1D(
		VoxelGrid< t_DataType> const & grid,
		Real const & neighborhood,
		bool check_relevance = false
	) const;

	/// @brief Has given object a neighbor in the grid within a specified distance
	/// @param obj: object to look for neighbors
	/// @param neighborhood: range in which to search for neighbors; must be <= resolution
	/// @param check_relevance: check and only return true if neighboring objects are relevant
	bool
	HasNeighbors1D(
		t_DataType const & obj,
		Real const & neighborhood,
		bool check_relevance = false
	) const;

private: // Data

	// Object assignments
	utility::vector1< utility::vector1< t_DataTypeVectorCOPtrPair> > Assignments_;

	// Neighboring voxels; used only for 2-3D grids, ignored for 1D.
	// This set is for directed edges; used when detected internal clashes
	utility::vector1< utility::vector1< utility::vector1< t_DataTypeVectorCOPtrPair> > > Edges_;

	// Neighboring voxels; used only for 2-3D grids, ignored for 1D
	// This set is for undirected edges; used when detected external clashes
	utility::vector1< utility::vector1< utility::vector1< t_DataTypeVectorCOPtrPair> > > EdgesComplete_;

	// Border flags, used if not caching edges
	utility::vector1< Size > BorderFlags_;

	// Number of bins in each direction
	Size NBinsX_;
	Size NBinsY_;
	Size NBinsZ_;

	// Requested resolution by the user
	Real Resolution_;

	// Should be true if the same grid will be used repeatedly (~20% speed improvement)
	// Edges should not be cached if the grid is only going to be used once, otherwise a
	// significant speed penalty can be expected
	bool CacheEdges_;

	// Resolution of the optimized grid in each direction
	Real ResolutionX_;
	Real ResolutionY_;
	Real ResolutionZ_;

	// Minimum resolution in any direction
	Real MinResolution_;

	// Grid origin
	Real MinX_;
	Real MinY_;
	Real MinZ_;

	// Bounding box of grid
	geometry::BoundingBox< Vector > Box_;

	// Number of dimensions expanded by the voxel grid
	Size NDimensional_;

	// Number of items
	Size NItems_;

	// Minimum number of elements this object should have; if smaller, reallocate the grid
	Size MinNumberElements_;
	// Maximum number of elements this object should have; if larger, reallocate the grid
	Size MaxNumberElements_;


};

} // namespace numeric

#endif // INCLUDED_numeric_VoxelGrid_hh
