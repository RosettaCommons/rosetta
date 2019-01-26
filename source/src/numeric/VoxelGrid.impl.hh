// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    numeric/VoxelGrid.impl.hh
/// @brief   Implementation file for VoxelGrid class
/// @author  Jeff Mendenhall (jeffrey.l.mendenhall@vanderbilt.edu)
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_numeric_VoxelGrid_impl_hh
#define INCLUDED_numeric_VoxelGrid_impl_hh

// Unit headers
#include <numeric/VoxelGrid.hh>

// Package headers
#include <utility/fixedsizearray1.hh>
#include <utility/exit.hh>

// C++ headers
#include <cmath>
#include <limits>

namespace numeric {

/// @brief default constructor
template< typename T >
VoxelGrid<T>::VoxelGrid
(
	Real const & resolution,
	bool const & cache_edges
) :
	NBinsX_( 0 ),
	NBinsY_( 0 ),
	NBinsZ_( 0 ),
	Resolution_( resolution ),
	CacheEdges_( cache_edges ),
	ResolutionX_( resolution ),
	ResolutionY_( resolution ),
	ResolutionZ_( resolution ),
	MinResolution_( resolution ),
	MinX_( 0 ),
	MinY_( 0 ),
	MinZ_( 0 ),
	NDimensional_( 0 ),
	NItems_( 0 ),
	MinNumberElements_( 0 ),
	MaxNumberElements_( 0 )
{}

/// @brief destructor
template< typename T >
VoxelGrid<T>::~VoxelGrid() {}

/// @brief Get the number of items in the voxel grid
template< typename T >
Size
VoxelGrid<T>::GetNumberItems() const
{
	return NItems_;
}

/// @brief Get the actual internally used dimension of the voxel grid
template< typename T >
Size
VoxelGrid<T>::GetDimension() const
{
	return NDimensional_;
}

/// @brief Get the requested resolution by the user.
template< typename T >
Real
VoxelGrid<T>::GetResolution() const
{
	return Resolution_;
}

/// @brief Get the bounding box of this grid
template< typename T >
geometry::BoundingBox< xyzVector<Real> > const &
VoxelGrid<T>::GetBoundingBox() const
{
	return Box_;
}

/// @brief Get the items of this grid
template< typename T >
utility::vector1< utility::vector1< std::pair< T const *, xyzVector<Real> const * > > > const &
VoxelGrid<T>::GetGridItems() const
{
	return Assignments_;
}

/// @brief Find neighbors within a specified distance of a given DATA
/// @param input: datapoint of interest
/// @param neighborhood: range in which to search for neighbors; must be <= resolution
/// @param check_relevance: check and only return neighboring objects which are relevant
template< typename T >
utility::vector1< std::pair< T const *, Real > >
VoxelGrid<T>::GetNeighbors(
	t_DataType const & input,
	Real const & neighborhood,
	bool check_relevance
) const
{
	return NDimensional_ <= Size( 1 ) ? GetNeighbors1D( input, neighborhood, check_relevance ) : GetNeighborsMultiDimensional( input, neighborhood, check_relevance );
}

/// @brief Has given DATA a neighbor in the grid within a specified distance
/// @param input: datapoint of interest
/// @param neighborhood: range in which to search for neighbors; must be <= resolution
/// @param check_relevance: check and only return true if neighboring objects are relevant
template< typename T >
bool
VoxelGrid<T>::HasNeighbors(
	t_DataType const & input,
	Real const & neighborhood,
	bool check_relevance
) const
{
	return NDimensional_ <= Size( 1 ) ? HasNeighbors1D( input, neighborhood, check_relevance ) : HasNeighborsMultiDimensional( input, neighborhood, check_relevance );
}

/// @brief Find all neighbor pairs within a specified distance
/// @param grid: other grid to search for neighbors with this grid
/// @param neighborhood: range in which to search for neighbors; must be <= resolution
/// @param check_relevance: check and only return neighboring objects which are relevant
template< typename T >
utility::vector1< std::pair< std::pair< T const *, T const * >, Real > >
VoxelGrid<T>::GetNeighborsIn(
	VoxelGrid< t_DataType > const & grid,
	Real const & neighborhood,
	bool check_relevance
) const
{
	return NDimensional_ <= Size( 1 ) ? GetNeighbors1D( grid, neighborhood, check_relevance ) : GetNeighborsMultiDimensional( grid, neighborhood, check_relevance );
}

/// @brief Find all neighbor pairs within a specified distance
/// @param neighborhood: range in which to search for neighbors; must be <= resolution
/// @param check_relevance: check and only return neighboring objects which are relevant
template< typename T >
utility::vector1< std::pair< std::pair< T const *, T const * >, Real> >
VoxelGrid<T>::GetNeighbors(
	Real const & neighborhood,
	bool check_relevance
) const
{
	return NDimensional_ <= Size( 1 ) ? GetNeighbors1D( neighborhood, check_relevance ) : GetNeighborsMultiDimensional( neighborhood, check_relevance );
}

/// @brief Update the VoxelGrid with new data
/// @param new_data: Vector of the new data
/// @param check_relevance: check if objects are relevant for neighbor calculation and only insert those
template< typename T >
void
VoxelGrid<T>::SetObjects(
	utility::vector1< t_DataType const * > const & new_data,
	bool check_relevance
)
{
	// track the boundaries of the cube enclosing all points given
	Box_.set_lower(Vector(std::numeric_limits< Real >::max()));
	Box_.set_upper(Vector(std::numeric_limits< Real >::min()));

	// store references to the objects and their coordinates
	utility::vector1< t_DataTypeVectorCOPtrPair > ref_vec;
	ref_vec.reserve( new_data.size() );
	for ( auto itr( new_data.begin() ), itr_end( new_data.end() ); itr != itr_end; ++itr ) {
		if ( !check_relevance || ( check_relevance && IsRelevantItem( **itr ) ) ) {
			Vector const * coord( ExtractPosition( **itr ) );
			if ( coord ) {
				Box_.add( *coord );
				ref_vec.push_back( t_DataTypeVectorCOPtrPair( *itr, coord ) );
			}
		}
	}
	// Determine optimal dimensionality for the grid and set it up
	SetupGrid( Box_.lower(), Box_.upper(), ref_vec.size() );

	// populate the grid
	for ( auto itr( ref_vec.begin() ), itr_end( ref_vec.end() ); itr != itr_end; ++itr ) {
		const Size index( GetIndex( *(itr->second) ) );
		Assignments_[ index ].push_back( *itr );
	}
	NItems_ = ref_vec.size();
}

/// @brief insert a new t_DataType object into the manager
/// @param new_item: reference to item to insert
/// @param check_relevance: check and only insert if object is relevant for neighbor calculation
template< typename T >
void
VoxelGrid<T>::InsertObject(
	t_DataType const & new_item,
	bool check_relevance
)
{
	if ( !check_relevance || ( check_relevance && IsRelevantItem( new_item ) ) ) {
		Vector const * coord( ExtractPosition( new_item ) );
		if ( coord ) {
			const Size index( GetIndex( *coord ) );
			Assignments_[ index ].push_back( t_DataTypeVectorCOPtrPair( &new_item, coord ) );
			++NItems_;
			Box_.add( *coord );
		}
	}
}

/// @brief inserts multiple t_DataType objects into the manager
/// @param new_items: vector of references to items to insert
/// @param check_relevance: check and only insert if objects are relevant for neighbor calculation
template< typename T >
void
VoxelGrid<T>::InsertObjects(
	utility::vector1< t_DataType const * > const & new_items,
	bool check_relevance
)
{
	for ( Size i = 1, sz = new_items.size(); i <= sz; ++i ) {
		InsertObject( *(new_items[i]), check_relevance );
	}
}

/// @brief remove a t_DataType object from the manager
/// @param item_to_remove: reference to the item to remove
template< typename T >
void
VoxelGrid<T>::RemoveObject(t_DataType const & item_to_remove)
{
	Vector const * pos( ExtractPosition( item_to_remove ));
	if ( pos && NBinsX_ ) {
		const Size index( GetIndex( *pos ) );
		utility::vector1< t_DataTypeVectorCOPtrPair > & vec( Assignments_[ index ] );

		for ( auto itr( vec.begin() ), itr_end( vec.end() ); itr !=itr_end; ++itr ) {
			if ( IsSameItem( *(itr->first), item_to_remove ) ) {
				vec.erase( itr );
				--NItems_;
				break;
			}
		}
	}
}

/// @brief remove multiple t_DataType objects from the manager
/// @param items_to_remove: vector of references to items to remove
template< typename T >
void
VoxelGrid<T>::RemoveObjects(utility::vector1< t_DataType const * > const & items_to_remove)
{
	for ( Size i = 1, sz = items_to_remove.size(); i <= sz; ++i ) {
		RemoveObject( *(items_to_remove[i]) );
	}
}

/// @brief Remove all objects from the voxel grid while keeping the grid itself intact
template< typename T >
void
VoxelGrid<T>::Clear()
{
	for ( Size i = 1, sz = Assignments_.size(); i <= sz; ++i ) {
		Assignments_[ i ].clear();
	}
	NItems_ = 0;
}

/// @brief Translate the grid
/// @param translation: amount to translate by
template< typename T >
void
VoxelGrid<T>::Translate(Vector const & translation)
{
	if ( NItems_ ) {
		MinX_ += translation.x();
		MinY_ += translation.y();
		MinZ_ += translation.z();
		Vector new_min( Box_.lower());
		new_min += translation;
		Vector new_max( Box_.upper());
		new_max += translation;
		Box_.set_lower(new_min);
		Box_.set_upper(new_max);
	}
}

/// @brief Get the position for the given coordinate vector.
///        Places out-of-bounds vector in the bin closest to where they belong.
/// @param X, Y, Z the position of interest
template< typename T >
Size
VoxelGrid<T>::GetIndex
(
	Real const & X,
	Real const & Y,
	Real const & Z
) const
{
	return Size( std::max( std::min( (X - MinX_) / ResolutionX_, NBinsX_ - 1.0 ), 0.0 ) )     // x bin
		+ NBinsX_ *
		(   Size( std::max( std::min( (Y - MinY_) / ResolutionY_, NBinsY_ - 1.0 ), 0.0) )      // y bin
		+ Size( std::max( std::min( (Z - MinZ_) / ResolutionZ_, NBinsZ_ - 1.0 ), 0.0) ) * NBinsY_  // z bin
		)
		+ 1;                       // 1-based indexing
}

/// @brief Get the position for the given coordinate vector.
///        Places out-of-bounds vector in the bin closest to where they belong.
/// @param pos: the position of interest
template< typename T >
Size
VoxelGrid<T>::GetIndex(Vector const & pos) const
{
	return GetIndex( pos.x(), pos.y(), pos.z() );
}

/// @brief Try to allocate a grid of a given size. If this size would
///        cause the object to definitely exceed the max memory allotment,
///        the grid will be truncated to include only the core. It is assumed
///        that Resolution_ has been setup prior to this.
/// @param min_vals: dimension of bounding box
/// @param max_vals: dimension of bounding box
template< typename T >
void
VoxelGrid<T>::SetupGrid
(
	Vector const & min_vals,
	Vector const & max_vals,
	Size const & expected_n_elements
)
{
	// always clear the grid of any objects
	Clear();

	if ( !expected_n_elements ) {
		return;
	}

	runtime_assert_msg(
		min_vals.x() <= max_vals.x() && min_vals.y() <= max_vals.y() &&
		min_vals.z() <= max_vals.z(), "Error in setting up VoxelGrid. Min vals and max vals incorrect."
	);

	// check whether we can continue using the existing grid
	if (
			min_vals.x() >= MinX_ && min_vals.y() >= MinY_ && min_vals.z() >= MinZ_
			&& max_vals.x() <= MinX_ + NBinsX_ * ResolutionX_
			&& max_vals.y() <= MinY_ + NBinsY_ * ResolutionY_
			&& max_vals.z() <= MinZ_ + NBinsZ_ * ResolutionZ_
			&& expected_n_elements >= MinNumberElements_
			&& expected_n_elements <= MaxNumberElements_
			) {
		return;
	}

	Vector diff( max_vals - min_vals);

	// Overallocation factor (>= 1.0). Grid dimensions overallocated by this amount to reduce allocation frequency
	static const Real s_OverAlloc( CacheEdges_ ? 1.1 : 1.0);
	ResolutionX_ = ResolutionY_ = ResolutionZ_ = Resolution_;

	NBinsX_ = std::max( Size( ( diff.x() + ResolutionX_) * s_OverAlloc / ResolutionX_), Size(1));
	NBinsY_ = std::max( Size( ( diff.y() + ResolutionY_) * s_OverAlloc / ResolutionY_), Size(1));
	NBinsZ_ = std::max( Size( ( diff.z() + ResolutionZ_) * s_OverAlloc / ResolutionZ_), Size(1));

	// two bins are always less efficient than 1 - results in the same number of comparisons but a lot more iteration
	// In numerous tests, 3 bins is also virtually always much slower owing to the additional iterations and fact that
	// the central bin will still be compared with bins on both sides, so only minimal comparisons are avoided
	if ( NBinsX_ <= 3 ) {
		NBinsX_ = 1;
		ResolutionX_ = std::numeric_limits< Real >::max() / 3.0;
	}
	if ( NBinsY_ <= 3 ) {
		NBinsY_ = 1;
		ResolutionY_ = std::numeric_limits< Real >::max() / 3.0;
	}
	if ( NBinsZ_ <= 3 ) {
		NBinsZ_ = 1;
		ResolutionZ_ = std::numeric_limits< Real >::max() / 3.0;
	}

	// Next, decide whether this should be 1D, 2D, or 3D Voxel grid
	// Generally, the goal is to keep the expected objects per bin and bins per object as close to 1 as possible
	// There is also some bonus to keeping the dimensions smaller since it means vastly fewer bin iterations
	// Empirically, the benefit is roughly 15 to go from 2D - 1D (owing to a space-oversampling optimization that is
	// fast at 1D), and another 50 to go from 3D -> 2D. In other words, 3D clustering only makes sense when we have, on average,
	// more than about 65 elements total in each 1D grid of equivalent dimension.
	const Real elements( expected_n_elements );
	utility::fixedsizearray1<Real, 7> dimension_bins;
	utility::fixedsizearray1<Real, 7> dimension_score;
	utility::fixedsizearray1<Real, 4> dimension_penalty;
	dimension_penalty[1] = 0;
	dimension_penalty[2] = 0;
	dimension_penalty[3] = 15;
	dimension_penalty[4] = 65;

	Real best_score( elements + 2);
	int best_pos(0);
	for ( int pos = 1; pos <= 7; ++pos ) {
		const int use_x( !!( pos & 1)), use_y( !!( pos & 2)), use_z( !!( pos & 4));
		const int dimension( use_x + use_y + use_z);
		const Real bins( ( use_x ? NBinsX_ : 1) * ( use_y ? NBinsY_ : 1) * ( use_z ? NBinsZ_ : 1) );
		dimension_bins[pos] = bins;
		// score is elements per bin + bins per element + dimension penalty
		const Real score( elements / bins + bins / elements + dimension_penalty[ dimension + 1 ]);
		dimension_score[pos] = score;
		if ( score < best_score ) {
			best_pos = pos;
			best_score = score;
		}
	}

	if ( !( best_pos & 1) ) {
		NBinsX_ = 1;
		ResolutionX_ = std::numeric_limits< Real >::max() / 3.0;
	}
	if ( !( best_pos & 2) ) {
		NBinsY_ = 1;
		ResolutionY_ = std::numeric_limits< Real >::max() / 3.0;
	}
	if ( !( best_pos & 4) ) {
		NBinsZ_ = 1;
		ResolutionZ_ = std::numeric_limits< Real >::max() / 3.0;
	}
	NDimensional_ = ( (NBinsX_ > 1) + (NBinsY_ > 1) + (NBinsZ_ > 1) );

	// while it can be computed exactly when the score would be better for a different grid allocation, an
	// approximate method is adequate and better suited to ensuring that we avoid re-allocating the grid due to
	// changes in density as much as possible
	MinNumberElements_ = expected_n_elements / 4;
	MaxNumberElements_ = NDimensional_ < 3 ? expected_n_elements * 4 : std::numeric_limits< Size >::max();

	// oversample if grid will be in 1D and there are sufficient elements.
	// This reduces the number of distance calculations and comparisons necessary
	if ( NDimensional_ == Size(1) ) {
		// for 1D grids, it's not much of a burden to allocate for much larger grid than is absolutely required
		const Real additional_over_alloc_factor( CacheEdges_ ? 1.25 : 1.0);
		if ( NBinsX_ > Size(1) ) {
			NBinsX_ = Size( Real(NBinsX_) * additional_over_alloc_factor );
			if ( NBinsX_ * 3 < expected_n_elements ) {
				NBinsX_ *= 3;
				ResolutionX_ /= 3;
			}
		} else if ( NBinsY_ > Size(1) ) {
			//  NBinsY_ = Size( Real(NBinsY_) * additional_over_alloc_factor );
			if ( NBinsY_ * 3 < expected_n_elements ) {
				NBinsY_ *= 3;
				ResolutionY_ /= 3;
			}
		} else { // if ( NBinsZ_ > Size(1))
			//  NBinsZ_ = Size( Real(NBinsZ_) * additional_over_alloc_factor );
			if ( NBinsZ_ * 3 < expected_n_elements ) {
				NBinsZ_ *= 3;
				ResolutionZ_ /= 3;
			}
		}
	}
	const Size total_requested_bins(NBinsX_ * NBinsY_ * NBinsZ_);
	MinResolution_ = std::min( ResolutionX_, std::min( ResolutionY_, ResolutionZ_) );

	// define origin
	MinX_ = min_vals.x();
	MinY_ = min_vals.y();
	MinZ_ = min_vals.z();

	// update origin shifted (to more negative values) due to overallocation and dimensional collapsing
	MinX_ -= ( MinX_ + ResolutionX_ * NBinsX_ - max_vals.x()) * 0.5;
	MinY_ -= ( MinY_ + ResolutionY_ * NBinsY_ - max_vals.y()) * 0.5;
	MinZ_ -= ( MinZ_ + ResolutionZ_ * NBinsZ_ - max_vals.z()) * 0.5;

	Assignments_.resize( total_requested_bins);

	if ( NDimensional_ > Size(1) ) {
		// # of bins to explore in each direction. This can be easily expanded to try different binning strategies
		const int bins_to_x( 1 ), bins_to_y( 1 ), bins_to_z( 1 );

		// last xyz bins, cached for convenience
		const int last_x_bin( NBinsX_ - 1 ), last_y_bin( NBinsY_ - 1 ), last_z_bin( NBinsZ_ - 1 );

		// compute edges for 2-3 dimensional grids
		if ( CacheEdges_ ) {
			Edges_.resize( total_requested_bins);
			EdgesComplete_.resize( total_requested_bins);

			// distance in memory between adjacent bins on Z-axis
			const Size z_stride( NBinsX_ * NBinsY_ );

			// determine side of edge vectors
			const Size edge_vec_size_internal( NDimensional_ == Size(2) ? 4 : 13 );
			const Size edge_vec_size_external( NDimensional_ == Size(2) ? 8 : 26 );

			// compute edges for multi-dimensional grids
			for ( int z = 0, pos = 0, nbz = NBinsZ_; z < nbz; ++z ) {

				const int min_z( -std::min( z, bins_to_z)), max_z( std::min( z + bins_to_z, last_z_bin) - z + 1);
				for ( int y = 0, nby =  NBinsY_; y < nby; ++y ) {

					const int min_y( -std::min( y, bins_to_y)), max_y( std::min( y + bins_to_y, last_y_bin) - y + 1);
					for ( int x = 0 , nbx = NBinsX_; x < nbx; ++x, ++pos ) {

						// determine extents in each direction
						const int min_x( -std::min( x, bins_to_x)), max_x( std::min( x + bins_to_x, last_x_bin) - x + 1);

						// reset existing edges
						utility::vector1< utility::vector1< t_DataTypeVectorCOPtrPair > > & edges( Edges_[ pos + 1 ] );
						edges.clear();
						edges.resize( edge_vec_size_internal );
						utility::vector1< utility::vector1< t_DataTypeVectorCOPtrPair > > & edgesc( EdgesComplete_[ pos + 1 ] );
						edgesc.clear();
						edgesc.resize( edge_vec_size_external);
						for ( int z_off = min_z, fz_off = min_z * z_stride; z_off < max_z; ++z_off, fz_off += z_stride ) {

							for ( int y_off = min_y, zy_off = fz_off + min_y * NBinsX_; y_off < max_y; ++y_off, zy_off += NBinsX_ ) {

								for ( int x_off = min_x, zyx_off = zy_off + min_x; x_off < max_x; ++x_off, ++zyx_off ) {

									if ( zyx_off > 0 ) {
										// directed edges
										edges.push_back( Assignments_[ pos + zyx_off + 1 ] );
									}
									if ( zyx_off != 0 ) {
										// undirected edges; just ensure that it isn't the same point
										edgesc.push_back( Assignments_[ pos + zyx_off + 1 ] );
									}
								}
							}
						}
					}
				}
			}
		} else {
			BorderFlags_.resize( total_requested_bins );
			// compute flags for multidimensional grids
			for ( int z = 0, pos = 0, nbz = NBinsZ_; z < nbz; ++z ) {

				int z_flag( int( z ? e_None : e_ZMin) | int( z < last_z_bin ? e_None : e_ZMax) );
				for ( int y = 0, nby = NBinsY_; y < nby; ++y ) {

					int zy_flag( z_flag | int( y ? e_None : e_YMin) | int( y < last_y_bin ? e_None : e_YMax) );
					for ( int x = 0, nbx = NBinsX_; x < nbx; ++x, ++pos ) {
						// determine extents in each direction
						BorderFlags_[ pos + 1 ] = zy_flag | int( x ? e_None : e_XMin) | int( x < last_x_bin ? e_None : e_XMax);
					}
				}
			}
		}
	}
}

/// @brief Compute distance from this point to the nearest edge of the bounding box
///        of all given points, provided the point is outside the box.
/// @param point: point of interest
template< typename T >
Real
VoxelGrid<T>::GetSqDistanceOutsideBoundingBox(Vector const & point) const
{
	Vector const & mn( Box_.lower() );
	Vector const & mx( Box_.upper() );
	Real oob_dist_sq( 0.0);
	if ( point.x() < mn.x() ) {
		oob_dist_sq += ( (point.x() - mn.x()) * (point.x() - mn.x()) );
	} else if ( point.x() > mx.x() ) {
		oob_dist_sq += ( (point.x() - mx.x()) * (point.x() - mx.x()) );
	}
	if ( point.y() < mn.y() ) {
		oob_dist_sq += ( (point.y() - mn.y()) * (point.y() - mx.y()) );
	} else if ( point.y() > mx.y() ) {
		oob_dist_sq += ( (point.y() - mx.y()) * (point.y() - mx.y()) );
	}
	if ( point.z() < mn.z() ) {
		oob_dist_sq += ( (point.z() - mn.z()) * (point.z() - mn.z()) );
	} else if ( point.z() > mx.z() ) {
		oob_dist_sq += ( (point.z() - mx.z()) * (point.z() - mx.z()) );
	}
	return oob_dist_sq;
}

/// @brief Compute distance from another grid's bounding box to the nearest edge of this
///        bounding box provided that the MIN and MAX is outside of this bounding box.
/// @param min_box: min of the grid of interest
/// @param min_box: max of the grid of interest
template< typename T >
Real
VoxelGrid<T>::GetSqDistanceOutsideBoundingBox
(
	Vector const & min_box,
	Vector const & max_box
) const
{
	Vector const & mn( Box_.lower());
	Vector const & mx( Box_.upper());
	Real oob_dist_sq( 0.0);
	if ( max_box.x() < mn.x() ) {
		oob_dist_sq += ( (max_box.x() - mn.x()) * (max_box.x() - mn.x()) );
	} else if ( min_box.x() > mx.x() ) {
		oob_dist_sq += ( (min_box.x() - mx.x()) * (min_box.x() - mx.x()) );
	}
	if ( max_box.y() < mn.y() ) {
		oob_dist_sq += ( (max_box.y() - mn.y()) * (max_box.y() - mn.y()) );
	} else if ( min_box.y() > mx.y() ) {
		oob_dist_sq += ( (min_box.y() - mx.y()) * (min_box.y() - mx.y()) );
	}
	if ( max_box.z() < mn.z() ) {
		oob_dist_sq += ( (max_box.z() - mn.z()) * (max_box.z() - mn.z()) );
	} else if ( min_box.z() > mx.z() ) {
		oob_dist_sq += ( (min_box.z() - mx.z()) * (min_box.z() - mx.z()) );
	}
	return oob_dist_sq;
}

/// @brief Find all neighbor pairs of this grid within a specified distance
/// @param neighborhood: range in which to search for neighbors; must be <= resolution
/// @param check_relevance: check and only return neighboring objects which are relevant
template< typename T >
utility::vector1< std::pair< std::pair< T const *, T const * >, Real > >
VoxelGrid<T>::GetNeighborsMultiDimensional
(
	Real const & neighborhood,
	bool check_relevance
) const
{
	runtime_assert_msg(neighborhood <= Resolution_, "Error: GetNeighbors not supported for voxel grids with resolution smaller than neighborhood distance.");
	runtime_assert_msg(neighborhood > 0.0, "Error: Neighborhood search radius must be greater 0.0");

	utility::vector1< std::pair< t_DataTypeCOPtrPair, Real > > neighbors_ret;
	const Real cutoff_sq( neighborhood * neighborhood );
	if ( CacheEdges_ ) {
		auto itr_edges( Edges_.begin() );

		for ( auto itr_vox_obj( Assignments_.begin() ), itr_vox_end( Assignments_.end() );
				itr_vox_obj != itr_vox_end; ++itr_vox_obj, ++itr_edges ) {

			if ( itr_vox_obj->empty() ) {
				continue;
			}
			auto itr_edges_begin( itr_edges->begin() ), itr_edges_end( itr_edges->end() );

			// iterate over all neighbors within the voxel.
			// For high-density voxels (~32+ items) in 3D, it would be more efficient to do an
			// additional labeling of all items as to their sub-octant. Items in the voxel (on-the-fly) or as an auxiliary grid).
			// Any two objects with the same octant number are automatically in the same neighborhood.
			for ( auto itr_vox_obj_a( itr_vox_obj->begin() ), itr_vox_obj_a_end( itr_vox_obj->end() );
					itr_vox_obj_a != itr_vox_obj_a_end; ++itr_vox_obj_a ) {
				if ( !check_relevance || (check_relevance && IsRelevantItem( *(itr_vox_obj_a->first) )) ) {
					Vector const & coord_a( *(itr_vox_obj_a->second) );

					for ( auto itr_vox_obj_b( itr_vox_obj->begin() ); itr_vox_obj_b != itr_vox_obj_a; ++itr_vox_obj_b ) {
						if ( !check_relevance || (check_relevance && IsRelevantItem( *(itr_vox_obj_b->first) )) ) {
							const Real sq_distance( coord_a.distance_squared( *(itr_vox_obj_b->second) ) );
							if ( sq_distance < cutoff_sq ) {
								neighbors_ret.push_back(std::make_pair(std::make_pair(itr_vox_obj_a->first, itr_vox_obj_b->first), std::sqrt(sq_distance)));
							}
						}
					}
					for ( auto itr_lists_b( itr_edges_begin ); itr_lists_b != itr_edges_end; ++itr_lists_b ) {

						for ( auto itr_vox_obj_b( itr_lists_b->begin() ), itr_vox_obj_b_end( itr_lists_b->end() );
								itr_vox_obj_b != itr_vox_obj_b_end; ++itr_vox_obj_b ) {
							if ( !check_relevance || (check_relevance && IsRelevantItem( *(itr_vox_obj_b->first) )) ) {
								const Real sq_distance( coord_a.distance_squared( *(itr_vox_obj_b->second) ) );
								if ( sq_distance < cutoff_sq ) {
									neighbors_ret.push_back(std::make_pair(std::make_pair(itr_vox_obj_a->first, itr_vox_obj_b->first), std::sqrt(sq_distance)));
								}
							}
						}
					}
				}
			}
		}
	} else {
		Size pos( 1 );
		const int z_stride( NBinsX_ * NBinsY_ ), nbx( NBinsX_);
		auto itr_flag( BorderFlags_.begin() );
		for ( auto itr_vox_obj( Assignments_.begin() ), itr_vox_end( Assignments_.end() );
				itr_vox_obj != itr_vox_end; ++itr_vox_obj, ++pos, ++itr_flag ) {
			if ( itr_vox_obj->empty() ) {
				continue;
			}
			auto itr_vox_obj_a_end( itr_vox_obj->end() );

			// iterate over all neighbors within the voxel.
			// For high-density voxels (~32+ items) in 3D, it would be more efficient to do an
			// additional labeling of all items as to their sub-octant. Items in the voxel (on-the-fly) or as an auxiliary grid).
			// Any two objects with the same octant number are automatically in the same neighborhood.
			for ( auto itr_vox_obj_a( itr_vox_obj->begin() ); itr_vox_obj_a != itr_vox_obj_a_end; ++itr_vox_obj_a ) {
				if ( !check_relevance || (check_relevance && IsRelevantItem( *(itr_vox_obj_a->first) )) ) {
					Vector const & coord_a( *(itr_vox_obj_a->second) );

					for ( auto itr_vox_obj_b( itr_vox_obj->begin()); itr_vox_obj_b != itr_vox_obj_a; ++itr_vox_obj_b ) {
						if ( !check_relevance || (check_relevance && IsRelevantItem( *(itr_vox_obj_b->first) )) ) {
							const Real sq_distance( coord_a.distance_squared( *(itr_vox_obj_b->second) ) );
							if ( sq_distance < cutoff_sq ) {
								neighbors_ret.push_back(std::make_pair(std::make_pair(itr_vox_obj_a->first, itr_vox_obj_b->first), std::sqrt(sq_distance)));
							}
						}
					}
				}
			}

			const int flag( *itr_flag );
			const int max_z( flag & e_ZMax ? 1 : 2);
			const int min_y( flag & e_YMin ? 0 : -1), max_y( flag & e_YMax ? 1 : 2);
			const int min_x( flag & e_XMin ? 0 : -1), max_x( flag & e_XMax ? 1 : 2);

			for ( int z = 0, zy_off_init = min_y * nbx, zy_off_end = max_y * nbx;
					z < max_z; ++z, zy_off_init += z_stride, zy_off_end += z_stride ) {
				if ( zy_off_end < 0 ) {
					continue;
				}

				for ( int zy_off = zy_off_init; zy_off < zy_off_end; zy_off += nbx ) {
					if ( zy_off < 0 ) {
						continue;
					}

					for ( int zyx_off = zy_off + min_x, zyx_off_mx = zy_off + max_x; zyx_off < zyx_off_mx; ++zyx_off ) {
						if ( zyx_off > 0 ) {
							utility::vector1< t_DataTypeVectorCOPtrPair > const & adjacent_vox( Assignments_[ pos + zyx_off ] );
							if ( adjacent_vox.empty() ) {
								continue;
							}

							auto itr_vox_obj_b_beg( adjacent_vox.begin() ), itr_vox_obj_b_end( adjacent_vox.end() );
							for ( auto itr_vox_obj_a( itr_vox_obj->begin() ); itr_vox_obj_a != itr_vox_obj_a_end; ++itr_vox_obj_a ) {
								if ( !check_relevance || (check_relevance && IsRelevantItem( *(itr_vox_obj_a->first) )) ) {
									Vector const &coord_a( *(itr_vox_obj_a->second) );

									for ( auto itr_vox_obj_b( itr_vox_obj_b_beg ); itr_vox_obj_b != itr_vox_obj_b_end; ++itr_vox_obj_b ) {
										if ( !check_relevance || (check_relevance && IsRelevantItem( *(itr_vox_obj_b->first) )) ) {
											const Real sq_distance( coord_a.distance_squared( *(itr_vox_obj_b->second) ) );
											if ( sq_distance < cutoff_sq ) {
												neighbors_ret.push_back(std::make_pair(std::make_pair(itr_vox_obj_a->first, itr_vox_obj_b->first), std::sqrt(sq_distance)));
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return neighbors_ret;
}

/// @brief Find all neighbors of this object within a specified distance
/// @param obj: object to look for neighbors
/// @param neighborhood: range in which to search for neighbors; must be <= resolution
/// @param check_relevance: check and only return neighboring objects which are relevant
template< typename T >
utility::vector1< std::pair< T const *, Real > >
VoxelGrid<T>::GetNeighborsMultiDimensional
(
	t_DataType const & obj,
	Real const & neighborhood,
	bool check_relevance
) const
{
	runtime_assert_msg(neighborhood <= Resolution_, "Error: GetNeighbors not supported for voxel grids with resolution smaller than neighborhood distance.");
	runtime_assert_msg(neighborhood > 0.0, "Error: Neighborhood search radius must be greater 0.0");

	utility::vector1< std::pair< t_DataType const *, Real > > neighbors_ret;
	const Real cutoff_sq( neighborhood * neighborhood );
	Vector const * coord_ptr( ExtractPosition( obj ) );
	if ( !coord_ptr ) {
		return neighbors_ret;
	}
	Vector const & coord(*coord_ptr);
	if ( GetSqDistanceOutsideBoundingBox( coord ) >= cutoff_sq ) {
		return neighbors_ret;
	}

	const Size pos( GetIndex( coord ) );
	auto itr_vox_obj( Assignments_.begin() + pos - 1);

	// iterate over all neighbors within the voxel
	for ( auto itr_vox_obj_a( itr_vox_obj->begin()), itr_vox_obj_vec_end( itr_vox_obj->end());
			itr_vox_obj_a != itr_vox_obj_vec_end; ++itr_vox_obj_a ) {
		if ( !check_relevance || (check_relevance && IsRelevantItem( *(itr_vox_obj_a->first) )) ) {
			const Real sq_distance( coord.distance_squared( *(itr_vox_obj_a->second) ) );
			if ( sq_distance < cutoff_sq && !IsSameItem( *(itr_vox_obj_a->first), obj) ) {
				neighbors_ret.push_back( std::make_pair(itr_vox_obj_a->first, std::sqrt(sq_distance)) );
			}
		}
	}

	if ( CacheEdges_ ) {
		for ( auto itr_lists_b( EdgesComplete_[ pos ].begin() ), itr_lists_b_end( EdgesComplete_[ pos ].end());
				itr_lists_b != itr_lists_b_end; ++itr_lists_b ) {

			for ( auto itr_vox_obj_b( itr_lists_b->begin() ), itr_vox_obj_b_end( itr_lists_b->end() );
					itr_vox_obj_b != itr_vox_obj_b_end; ++itr_vox_obj_b ) {
				if ( !check_relevance || (check_relevance && IsRelevantItem( *(itr_vox_obj_b->first) )) ) {
					const Real sq_distance( coord.distance_squared( *(itr_vox_obj_b->second) ) );
					if ( sq_distance < cutoff_sq ) {
						neighbors_ret.push_back( std::make_pair(itr_vox_obj_b->first, std::sqrt(sq_distance)) );
					}
				}
			}
		}
	} else {
		// determine extents in each direction
		const int flag( BorderFlags_[ pos ] );
		const int min_z( flag & e_ZMin ? 0 : -1), max_z( flag & e_ZMax ? 1 : 2);
		const int min_y( flag & e_YMin ? 0 : -1), max_y( flag & e_YMax ? 1 : 2);
		const int min_x( flag & e_XMin ? 0 : -1), max_x( flag & e_XMax ? 1 : 2);
		const int z_stride( NBinsX_ * NBinsY_ );

		for ( int z = min_z, zy_off_init = min_z * z_stride + min_y * NBinsX_, zy_off_end = min_z * z_stride + max_y * NBinsX_;
				z < max_z; ++z, zy_off_init += z_stride, zy_off_end += z_stride ) {

			for ( int zy_off = zy_off_init; zy_off < zy_off_end; zy_off += NBinsX_ ) {

				for ( int zyx_off = zy_off + min_x, zyx_off_mx = zy_off + max_x; zyx_off < zyx_off_mx; ++zyx_off ) {
					if ( zyx_off ) {
						utility::vector1< t_DataTypeVectorCOPtrPair > const & adjacent_vox( Assignments_[ pos + zyx_off ] );

						for ( auto itr_vox_obj_b( adjacent_vox.begin() ), itr_vox_obj_b_end( adjacent_vox.end() );
								itr_vox_obj_b != itr_vox_obj_b_end; ++itr_vox_obj_b ) {
							if ( !check_relevance || (check_relevance && IsRelevantItem( *(itr_vox_obj_b->first) )) ) {
								const Real sq_distance( coord.distance_squared( *(itr_vox_obj_b->second) ) );
								if ( sq_distance < cutoff_sq ) {
									neighbors_ret.push_back( std::make_pair(itr_vox_obj_b->first, std::sqrt(sq_distance)) );
								}
							}
						}
					}
				}
			}
		}
	}
	return neighbors_ret;
}

/// @brief Find all neighbors between the input and this grid within a specified distance
/// @param grid: input grid in which to look for objects that are neighbors of this grid
/// @param neighborhood: range in which to search for neighbors; must be <= resolution
/// @param check_relevance: check and only return neighboring objects which are relevant
template< typename T >
utility::vector1< std::pair< std::pair< T const *, T const * >, Real > >
VoxelGrid<T>::GetNeighborsMultiDimensional
(
	VoxelGrid< t_DataType > const & grid,
	Real const & neighborhood,
	bool check_relevance
) const
{
	runtime_assert_msg(neighborhood <= Resolution_, "Error: GetNeighbors not supported for voxel grids with resolution smaller than neighborhood distance.");
	runtime_assert_msg(neighborhood > 0.0, "Error: Neighborhood search radius must be greater 0.0");

	utility::vector1< std::pair< t_DataTypeCOPtrPair, Real > > neighbors_ret;
	const Real cutoff_sq( neighborhood * neighborhood );

	// test whether it is even possible that any points on the grids are within the cutoff distance of one another
	if ( GetSqDistanceOutsideBoundingBox( grid.GetBoundingBox().lower(), grid.GetBoundingBox().upper()) > cutoff_sq ) {
		return neighbors_ret;
	}

	for ( auto itr_other_grid( grid.GetGridItems().begin() ), itr_other_grid_end( grid.GetGridItems().end() );
			itr_other_grid != itr_other_grid_end; ++itr_other_grid ) {

		for ( auto itr_obj_a( itr_other_grid->begin() ), itr_obj_a_end( itr_other_grid->end() );
				itr_obj_a != itr_obj_a_end; ++itr_obj_a ) {
			if ( !check_relevance || (check_relevance && IsRelevantItem( *(itr_obj_a->first) )) ) {

				Vector const & coord( *(itr_obj_a->second) );
				if ( GetSqDistanceOutsideBoundingBox(coord) >= cutoff_sq ) {
					continue;
				}

				t_DataType const * obj_a( itr_obj_a->first );
				const Size pos( GetIndex( coord ) );
				auto itr_vox_obj( Assignments_.begin() + pos - 1 );

				// iterate over all neighbors within the voxel
				for ( auto itr_vox_obj_a( itr_vox_obj->begin() ), itr_vox_obj_vec_end( itr_vox_obj->end() );
						itr_vox_obj_a != itr_vox_obj_vec_end; ++itr_vox_obj_a ) {
					if ( !check_relevance || (check_relevance && IsRelevantItem( *(itr_vox_obj_a->first) )) ) {
						const Real sq_distance( coord.distance_squared( *(itr_vox_obj_a->second) ) );
						if ( sq_distance < cutoff_sq && !IsSameItem( *(itr_vox_obj_a->first), *obj_a ) ) {
							neighbors_ret.push_back( std::make_pair(std::make_pair(obj_a, itr_vox_obj_a->first), std::sqrt(sq_distance)) );
						}
					}
				}

				if ( CacheEdges_ ) {
					for ( auto itr_lists_b( EdgesComplete_[ pos ].begin() ), itr_lists_b_end( EdgesComplete_[ pos ].end() );
							itr_lists_b != itr_lists_b_end; ++itr_lists_b ) {

						for ( auto itr_vox_obj_b( itr_lists_b->begin() ), itr_vox_obj_b_end( itr_lists_b->end() );
								itr_vox_obj_b != itr_vox_obj_b_end; ++itr_vox_obj_b ) {
							if ( !check_relevance || (check_relevance && IsRelevantItem( *(itr_vox_obj_b->first) )) ) {
								const Real sq_distance( coord.distance_squared( *(itr_vox_obj_b->second) ) );
								if ( sq_distance < cutoff_sq ) {
									neighbors_ret.push_back( std::make_pair(std::make_pair(obj_a, itr_vox_obj_b->first), std::sqrt(sq_distance)) );
								}
							}
						}
					}
				} else {
					// determine extents in each direction
					const int flag( BorderFlags_[ pos ] );
					const int min_z( flag & e_ZMin ? 0 : -1), max_z( flag & e_ZMax ? 1 : 2);
					const int min_y( flag & e_YMin ? 0 : -1), max_y( flag & e_YMax ? 1 : 2);
					const int min_x( flag & e_XMin ? 0 : -1), max_x( flag & e_XMax ? 1 : 2);
					const int z_stride( NBinsX_ * NBinsY_);

					for ( int z = min_z, zy_off_init = min_z * z_stride + min_y * NBinsX_, zy_off_end = min_z * z_stride + max_y * NBinsX_;
							z < max_z; ++z, zy_off_init += z_stride, zy_off_end += z_stride ) {

						for ( int zy_off = zy_off_init; zy_off < zy_off_end; zy_off += NBinsX_ ) {

							for ( int zyx_off = zy_off + min_x, zyx_off_mx = zy_off + max_x; zyx_off < zyx_off_mx; ++zyx_off ) {

								if ( zyx_off ) {
									utility::vector1< t_DataTypeVectorCOPtrPair >const & adjacent_vox( Assignments_[ pos + zyx_off ] );

									for ( auto itr_vox_obj_b( adjacent_vox.begin() ), itr_vox_obj_b_end( adjacent_vox.end() );
											itr_vox_obj_b != itr_vox_obj_b_end; ++itr_vox_obj_b ) {
										if ( !check_relevance || (check_relevance && IsRelevantItem( *(itr_vox_obj_b->first) )) ) {
											const Real sq_distance( coord.distance_squared( *(itr_vox_obj_b->second) ) );
											if ( sq_distance < cutoff_sq ) {
												neighbors_ret.push_back( std::make_pair(std::make_pair(obj_a, itr_vox_obj_b->first), std::sqrt(sq_distance)) );
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return neighbors_ret;
}

/// @brief Has given object a neighbor in the grid within a specified distance
/// @param obj: object to look for neighbors
/// @param neighborhood: range in which to search for neighbors; must be <= resolution
/// @param check_relevance: check and only return true if neighboring objects are relevant
template< typename T >
bool
VoxelGrid<T>::HasNeighborsMultiDimensional
(
	t_DataType const & obj,
	Real const & neighborhood,
	bool check_relevance
) const
{
	runtime_assert_msg(neighborhood <= Resolution_, "Error: HasNeighbors not supported for voxel grids with resolution smaller than neighborhood distance.");
	runtime_assert_msg(neighborhood > 0.0, "Error: Neighborhood search radius must be greater 0.0");

	const Real cutoff_sq( neighborhood * neighborhood );
	Vector const * coord_ptr( ExtractPosition( obj ) );
	if ( !coord_ptr ) {
		return false;
	}
	Vector const & coord(*coord_ptr);
	if ( GetSqDistanceOutsideBoundingBox( coord ) >= cutoff_sq ) {
		return false;
	}

	const Size pos( GetIndex( coord ) );
	auto itr_vox_obj( Assignments_.begin() + pos - 1);

	// iterate over all neighbors within the voxel
	for ( auto itr_vox_obj_a( itr_vox_obj->begin()), itr_vox_obj_vec_end( itr_vox_obj->end());
			itr_vox_obj_a != itr_vox_obj_vec_end; ++itr_vox_obj_a ) {
		if ( !check_relevance || (check_relevance && IsRelevantItem( *(itr_vox_obj_a->first) )) ) {
			const Real sq_distance( coord.distance_squared( *(itr_vox_obj_a->second) ) );
			if ( sq_distance < cutoff_sq && !IsSameItem( *(itr_vox_obj_a->first), obj) ) {
				return true;
			}
		}
	}

	if ( CacheEdges_ ) {
		for ( auto itr_lists_b( EdgesComplete_[ pos ].begin() ), itr_lists_b_end( EdgesComplete_[ pos ].end());
				itr_lists_b != itr_lists_b_end; ++itr_lists_b ) {

			for ( auto itr_vox_obj_b( itr_lists_b->begin() ), itr_vox_obj_b_end( itr_lists_b->end() );
					itr_vox_obj_b != itr_vox_obj_b_end; ++itr_vox_obj_b ) {
				if ( !check_relevance || (check_relevance && IsRelevantItem( *(itr_vox_obj_b->first) )) ) {
					const Real sq_distance( coord.distance_squared( *(itr_vox_obj_b->second) ) );
					if ( sq_distance < cutoff_sq ) {
						return true;
					}
				}
			}
		}
	} else {
		// determine extents in each direction
		const int flag( BorderFlags_[ pos ] );
		const int min_z( flag & e_ZMin ? 0 : -1), max_z( flag & e_ZMax ? 1 : 2);
		const int min_y( flag & e_YMin ? 0 : -1), max_y( flag & e_YMax ? 1 : 2);
		const int min_x( flag & e_XMin ? 0 : -1), max_x( flag & e_XMax ? 1 : 2);
		const int z_stride( NBinsX_ * NBinsY_ );

		for ( int z = min_z, zy_off_init = min_z * z_stride + min_y * NBinsX_, zy_off_end = min_z * z_stride + max_y * NBinsX_;
				z < max_z; ++z, zy_off_init += z_stride, zy_off_end += z_stride ) {

			for ( int zy_off = zy_off_init; zy_off < zy_off_end; zy_off += NBinsX_ ) {

				for ( int zyx_off = zy_off + min_x, zyx_off_mx = zy_off + max_x; zyx_off < zyx_off_mx; ++zyx_off ) {
					if ( zyx_off ) {
						utility::vector1< t_DataTypeVectorCOPtrPair > const & adjacent_vox( Assignments_[ pos + zyx_off ] );

						for ( auto itr_vox_obj_b( adjacent_vox.begin() ), itr_vox_obj_b_end( adjacent_vox.end() );
								itr_vox_obj_b != itr_vox_obj_b_end; ++itr_vox_obj_b ) {
							if ( !check_relevance || (check_relevance && IsRelevantItem( *(itr_vox_obj_b->first) )) ) {
								const Real sq_distance( coord.distance_squared( *(itr_vox_obj_b->second) ) );
								if ( sq_distance < cutoff_sq ) {
									return true;
								}
							}
						}
					}
				}
			}
		}
	}
	return false;
}

/// @brief Find all neighbor pairs of this 1D grid within a specified distance
/// @param neighborhood: range in which to search for neighbors; must be <= resolution
/// @param check_relevance: check and only return neighboring objects which are relevant
template< typename T >
utility::vector1< std::pair< std::pair< T const *, T const * >, Real > >
VoxelGrid<T>::GetNeighbors1D
(
	Real const & neighborhood,
	bool check_relevance
) const
{
	runtime_assert_msg(neighborhood <= Resolution_, "Error: GetNeighbors not supported for voxel grids with resolution smaller than neighborhood distance.");
	runtime_assert_msg(neighborhood > 0.0, "Error: Neighborhood search radius must be greater 0.0");

	utility::vector1< std::pair< t_DataTypeCOPtrPair, Real > > neighbors_ret;

	const Real cutoff_sq( neighborhood * neighborhood );
	Size index_z_end( std::ceil(neighborhood / MinResolution_) + 1);
	for ( auto itr_zvec( Assignments_.begin() ), itr_zvec_end( Assignments_.end() );
			itr_zvec != itr_zvec_end; ++itr_zvec, ++index_z_end ) {

		auto itr_zvec_b_end( index_z_end <= Assignments_.size() ? Assignments_.begin() + index_z_end : itr_zvec_end );
		for ( auto itr_zobj( itr_zvec->begin() ), itr_zobj_end( itr_zvec->end() );
				itr_zobj != itr_zobj_end; ++itr_zobj ) {

			if ( !check_relevance || (check_relevance && IsRelevantItem( *(itr_zobj->first) )) ) {
				// for each object
				t_DataType const * obj_ptr( itr_zobj->first );
				Vector const & coord_ref( *(itr_zobj->second) );
				auto itr_zvec_b( itr_zvec );
				for ( auto itr_short( itr_zvec_b->begin() ); itr_short != itr_zobj; ++itr_short ) {
					if ( !check_relevance || (check_relevance && IsRelevantItem( *(itr_short->first) )) ) {
						const Real dist_sq( coord_ref.distance_squared( *(itr_short->second) ) );
						if ( dist_sq < cutoff_sq ) {
							neighbors_ret.push_back( std::make_pair(std::make_pair(obj_ptr, itr_short->first), std::sqrt(dist_sq)) );
						}
					}
				}
				for ( ++itr_zvec_b; itr_zvec_b != itr_zvec_b_end; ++itr_zvec_b ) {

					for ( auto itr_zobj_b( itr_zvec_b->begin() ), itr_zobj_b_end( itr_zvec_b->end() );
							itr_zobj_b != itr_zobj_b_end; ++itr_zobj_b ) {
						if ( !check_relevance || (check_relevance && IsRelevantItem( *(itr_zobj_b->first) )) ) {
							const Real dist_sq( coord_ref.distance_squared( *(itr_zobj_b->second) ) );
							if ( dist_sq < cutoff_sq ) {
								neighbors_ret.push_back( std::make_pair(std::make_pair(obj_ptr, itr_zobj_b->first), std::sqrt(dist_sq)) );
							}
						}
					}
				}
			}
		}
	}
	return neighbors_ret;
}

/// @brief Find all neighbors of this object within a specified distance
/// @param obj: object to look for neighbors
/// @param neighborhood: range in which to search for neighbors; must be <= resolution
/// @param check_relevance: check and only return neighboring objects which are relevant
template< typename T >
utility::vector1< std::pair< T const *, Real > >
VoxelGrid<T>::GetNeighbors1D
(
	t_DataType const & obj,
	Real const & neighborhood,
	bool check_relevance
) const
{
	runtime_assert_msg(neighborhood <= Resolution_, "Error: GetNeighbors not supported for voxel grids with resolution smaller than neighborhood distance.");
	runtime_assert_msg(neighborhood > 0.0, "Error: Neighborhood search radius must be greater 0.0");

	utility::vector1< std::pair< t_DataType const *, Real > > neighbors_ret;
	Vector const * coord_ptr( ExtractPosition(obj) );
	if ( !coord_ptr ) {
		return neighbors_ret;
	}
	Vector const & coord(*coord_ptr);
	const Real cutoff_sq( neighborhood * neighborhood );
	if ( GetSqDistanceOutsideBoundingBox( coord ) >= cutoff_sq ) {
		return neighbors_ret;
	}

	const Size n_bins_1d( std::ceil(neighborhood / MinResolution_) );
	const Size pos( GetIndex(coord) - 1 );
	auto itr_zvec_b( Assignments_.begin() + pos ), itr_zvec_b_base( Assignments_.begin() + pos );
	// We can't use  std::min( Assignments_.begin() + pos + n_bins_1d + 1, Assignments_.end() )
	auto itr_zvec_b_end( pos + n_bins_1d >= Assignments_.size() ? Assignments_.end() : Assignments_.begin() + pos + n_bins_1d + 1 );

	for ( auto itr_zvec_a( itr_zvec_b->begin() ), itr_zvec_a_end( itr_zvec_b->end() );
			itr_zvec_a != itr_zvec_a_end; ++itr_zvec_a ) {
		if ( !check_relevance || (check_relevance && IsRelevantItem( *(itr_zvec_a->first) )) ) {
			const Real dist_sq( coord.distance_squared( *(itr_zvec_a->second) ) );
			if ( dist_sq < cutoff_sq && !IsSameItem( *(itr_zvec_a->first), obj) ) {
				neighbors_ret.push_back( std::make_pair(itr_zvec_a->first, std::sqrt(dist_sq)) );
			}
		}
	}

	// We can't use this former construction here because if Assignments_.beign() - n_bins_1d + pos
	// is less than Assignments_.begin(), it is nonsense.
	//  std::max( Assignments_.begin() - n_bins_1d + pos, Assignments_.begin() )
	for ( itr_zvec_b = ( ( pos < n_bins_1d ) ? Assignments_.begin() : Assignments_.begin() + ( pos - n_bins_1d ) );
			itr_zvec_b != itr_zvec_b_end; ++itr_zvec_b ) {
		if ( itr_zvec_b == itr_zvec_b_base ) {
			continue;
		}
		for ( auto itr_zobj_b( itr_zvec_b->begin() ), itr_zobj_b_end( itr_zvec_b->end() );
				itr_zobj_b != itr_zobj_b_end; ++itr_zobj_b ) {
			if ( !check_relevance || (check_relevance && IsRelevantItem( *(itr_zobj_b->first) )) ) {
				const Real dist_sq( coord.distance_squared( *(itr_zobj_b->second) ) );
				if ( dist_sq < cutoff_sq ) {
					neighbors_ret.push_back( std::make_pair(itr_zobj_b->first, std::sqrt(dist_sq)) );
				}
			}
		}
	}
	return neighbors_ret;
}


/// @brief Find all neighbors between the input and this grid within a specified distance
/// @param grid: input grid in which to look for objects that are neighbors of this grid
/// @param neighborhood: range in which to search for neighbors; must be <= resolution
/// @param check_relevance: check and only return neighboring objects which are relevant
template< typename T >
utility::vector1< std::pair< std::pair< T const *, T const * >, Real > >
VoxelGrid<T>::GetNeighbors1D
(
	VoxelGrid< t_DataType> const & grid,
	Real const & neighborhood,
	bool check_relevance
) const
{
	runtime_assert_msg(neighborhood <= Resolution_, "Error: GetNeighbors not supported for voxel grids with resolution smaller than neighborhood distance.");
	runtime_assert_msg(neighborhood > 0.0, "Error: Neighborhood search radius must be greater 0.0");

	utility::vector1< std::pair< t_DataTypeCOPtrPair, Real > > neighbors_ret;
	const Real cutoff_sq( neighborhood * neighborhood );
	// test whether it is even possible that any points on the grids are within the cutoff distance of one another
	if ( GetSqDistanceOutsideBoundingBox(grid.GetBoundingBox().lower(), grid.GetBoundingBox().upper()) > cutoff_sq ) {
		return neighbors_ret;
	}

	const Size n_bins_1d( std::ceil( neighborhood / MinResolution_ ) );
	for ( auto itr_zvec( grid.GetGridItems().begin() ), itr_zvec_end( grid.GetGridItems().end() );
			itr_zvec != itr_zvec_end; ++itr_zvec ) {

		for ( auto itr_obj_a( itr_zvec->begin() ), itr_obj_a_end( itr_zvec->end() );
				itr_obj_a != itr_obj_a_end; ++itr_obj_a ) {
			if ( !check_relevance || (check_relevance && IsRelevantItem( *(itr_obj_a->first) )) ) {
				Vector const & coord( *(itr_obj_a->second) );
				if ( GetSqDistanceOutsideBoundingBox( coord ) >= cutoff_sq ) {
					continue;
				}

				t_DataType const * obj_a( itr_obj_a->first );
				const Size pos( GetIndex(coord) - 1 );
				auto itr_zvec_b( Assignments_.begin() + pos ), itr_zvec_b_base( Assignments_.begin() + pos );
				// can't do std::min( Assignments_.begin() + pos + n_bins_1d + 1, Assignments_.end() )
				auto itr_zvec_b_end( pos + n_bins_1d >= Assignments_.size() ? Assignments_.end() : Assignments_.begin() + pos + n_bins_1d + 1 );

				for ( auto itr_zvec_a( itr_zvec_b->begin() ), itr_zvec_a_end( itr_zvec_b->end() );
						itr_zvec_a != itr_zvec_a_end; ++itr_zvec_a ) {
					if ( !check_relevance || (check_relevance && IsRelevantItem( *(itr_zvec_a->first) )) ) {
						const Real dist_sq( coord.distance_squared( *(itr_zvec_a->second) ) );
						if ( dist_sq < cutoff_sq && !IsSameItem( *(itr_zvec_a->first), *obj_a) ) {
							neighbors_ret.push_back( std::make_pair(std::make_pair(itr_zvec_a->first, obj_a), std::sqrt(dist_sq)) );
						}
					}
				}

				// can't do std::max( Assignments_.begin() - n_bins_1d + pos, Assignments_.begin()
				for ( itr_zvec_b = ( ( pos >= n_bins_1d ) ? Assignments_.begin() + ( pos - n_bins_1d ) : Assignments_.begin() );
						itr_zvec_b != itr_zvec_b_end; ++itr_zvec_b ) {
					if ( itr_zvec_b == itr_zvec_b_base ) {
						continue;
					}

					for ( auto itr_zobj_b( itr_zvec_b->begin() ), itr_zobj_b_end( itr_zvec_b->end() );
							itr_zobj_b != itr_zobj_b_end; ++itr_zobj_b ) {
						if ( !check_relevance || (check_relevance && IsRelevantItem( *(itr_zobj_b->first) )) ) {
							const Real dist_sq( coord.distance_squared( *(itr_zobj_b->second) ) );
							if ( dist_sq < cutoff_sq ) {
								neighbors_ret.push_back( std::make_pair(std::make_pair(itr_zobj_b->first, obj_a), std::sqrt(dist_sq)) );
							}
						}
					}
				}
			}
		}
	}
	return neighbors_ret;
}

/// @brief Has given object a neighbor in the grid within a specified distance
/// @param obj: object to look for neighbors
/// @param neighborhood: range in which to search for neighbors; must be <= resolution
/// @param check_relevance: check and only return true if neighboring objects are relevant
template< typename T >
bool
VoxelGrid<T>::HasNeighbors1D
(
	t_DataType const & obj,
	Real const & neighborhood,
	bool check_relevance
) const
{
	runtime_assert_msg(neighborhood <= Resolution_, "Error: HasNeighbors not supported for voxel grids with resolution smaller than neighborhood distance.");
	runtime_assert_msg(neighborhood > 0.0, "Error: Neighborhood search radius must be greater 0.0");

	Vector const * coord_ptr( ExtractPosition(obj) );
	if ( !coord_ptr ) {
		return false;
	}

	Vector const & coord(*coord_ptr);
	const Real cutoff_sq( neighborhood * neighborhood );
	if ( GetSqDistanceOutsideBoundingBox( coord ) >= cutoff_sq ) {
		return false;
	}

	const Size n_bins_1d( std::ceil(neighborhood / MinResolution_) );
	const Size pos( GetIndex(coord) - 1 );
	auto itr_zvec_b( Assignments_.begin() + pos ), itr_zvec_b_base( Assignments_.begin() + pos );
	auto itr_zvec_b_end( pos + n_bins_1d >= Assignments_.size() ? Assignments_.end() : Assignments_.begin() + pos + n_bins_1d + 1 );

	// Loop through the object's own voxel
	for ( auto itr_zvec_a( itr_zvec_b->begin() ), itr_zvec_a_end( itr_zvec_b->end() );
			itr_zvec_a != itr_zvec_a_end; ++itr_zvec_a ) {
		if ( !check_relevance || (check_relevance && IsRelevantItem( *(itr_zvec_a->first) )) ) {
			const Real dist_sq( coord.distance_squared( *(itr_zvec_a->second) ) );
			if ( dist_sq < cutoff_sq && !IsSameItem( *(itr_zvec_a->first), obj) ) {
				return true;
			}
		}
	}

	// Look through neighboring voxels
	for ( itr_zvec_b = ( ( pos >= n_bins_1d ) ? Assignments_.begin() - n_bins_1d + pos : Assignments_.begin() );
			itr_zvec_b != itr_zvec_b_end; ++itr_zvec_b ) {
		if ( itr_zvec_b == itr_zvec_b_base ) {
			continue;
		}
		for ( auto itr_zobj_b( itr_zvec_b->begin() ), itr_zobj_b_end( itr_zvec_b->end() );
				itr_zobj_b != itr_zobj_b_end; ++itr_zobj_b ) {
			if ( !check_relevance || (check_relevance && IsRelevantItem( *(itr_zobj_b->first) )) ) {
				const Real dist_sq( coord.distance_squared( *(itr_zobj_b->second) ) );
				if ( dist_sq < cutoff_sq ) {
					return true;
				}
			}
		}
	}
	return false;
}

} // namespace numeric

#endif // INCLUDED_numeric_VoxelGrid_impl_hh
