// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/BumpGrid.cc
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini


#include <protocols/match/BumpGrid.hh> // typedefs for Bin3D.

#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace match {

static thread_local basic::Tracer TR( "protocols.BumpGrid" );

/// Creation and initialization
Bool3DGrid::Bool3DGrid() :
	bb_( Vector(0.0) ),
	bb_extended_( Vector(0.0) ),
	bin_width_( 1.0 ),
	bin_width_2x_( 2.0 ),
	dimsizes_( /* 0 */ ),
	dimprods_( /* 0 */ ),
	halfdimsizes_( /* 0 */ ),
	halfdimprods_( /* 0 */ )
{}

Bool3DGrid::~Bool3DGrid() {}

void
Bool3DGrid::set_bounding_box( BoundingBox const & bb ) {
	bb_ = bb;
	reset_grid();
}

void Bool3DGrid::set_bin_width( Real width ) {
	bin_width_ = width;
	reset_grid();
}

Bool3DGrid
Bool3DGrid::create_grid_for_sphere( Vector const & center, Real radius ) const
{
	Bool3DGrid newgrid;
	newgrid.set_bin_width( bin_width_ );

	/// Make the axes ling up with this grid.
	Vector local_center = center - bb_.lower();
	Vector low_corner = local_center - radius;

	assert( low_corner.x() >= 0.0 && low_corner.y() >= 0.0 && low_corner.z() >= 0.0 );

	Bin3D halfgrid_of_lowcorner;
	halfgrid_of_lowcorner[ 1 ] = static_cast< Size > ( low_corner.x() / ( 2 * bin_width_ ) );
	halfgrid_of_lowcorner[ 2 ] = static_cast< Size > ( low_corner.y() / ( 2 * bin_width_ ) );
	halfgrid_of_lowcorner[ 3 ] = static_cast< Size > ( low_corner.z() / ( 2 * bin_width_ ) );

	Vector aligned_low_corner( bb_.lower() );
	aligned_low_corner.x() += bin_width_2x_ * halfgrid_of_lowcorner[ 1 ];
	aligned_low_corner.y() += bin_width_2x_ * halfgrid_of_lowcorner[ 2 ];
	aligned_low_corner.z() += bin_width_2x_ * halfgrid_of_lowcorner[ 3 ];

	Vector upper_corner = local_center + radius;
	newgrid.set_bounding_box( BoundingBox( aligned_low_corner, upper_corner ) );
	return newgrid;
}

/// @brief create a grid for the input bounding box that aligns to this grid
Bool3DGrid
Bool3DGrid::create_grid_for_bb( BoundingBox const & bb )
{
	Bool3DGrid newgrid;
	newgrid.set_bin_width( bin_width_ );

	/// Make the axes line up with this grid.
	Vector local_center = bb.lower() - bb_.lower();

	utility::fixedsizearray1< int, 3 > floored_lower;
	floored_lower[ 1 ] = static_cast< Size > ( local_center.x() / bin_width_ );
	floored_lower[ 2 ] = static_cast< Size > ( local_center.y() / bin_width_ );
	floored_lower[ 3 ] = static_cast< Size > ( local_center.z() / bin_width_ );

	Vector aligned_lower(
		bb_.lower().x() + floored_lower[ 1 ] * bin_width_,
		bb_.lower().y() + floored_lower[ 2 ] * bin_width_,
		bb_.lower().z() + floored_lower[ 3 ] * bin_width_ );

	newgrid.set_bounding_box( BoundingBox( aligned_lower, bb.upper() ) );
	return newgrid;
}


/// @details Sets the value of a voxel to true if all eight corners
/// of the voxel are contained by the sphere.
void
Bool3DGrid::or_by_sphere_conservative( Vector const & center, Real radius )
{
	using namespace core;

	Real const rad2 = radius * radius;
	Real const rad_minus_cendis =  radius - half_bin_width_root_three_;

	/// If the radius is so small that no voxel could be completly covered by a sphere,
	/// then find the handful of voxels that might be overlapped by this tiny sphere,
	/// and process them separately than if rad_minus_cendis were positive.
	if ( rad_minus_cendis < 0 ) {
		/// unimplemented at the moment.
		return;
	}

	Real rad_minus_cendis2 = rad_minus_cendis * rad_minus_cendis;
	Real const rad_plus_cendis = half_bin_width_root_three_ + radius;
	Real const rad_plus_cendis2 = rad_plus_cendis * rad_plus_cendis;

	Bin3D bin;

	Vector relpos = center - bb_extended_.lower();
	SSize xlo( static_cast< SSize > ((relpos.x() - radius) / bin_width_ ));
	SSize xhi( static_cast< SSize > ((relpos.x() + radius) / bin_width_ ) + 1);
	SSize ylo( static_cast< SSize > ((relpos.y() - radius) / bin_width_ ));
	SSize yhi( static_cast< SSize > ((relpos.y() + radius) / bin_width_ ) + 1);
	SSize zlo( static_cast< SSize > ((relpos.z() - radius) / bin_width_ ));
	SSize zhi( static_cast< SSize > ((relpos.z() + radius) / bin_width_ ) + 1);

	Size uxlo = xlo < 0 ? 0 : xlo;
	Size uxhi = xhi < 0 ? 0 : ( xhi > SSize( dimsizes_[ 1 ] ) ? dimsizes_[ 1 ] : xhi ) ;
	Size uylo = ylo < 0 ? 0 : ylo;
	Size uyhi = yhi < 0 ? 0 : ( yhi > SSize( dimsizes_[ 2 ] ) ? dimsizes_[ 2 ] : yhi ) ;
	Size uzlo = zlo < 0 ? 0 : zlo;
	Size uzhi = zhi < 0 ? 0 : ( zhi > SSize( dimsizes_[ 3 ] ) ? dimsizes_[ 3 ] : zhi ) ;

	for ( Size ii = uxlo; ii < uxhi; ++ii ) {
		bin[ 1 ] = ii;
		for ( Size jj = uylo; jj < uyhi; ++jj ) {
			bin[ 2 ] = jj;
			for ( Size kk = uzlo; kk < uzhi; ++kk ) {
				bin[ 3 ] = kk;

				Vector bin_center_point( bin_center( bin ));
				Real sphere_cent_to_voxel_cent2 = bin_center_point.distance_squared( center );
				if (  sphere_cent_to_voxel_cent2 > rad_plus_cendis2 ) continue;

				bool all_contained( true );
				if ( sphere_cent_to_voxel_cent2 > rad_minus_cendis2 ) { // Only look at all the corners if it's possible not all corners are covered by the sphere.
					utility::fixedsizearray1< Vector, 8 > grid_corners = corners( bin );
					for ( Size ll = 1; ll <= 8; ++ll ) {
						Real d2 = center.distance_squared( grid_corners[ ll ] );
						if ( d2 > rad2 ) { all_contained = false; break; }
					}
				}
				if ( all_contained ) {
					set_value_for_bin( bin, true );
				}
			}
		}
	}

}

/// @details Sets the value of a voxel to true if any volume of a voxel is inside
/// the sphere.  Handles edge cases properly where a voxel is glanced by a sphere
/// at a face, but where the corners of the voxel are not contained by the sphere,
/// and when a sphere is entirely contained within a voxel.
void
Bool3DGrid::or_by_sphere_liberal( Vector const & center, Real radius )
{
/*	Real rad2 = radius * radius;
	Bin3D bin;
	for ( Size ii = 0; ii < dimsizes_[ 1 ]; ++ii ) {
		bin[ 1 ] = ii;
		for ( Size jj = 0; jj < dimsizes_[ 2 ]; ++jj ) {
			bin[ 2 ] = jj;
			for ( Size kk = 0; kk < dimsizes_[ 3 ]; ++kk ) {
				bin[ 3 ] = kk;
				utility::fixedsizearray1< Vector, 8 > grid_corners = corners( bin );
				bool any_contained( false );
				for ( Size ll = 1; ll <= 8; ++ll ) {
					Real d2 = center.distance_squared( grid_corners[ ll ] );
					if ( d2 <= rad2 ) { any_contained = true; break; }
				}
				if ( any_contained ) {
					set_value_for_bin( bin, true );
				}
			}
		}
	}
*/

	using namespace core;

	Real const rad2 = radius * radius;
	Real const rad_minus_cendis =  radius - half_bin_width_root_three_;

	Real rad_minus_cendis2 = rad_minus_cendis * rad_minus_cendis;
	Real const rad_plus_cendis = half_bin_width_root_three_ + radius;
	Real const rad_plus_cendis2 = rad_plus_cendis * rad_plus_cendis;

	Bin3D bin;
	Vector relpos = center - bb_extended_.lower();
	SSize xlo( static_cast< SSize > ((relpos.x() - radius) / bin_width_ ));
	SSize xhi( static_cast< SSize > ((relpos.x() + radius) / bin_width_ ) + 1);
	SSize ylo( static_cast< SSize > ((relpos.y() - radius) / bin_width_ ));
	SSize yhi( static_cast< SSize > ((relpos.y() + radius) / bin_width_ ) + 1);
	SSize zlo( static_cast< SSize > ((relpos.z() - radius) / bin_width_ ));
	SSize zhi( static_cast< SSize > ((relpos.z() + radius) / bin_width_ ) + 1);

	Size uxlo = xlo < 0 ? 0 : xlo;
	Size uxhi = xhi < 0 ? 0 : ( xhi > SSize( dimsizes_[ 1 ] ) ? dimsizes_[ 1 ] : xhi ) ;
	Size uylo = ylo < 0 ? 0 : ylo;
	Size uyhi = yhi < 0 ? 0 : ( yhi > SSize( dimsizes_[ 2 ] ) ? dimsizes_[ 2 ] : yhi ) ;
	Size uzlo = zlo < 0 ? 0 : zlo;
	Size uzhi = zhi < 0 ? 0 : ( zhi > SSize( dimsizes_[ 3 ] ) ? dimsizes_[ 3 ] : zhi ) ;

	for ( Size ii = uxlo; ii < uxhi; ++ii ) {
		bin[ 1 ] = ii;
		for ( Size jj = uylo; jj < uyhi; ++jj ) {
			bin[ 2 ] = jj;
			for ( Size kk = uzlo; kk < uzhi; ++kk ) {
				bin[ 3 ] = kk;

				Vector bin_center_point( bin_center( bin ));
				Real sphere_cent_to_voxel_cent2 = bin_center_point.distance_squared( center );

				/// If the distance between the sphere center and the voxel center is greater
				/// than the radius of the sphere plus the distance from the voxel center and
				/// one of it's corners, then it is impossible for any part of this voxel to
				/// overlap with the sphere.  Proceed to the next voxel.
				if (  sphere_cent_to_voxel_cent2 > rad_plus_cendis2 ) continue;

				/// if the voxel center is contained within the sphere, mark this voxel as occupied and move on.
				bool any_contained( sphere_cent_to_voxel_cent2 < rad2 );

				if ( !any_contained && sphere_cent_to_voxel_cent2 > rad_minus_cendis2 ) {
					// Possible that some corners are covered by the sphere, but not guaranteed.
					utility::fixedsizearray1< Vector, 8 > grid_corners = corners( bin );
					for ( Size ll = 1; ll <= 8; ++ll ) {
						Real d2 = center.distance_squared( grid_corners[ ll ] );
						if ( d2 < rad2 ) { any_contained = true; break; }
					}
					if ( ! any_contained ) {
						/// look for a glancing contact between a sphere and one face of this box,
						/// or containment of the sphere within the box
						bool x_in_range = ( center.x() >= grid_corners[ 1 ].x() && center.x() <= grid_corners[ 8 ].x() );
						bool y_in_range = ( center.y() >= grid_corners[ 1 ].y() && center.y() <= grid_corners[ 8 ].y() );
						bool z_in_range = ( center.z() >= grid_corners[ 1 ].z() && center.z() <= grid_corners[ 8 ].z() );
						Size n_in_range = (Size) x_in_range + ( Size ) y_in_range + ( Size ) z_in_range;

						if ( n_in_range == 2 ) {
							/// glancing contact possible.  Distance from sphere center to the voxel face is
							/// the distance along the dimension that's not in range.
							Size dim_out_of_range( 1 );
							if ( ! x_in_range ) {
								dim_out_of_range = 1;
							} else if ( ! y_in_range ) {
								dim_out_of_range = 2;
							} else { //if ( ! z_in_range ) {
								assert( ! z_in_range );
								dim_out_of_range = 3;
							}

							Real dim_edge;
							if ( center( dim_out_of_range ) < grid_corners[ 1 ]( dim_out_of_range )) {
								dim_edge = grid_corners[ 1 ]( dim_out_of_range );
							} else {
								dim_edge = grid_corners[ 8 ]( dim_out_of_range );
							}
							Real d = center( dim_out_of_range ) - dim_edge;
							if ( d*d < rad2 ) {
								any_contained = true;
							}
						}
					}
				}

				if ( any_contained ) {
					set_value_for_bin( bin, true );
				}
			}
		}
	}

}

/// @details The sphere list should describe spheres by the centers and by their
/// square radii, not by their radii.
void
Bool3DGrid::or_by_spheres_conservative(
	utility::vector1< std::pair< Vector, Real > > const & spheres
)
{
	Bin3D bin;
	for ( Size ii = 0; ii < dimsizes_[ 1 ]; ++ii ) {
		bin[ 1 ] = ii;
		for ( Size jj = 0; jj < dimsizes_[ 2 ]; ++jj ) {
			bin[ 2 ] = jj;
			for ( Size kk = 0; kk < dimsizes_[ 3 ]; ++kk ) {
				bin[ 3 ] = kk;
				CornerPoints grid_corners = corners( bin );
				bool all_corners_contained( true );
				for ( Size ll = 1; ll <= 8; ++ll ) {
					bool contained_anywhere( false );
					for ( Size mm = 1; mm <= spheres.size(); ++mm ) {
						Real d2 = spheres[ mm ].first.distance_squared( grid_corners[ ll ] );
						if ( d2 <= spheres[ mm ].second ) { contained_anywhere = true; break; }
					}
					if ( ! contained_anywhere ) {
						all_corners_contained = false;
						break;
					}
				}
				if ( all_corners_contained ) {
					set_value_for_bin( bin, true );
				}
			}
		}
	}

}

//// @brief Turn the values of all the bins that overlap with the
/// volume in this bounding box to true.
void
Bool3DGrid::or_by_box_liberal(
	BoundingBox const & bb
)
{
	Vector clipped_lower( bb.lower() ), clipped_upper( bb.upper() );
	for ( Size ii = 1; ii <= 3; ++ii ) if ( clipped_lower( ii ) < bb_.lower()( ii ) ) clipped_lower( ii ) = bb_.lower()( ii );
	for ( Size ii = 1; ii <= 3; ++ii ) if ( clipped_upper( ii ) > bb_.upper()( ii ) ) clipped_upper( ii ) = bb_.upper()( ii );

	Bin3D const bin_lo = bin_for_point( clipped_lower );
	Bin3D const bin_hi = bin_for_point( clipped_upper );

	Bin3D bin;
	for ( Size ii = bin_lo[ 1 ]; ii <= bin_hi[ 1 ]; ++ii ) {
		bin[ 1 ] = ii;
		for ( Size jj = bin_lo[ 2 ]; jj <= bin_hi[ 2 ]; ++jj ) {
			bin[ 2 ] = jj;
			for ( Size kk = bin_lo[ 3 ]; kk <= bin_hi[ 3 ]; ++kk ) {
				bin[ 3 ] = kk;
				set_value_for_bin( bin, true );
			}
		}
	}
}

/// Accessors
Bool3DGrid::Bin3D Bool3DGrid::dimsizes() const {
	return dimsizes_;
}

Bool3DGrid::CornerPoints
Bool3DGrid::corners( Bin3D const & bin ) const {
	CornerPoints corns;
	Vector lower_corner = bb_.lower();
	lower_corner.x() += bin_width_ * bin[ 1 ];
	lower_corner.y() += bin_width_ * bin[ 2 ];
	lower_corner.z() += bin_width_ * bin[ 3 ];

	//corns[ 1 ] = lower_corner;
	//corns[ 2 ] = lower_corner + bin_width_ * Vector( 0, 0, 1 );
	//corns[ 3 ] = lower_corner + bin_width_ * Vector( 0, 1, 0 );
	//corns[ 4 ] = lower_corner + bin_width_ * Vector( 0, 1, 1 );
	//corns[ 5 ] = lower_corner + bin_width_ * Vector( 1, 0, 0 );
	//corns[ 6 ] = lower_corner + bin_width_ * Vector( 1, 0, 1 );
	//corns[ 7 ] = lower_corner + bin_width_ * Vector( 1, 1, 0 );
	//corns[ 8 ] = lower_corner + bin_width_ * Vector( 1, 1, 1 );

	Real const xlo( lower_corner.x() ), xhi( lower_corner.x() + bin_width_ );
	Real const ylo( lower_corner.y() ), yhi( lower_corner.y() + bin_width_ );
	Real const zlo( lower_corner.z() ), zhi( lower_corner.z() + bin_width_ );

	corns[ 1 ] = Vector( xlo, ylo, zlo );
	corns[ 2 ] = Vector( xlo, ylo, zhi );
	corns[ 3 ] = Vector( xlo, yhi, zlo );
	corns[ 4 ] = Vector( xlo, yhi, zhi );
	corns[ 5 ] = Vector( xhi, ylo, zlo );
	corns[ 6 ] = Vector( xhi, ylo, zhi );
	corns[ 7 ] = Vector( xhi, yhi, zlo );
	corns[ 8 ] = Vector( xhi, yhi, zhi );


	//corns[ 1 ] = lower_corner;
	//corns[ 2 ] = lower_corner; corns[ 2 ].z() += bin_width_; // Vector( 0, 0, 1 );
	//corns[ 3 ] = lower_corner; corns[ 3 ].y() += bin_width_; // Vector( 0, 1, 0 );
	//corns[ 4 ] = lower_corner; corns[ 4 ].y() += bin_width_; corns[ 4 ].z() += bin_width_; // Vector( 0, 1, 1 );
	//corns[ 5 ] = lower_corner; corns[ 5 ].x() += bin_width_; // Vector( 1, 0, 0 );
	//corns[ 6 ] = corns[ 5 ];   corns[ 6 ].z() += bin_width_; // Vector( 1, 0, 1 );
	//corns[ 7 ] = corns[ 5 ];   corns[ 7 ].y() += bin_width_; // Vector( 1, 1, 0 );
	//corns[ 8 ] = corns[ 5 ];   corns[ 8 ].y() += bin_width_; corns[ 8 ].z() += bin_width_;// Vector( 1, 1, 1 );

	return corns;
}

Bool3DGrid::BoundingBox
Bool3DGrid::bin_extrema( Bin3D const & bin ) const
{
	Vector lower_corner = bb_.lower();
	lower_corner.x() += bin_width_ * bin[ 1 ];
	lower_corner.y() += bin_width_ * bin[ 2 ];
	lower_corner.z() += bin_width_ * bin[ 3 ];

	Vector upper_corner = lower_corner + bin_width_;

	return BoundingBox( lower_corner, upper_corner );
}

Bool3DGrid::Vector
Bool3DGrid::bin_center( Bin3D const & bin ) const {

	return Vector(
		bin_width_ * bin[ 1 ] + bin_width_ / 2 ,
		bin_width_ * bin[ 2 ] + bin_width_ / 2 ,
		bin_width_ * bin[ 3 ] + bin_width_ / 2 ) + bb_extended_.lower();
}

bool Bool3DGrid::occupied( Vector const & point ) const
{
	if ( ! bb_extended_.contains( point ) ) return false;

	index_mask_pair imp = index_and_mask_for_point( point );
	//std::cout << "occupied?: " << imp.first << " " << (int) grid_[ imp.first ] << " " << (int) imp.second << std::endl;
	return grid_[ imp.first ] & imp.second;
}

bool Bool3DGrid::occupied( Bin3D const & bin ) const
{
	for ( Size ii = 1; ii <= 3; ++ii ) {
		if ( bin[ ii ] >= dimsizes_[ ii ] ) return false;
	}

	index_mask_pair imp = index_and_mask_for_bin( bin );
	//std::cout << "occupied?: " << imp.first << " " << (int) grid_[ imp.first ] << " " << (int) imp.second << std::endl;
	return grid_[ imp.first ] & imp.second;

}

// Mutators
void Bool3DGrid::set_value_for_bin( Bin3D const & bin, bool setting )
{
	//std::cout << "setting value for bin " << bin[ 1 ] << " " << bin[ 2 ] << " " << bin[ 3 ] << " = " << setting << std::endl;

	Bin3D halfbin;
	halfbin[ 1 ] = bin[ 1 ] / 2; Size xmod2 = bin[ 1 ] % 2;
	halfbin[ 2 ] = bin[ 2 ] / 2; Size ymod2 = bin[ 2 ] % 2;
	halfbin[ 3 ] = bin[ 3 ] / 2; Size zmod2 = bin[ 3 ] % 2;

	Size index = byte_index_from_doublebin( halfbin );
	unsigned char & voxel_bits( grid_[ index ] );

	//unsigned char original_bits( grid_[ index ] ); /// for debugging -- commented out code below.


	if ( setting ) {
		// change only the value of the bit in question
		unsigned char voxbit = mask_from_offsets( xmod2, ymod2, zmod2 );
		voxel_bits |= voxbit;
	} else {
		// keep everything in the voxel bit in it's current state, except the bit in question
		unsigned char mask = negmask_from_offsets( xmod2, ymod2, zmod2 );
		voxel_bits &= mask;
	}

	/*std::cout << "Setting " << setting << " for bin " << bin[ 1 ] << " " << bin[ 2 ] << " " << bin[ 3 ] << " (index " << index << ")";
	std::cout << " with offset " << xmod2 << " " << ymod2 << " " << zmod2 << " original: ";
	for ( Size ii = 0; ii <= 1; ++ii ) {
		for ( Size jj = 0; jj <= 1; ++jj ) {
			for ( Size kk = 0; kk <= 1; ++kk ) {
			  std::cout << ( mask_from_offsets( ii, jj, kk ) & original_bits ? 1 : 0 );
			}
		}
	}
	std::cout << " becomes ";
	for ( Size ii = 0; ii <= 1; ++ii ) {
		for ( Size jj = 0; jj <= 1; ++jj ) {
			for ( Size kk = 0; kk <= 1; ++kk ) {
			  std::cout << ( mask_from_offsets( ii, jj, kk ) & voxel_bits ? 1 : 0  );
			}
		}
	}
	std::cout << std::endl;*/
}

void Bool3DGrid::or_with( Bool3DGrid const & other )
{
	assert( bin_width_ == other.bin_width_ );

	/// 1. Compute the bin overlap
	Vector overlap_low(
		std::max( bb_.lower().x(), other.bb_.lower().x() ),
		std::max( bb_.lower().y(), other.bb_.lower().y() ),
		std::max( bb_.lower().z(), other.bb_.lower().z() ) );

	Vector overlap_high(
		std::min( bb_.upper().x(), other.bb_.upper().x() ),
		std::min( bb_.upper().y(), other.bb_.upper().y() ),
		std::min( bb_.upper().z(), other.bb_.upper().z() ) );

	/// Check for a situation in which the grids have no overlap.
	if ( overlap_low.x() >= overlap_high.x() ||
			overlap_low.y() >= overlap_high.y() ||
			overlap_low.z() >= overlap_high.z() ) return;

	/// Assert that the two grids are compatible: that the region of overlap perfectly
	/// contains all the voxels inside of it (no voxels spill out).
	assert( std::abs( (overlap_low.x() - bb_.lower().x()) - bin_width_ * static_cast< int > ((overlap_low.x() - bb_.lower().x()) / bin_width_ )) < 1e-6 );
	assert( std::abs( (overlap_low.y() - bb_.lower().y()) - bin_width_ * static_cast< int > ((overlap_low.y() - bb_.lower().y()) / bin_width_ )) < 1e-6 );
	assert( std::abs( (overlap_low.z() - bb_.lower().z()) - bin_width_ * static_cast< int > ((overlap_low.z() - bb_.lower().z()) / bin_width_ )) < 1e-6 );

	assert( std::abs( (overlap_low.x() - other.bb_.lower().x()) - bin_width_ * static_cast< int > ((overlap_low.x() - other.bb_.lower().x()) / bin_width_ )) < 1e-6 );
	assert( std::abs( (overlap_low.y() - other.bb_.lower().y()) - bin_width_ * static_cast< int > ((overlap_low.y() - other.bb_.lower().y()) / bin_width_ )) < 1e-6 );
	assert( std::abs( (overlap_low.z() - other.bb_.lower().z()) - bin_width_ * static_cast< int > ((overlap_low.z() - other.bb_.lower().z()) / bin_width_ )) < 1e-6 );

	/// Iterate over the grid indices of this and other in the overlapping region.
	Bin3D n_overlap;
	Real const perturb = bin_width_ / 2;
	n_overlap[ 1 ] = static_cast< Size > (( overlap_high.x() - overlap_low.x() + perturb  ) / bin_width_ );
	n_overlap[ 2 ] = static_cast< Size > (( overlap_high.y() - overlap_low.y() + perturb  ) / bin_width_ );
	n_overlap[ 3 ] = static_cast< Size > (( overlap_high.z() - overlap_low.z() + perturb  ) / bin_width_ );

	Bin3D my_grid_start;
	my_grid_start[ 1 ] = static_cast< Size > (( overlap_low.x() - bb_.lower().x() + perturb ) / bin_width_ );
	my_grid_start[ 2 ] = static_cast< Size > (( overlap_low.y() - bb_.lower().y() + perturb ) / bin_width_ );
	my_grid_start[ 3 ] = static_cast< Size > (( overlap_low.z() - bb_.lower().z() + perturb ) / bin_width_ );

	Bin3D other_grid_start;
	other_grid_start[ 1 ] = static_cast< Size > (( overlap_low.x() - other.bb_.lower().x() + perturb ) / bin_width_ );
	other_grid_start[ 2 ] = static_cast< Size > (( overlap_low.y() - other.bb_.lower().y() + perturb ) / bin_width_ );
	other_grid_start[ 3 ] = static_cast< Size > (( overlap_low.z() - other.bb_.lower().z() + perturb ) / bin_width_ );

	Bin3D my_pos;
	Bin3D other_pos;

	for ( Size ii = 0; ii < n_overlap[ 1 ]; ++ii ) {
		my_pos[    1 ] = my_grid_start[    1 ] + ii;
		other_pos[ 1 ] = other_grid_start[ 1 ] + ii;
		for ( Size jj = 0; jj < n_overlap[ 2 ]; ++jj ) {
			my_pos[    2 ] = my_grid_start[    2 ] + jj;
			other_pos[ 2 ] = other_grid_start[ 2 ] + jj;
			for ( Size kk = 0; kk < n_overlap[ 3 ]; ++kk ) {
				my_pos[    3 ] = my_grid_start[    3 ] + kk;
				other_pos[ 3 ] = other_grid_start[ 3 ] + kk;

				index_mask_pair my_indmask    =       index_and_mask_for_bin( my_pos    );
				index_mask_pair other_indmask = other.index_and_mask_for_bin( other_pos );
				if ( grid_[ my_indmask.first ] & my_indmask.second ) continue;
				if ( other.grid_[ other_indmask.first ] & other_indmask.second ) {
					grid_[ my_indmask.first ] |= my_indmask.second;
				}
			}
		}
	}

}

void Bool3DGrid::and_with( Bool3DGrid const & other )
{
		assert( bin_width_ == other.bin_width_ );

	/// 1. Compute the bin overlap
	Vector overlap_low(
		std::max( bb_.lower().x(), other.bb_.lower().x() ),
		std::max( bb_.lower().y(), other.bb_.lower().y() ),
		std::max( bb_.lower().z(), other.bb_.lower().z() ) );

	Vector overlap_high(
		std::min( bb_.upper().x(), other.bb_.upper().x() ),
		std::min( bb_.upper().y(), other.bb_.upper().y() ),
		std::min( bb_.upper().z(), other.bb_.upper().z() ) );

	/// Check for a situation in which the grids have no overlap.
	if ( overlap_low.x() >= overlap_high.x() ||
			overlap_low.y() >= overlap_high.y() ||
			overlap_low.z() >= overlap_high.z() ) return;

	/// Assert that the two grids are compatible: that the region of overlap perfectly
	/// contains all the voxels inside of it (no voxels spill out).
	assert( std::abs( (overlap_low.x() - bb_.lower().x()) - bin_width_ * static_cast< int > ((overlap_low.x() - bb_.lower().x()) / bin_width_ )) < 1e-6 );
	assert( std::abs( (overlap_low.y() - bb_.lower().y()) - bin_width_ * static_cast< int > ((overlap_low.y() - bb_.lower().y()) / bin_width_ )) < 1e-6 );
	assert( std::abs( (overlap_low.z() - bb_.lower().z()) - bin_width_ * static_cast< int > ((overlap_low.z() - bb_.lower().z()) / bin_width_ )) < 1e-6 );

	assert( std::abs( (overlap_low.x() - other.bb_.lower().x()) - bin_width_ * static_cast< int > ((overlap_low.x() - other.bb_.lower().x()) / bin_width_ )) < 1e-6 );
	assert( std::abs( (overlap_low.y() - other.bb_.lower().y()) - bin_width_ * static_cast< int > ((overlap_low.y() - other.bb_.lower().y()) / bin_width_ )) < 1e-6 );
	assert( std::abs( (overlap_low.z() - other.bb_.lower().z()) - bin_width_ * static_cast< int > ((overlap_low.z() - other.bb_.lower().z()) / bin_width_ )) < 1e-6 );

	/// Iterate over the grid indices of this and other in the overlapping region.
	Bin3D n_overlap;
	Real const perturb = bin_width_ / 2;
	n_overlap[ 1 ] = static_cast< Size > (( overlap_high.x() - overlap_low.x() + perturb  ) / bin_width_ );
	n_overlap[ 2 ] = static_cast< Size > (( overlap_high.y() - overlap_low.y() + perturb  ) / bin_width_ );
	n_overlap[ 3 ] = static_cast< Size > (( overlap_high.z() - overlap_low.z() + perturb  ) / bin_width_ );

	Bin3D my_grid_start;
	my_grid_start[ 1 ] = static_cast< Size > (( overlap_low.x() - bb_.lower().x() + perturb ) / bin_width_ );
	my_grid_start[ 2 ] = static_cast< Size > (( overlap_low.y() - bb_.lower().y() + perturb ) / bin_width_ );
	my_grid_start[ 3 ] = static_cast< Size > (( overlap_low.z() - bb_.lower().z() + perturb ) / bin_width_ );

	Bin3D other_grid_start;
	other_grid_start[ 1 ] = static_cast< Size > (( overlap_low.x() - other.bb_.lower().x() + perturb ) / bin_width_ );
	other_grid_start[ 2 ] = static_cast< Size > (( overlap_low.y() - other.bb_.lower().y() + perturb ) / bin_width_ );
	other_grid_start[ 3 ] = static_cast< Size > (( overlap_low.z() - other.bb_.lower().z() + perturb ) / bin_width_ );

	Bin3D my_pos;
	Bin3D other_pos;

	for ( Size ii = 0; ii < n_overlap[ 1 ]; ++ii ) {
		my_pos[    1 ] = my_grid_start[    1 ] + ii;
		other_pos[ 1 ] = other_grid_start[ 1 ] + ii;
		for ( Size jj = 0; jj < n_overlap[ 2 ]; ++jj ) {
			my_pos[    2 ] = my_grid_start[    2 ] + jj;
			other_pos[ 2 ] = other_grid_start[ 2 ] + jj;
			for ( Size kk = 0; kk < n_overlap[ 3 ]; ++kk ) {
				my_pos[    3 ] = my_grid_start[    3 ] + kk;
				other_pos[ 3 ] = other_grid_start[ 3 ] + kk;

				index_mask_pair my_indmask    =       index_and_mask_for_bin( my_pos    );
				index_mask_pair other_indmask = other.index_and_mask_for_bin( other_pos );
				if ( ! ( grid_[ my_indmask.first ] & my_indmask.second ) ) continue;
				if ( other.grid_[ other_indmask.first ] & other_indmask.second ) continue;
				grid_[ my_indmask.first ] &= ~my_indmask.second;  // set the bit to 0.
			}
		}
	}


}

/// @details Set all the values in this grid to "false" that are true in the other grid.
void Bool3DGrid::subtract( Bool3DGrid const & other )
{
	assert( bin_width_ == other.bin_width_ );

	/// 1. Compute the bin overlap
	Vector overlap_low(
		std::max( bb_.lower().x(), other.bb_.lower().x() ),
		std::max( bb_.lower().y(), other.bb_.lower().y() ),
		std::max( bb_.lower().z(), other.bb_.lower().z() ) );

	Vector overlap_high(
		std::min( bb_.upper().x(), other.bb_.upper().x() ),
		std::min( bb_.upper().y(), other.bb_.upper().y() ),
		std::min( bb_.upper().z(), other.bb_.upper().z() ) );

	/// Check for a situation in which the grids have no overlap.
	if ( overlap_low.x() >= overlap_high.x() ||
			overlap_low.y() >= overlap_high.y() ||
			overlap_low.z() >= overlap_high.z() ) return;

	/// Assert that the two grids are compatible: that the region of overlap perfectly
	/// contains all the voxels inside of it (no voxels spill out).
	assert( std::abs( (overlap_low.x() - bb_.lower().x()) - bin_width_ * static_cast< int > ((overlap_low.x() - bb_.lower().x()) / bin_width_ )) < 1e-6 );
	assert( std::abs( (overlap_low.y() - bb_.lower().y()) - bin_width_ * static_cast< int > ((overlap_low.y() - bb_.lower().y()) / bin_width_ )) < 1e-6 );
	assert( std::abs( (overlap_low.z() - bb_.lower().z()) - bin_width_ * static_cast< int > ((overlap_low.z() - bb_.lower().z()) / bin_width_ )) < 1e-6 );

	assert( std::abs( (overlap_low.x() - other.bb_.lower().x()) - bin_width_ * static_cast< int > ((overlap_low.x() - other.bb_.lower().x()) / bin_width_ )) < 1e-6 );
	assert( std::abs( (overlap_low.y() - other.bb_.lower().y()) - bin_width_ * static_cast< int > ((overlap_low.y() - other.bb_.lower().y()) / bin_width_ )) < 1e-6 );
	assert( std::abs( (overlap_low.z() - other.bb_.lower().z()) - bin_width_ * static_cast< int > ((overlap_low.z() - other.bb_.lower().z()) / bin_width_ )) < 1e-6 );

	/// Iterate over the grid indices of this and other in the overlapping region.
	Bin3D n_overlap;
	Real const perturb = bin_width_ / 2;
	n_overlap[ 1 ] = static_cast< Size > (( overlap_high.x() - overlap_low.x() + perturb  ) / bin_width_ );
	n_overlap[ 2 ] = static_cast< Size > (( overlap_high.y() - overlap_low.y() + perturb  ) / bin_width_ );
	n_overlap[ 3 ] = static_cast< Size > (( overlap_high.z() - overlap_low.z() + perturb  ) / bin_width_ );

	Bin3D my_grid_start;
	my_grid_start[ 1 ] = static_cast< Size > (( overlap_low.x() - bb_.lower().x() + perturb ) / bin_width_ );
	my_grid_start[ 2 ] = static_cast< Size > (( overlap_low.y() - bb_.lower().y() + perturb ) / bin_width_ );
	my_grid_start[ 3 ] = static_cast< Size > (( overlap_low.z() - bb_.lower().z() + perturb ) / bin_width_ );

	Bin3D other_grid_start;
	other_grid_start[ 1 ] = static_cast< Size > (( overlap_low.x() - other.bb_.lower().x() + perturb ) / bin_width_ );
	other_grid_start[ 2 ] = static_cast< Size > (( overlap_low.y() - other.bb_.lower().y() + perturb ) / bin_width_ );
	other_grid_start[ 3 ] = static_cast< Size > (( overlap_low.z() - other.bb_.lower().z() + perturb ) / bin_width_ );

	Bin3D my_pos;
	Bin3D other_pos;

	for ( Size ii = 0; ii < n_overlap[ 1 ]; ++ii ) {
		my_pos[    1 ] = my_grid_start[    1 ] + ii;
		other_pos[ 1 ] = other_grid_start[ 1 ] + ii;
		for ( Size jj = 0; jj < n_overlap[ 2 ]; ++jj ) {
			my_pos[    2 ] = my_grid_start[    2 ] + jj;
			other_pos[ 2 ] = other_grid_start[ 2 ] + jj;
			for ( Size kk = 0; kk < n_overlap[ 3 ]; ++kk ) {
				my_pos[    3 ] = my_grid_start[    3 ] + kk;
				other_pos[ 3 ] = other_grid_start[ 3 ] + kk;

				index_mask_pair my_indmask    =       index_and_mask_for_bin( my_pos    );
				index_mask_pair other_indmask = other.index_and_mask_for_bin( other_pos );
				if ( ! ( grid_[ my_indmask.first ] & my_indmask.second ) ) continue;

				if ( ! ( other.grid_[ other_indmask.first ] & other_indmask.second ) ) continue;

				grid_[ my_indmask.first ] &= ~my_indmask.second;  // set the bit to 0.
			}
		}
	}


}


void Bool3DGrid::clear()
{
	for ( Size ii = 0; ii < grid_.size(); ++ii ) {
		grid_[ ii ] = 0;
	}
}


Bool3DGrid::index_mask_pair
Bool3DGrid::index_and_mask_for_point( Vector const & point ) const
{
	Vector local = point - bb_.lower();
	Bin3D halfbin;
	halfbin[ 1 ] = static_cast< Size > ( local.x() / bin_width_2x_ );
	halfbin[ 2 ] = static_cast< Size > ( local.y() / bin_width_2x_ );
	halfbin[ 3 ] = static_cast< Size > ( local.z() / bin_width_2x_ );
	Size xmod2 = local.x() - halfbin[ 1 ] * bin_width_2x_ >= bin_width_ ? 1 : 0;
	Size ymod2 = local.y() - halfbin[ 2 ] * bin_width_2x_ >= bin_width_ ? 1 : 0;
	Size zmod2 = local.z() - halfbin[ 3 ] * bin_width_2x_ >= bin_width_ ? 1 : 0;

	//std::cout << "index for point " << point.x() << " " << point.y() << " " << point.z() << " ";
	//std::cout << "local point: " << local.x() << " " << local.y() << " " << local.z() << " ";
	//std::cout << "lower bound: " << bb_.lower().x() << " " << bb_.lower().y() << " " << bb_.lower().z() << " halfbin: ";
	//std::cout << halfbin[ 1 ] << " " << halfbin[ 2 ] << " " << halfbin[ 3 ] << " halfdimprods:";
	//std::cout << halfdimprods_[ 1 ] << " " << halfdimprods_[ 2 ] << " " << halfdimprods_[ 3 ] << " " << std::endl;

	Size index = byte_index_from_doublebin( halfbin );
	unsigned char mask = mask_from_offsets( xmod2, ymod2, zmod2 );

	return std::make_pair( index, mask );
}

Bool3DGrid::index_mask_pair
Bool3DGrid::index_and_mask_for_bin( Bin3D const & bin ) const
{
	Bin3D halfbin;
	halfbin[ 1 ] = bin[ 1 ] / 2;
	halfbin[ 2 ] = bin[ 2 ] / 2;
	halfbin[ 3 ] = bin[ 3 ] / 2;
	Size xmod2 = bin[ 1 ] % 2;
	Size ymod2 = bin[ 2 ] % 2;
	Size zmod2 = bin[ 3 ] % 2;

	//std::cout << "index for point " << point.x() << " " << point.y() << " " << point.z() << " ";
	//std::cout << "local point: " << local.x() << " " << local.y() << " " << local.z() << " ";
	//std::cout << "lower bound: " << bb_.lower().x() << " " << bb_.lower().y() << " " << bb_.lower().z() << " halfbin: ";
	//std::cout << halfbin[ 1 ] << " " << halfbin[ 2 ] << " " << halfbin[ 3 ] << " halfdimprods:";
	//std::cout << halfdimprods_[ 1 ] << " " << halfdimprods_[ 2 ] << " " << halfdimprods_[ 3 ] << " " << std::endl;

	Size index = byte_index_from_doublebin( halfbin );
	unsigned char mask = mask_from_offsets( xmod2, ymod2, zmod2 );

	return std::make_pair( index, mask );
}

Bool3DGrid::Bin3D
Bool3DGrid::bin_for_point( Vector const & point ) const
{
	assert( bb_.contains( point ) );
	Vector local = point - bb_.lower();
	Bin3D bin;
	bin[ 1 ] = static_cast< Size > ( local.x() / bin_width_ );
	bin[ 2 ] = static_cast< Size > ( local.y() / bin_width_ );
	bin[ 3 ] = static_cast< Size > ( local.z() / bin_width_ );
	if ( bin[ 1 ] >= dimsizes_[ 1 ] ) bin[ 1 ] = dimsizes_[ 1 ] - 1;
	if ( bin[ 2 ] >= dimsizes_[ 2 ] ) bin[ 2 ] = dimsizes_[ 2 ] - 1;
	if ( bin[ 3 ] >= dimsizes_[ 3 ] ) bin[ 3 ] = dimsizes_[ 3 ] - 1;
	return bin;
}

void
Bool3DGrid::reset_grid() {
	Vector local_upper_corner = bb_.upper() - bb_.lower();
	Bin3D nbins;
	Real const perturb( bin_width_ * 0.5 );
	nbins[ 1 ] = static_cast< Size > ( (local_upper_corner.x() + perturb ) / bin_width_ );
	nbins[ 2 ] = static_cast< Size > ( (local_upper_corner.y() + perturb ) / bin_width_ );
	nbins[ 3 ] = static_cast< Size > ( (local_upper_corner.z() + perturb ) / bin_width_ );

	if ( nbins[ 1 ] * bin_width_ < local_upper_corner.x() ) ++nbins[ 1 ];
	if ( nbins[ 2 ] * bin_width_ < local_upper_corner.y() ) ++nbins[ 2 ];
	if ( nbins[ 3 ] * bin_width_ < local_upper_corner.z() ) ++nbins[ 3 ];

	for ( Size ii = 1; ii <= 3; ++ii ) {
		// make nbins a multiple of the number of bins per super voxel
		// The number of double-bins per supervoxel may be any value, however,
		// the number of bins per supervoxel must be even.
		nbins[ ii ] += nbins[ ii ] % ( 2 * n_doublebins_per_supervoxel ) == 0 ?
			0 : 2 * n_doublebins_per_supervoxel - nbins[ ii ] % ( 2 * n_doublebins_per_supervoxel );
		dimsizes_[ ii ] = nbins[ ii ];
		halfdimsizes_[ ii ] = nbins[ ii ] / 2;
		supervoxel_dimsizes_[ ii ] = halfdimsizes_[ ii ] / n_doublebins_per_supervoxel;
	}


	Vector new_upper( bb_.lower() );
	new_upper.x() += bin_width_ * ( dimsizes_[ 1 ] );
	new_upper.y() += bin_width_ * ( dimsizes_[ 2 ] );
	new_upper.z() += bin_width_ * ( dimsizes_[ 3 ] );

	bb_extended_ = BoundingBox( bb_.lower(), new_upper );

	dimprods_[ 3 ] = 1;
	halfdimprods_[ 3 ] = 1;
	supervoxel_dimprods_[ 3 ] = 1;

	dimprods_[ 2 ] = dimsizes_[ 3 ];
	halfdimprods_[ 2 ] = halfdimsizes_[ 3 ];
	supervoxel_dimprods_[ 2 ] = supervoxel_dimsizes_[ 3 ];

	dimprods_[ 1 ] = dimprods_[ 2 ] * dimsizes_[ 2 ];
	halfdimprods_[ 1 ] = halfdimprods_[ 2 ] * halfdimsizes_[ 2 ];
	supervoxel_dimprods_[ 1 ] = supervoxel_dimprods_[ 2 ] * supervoxel_dimsizes_[ 2 ];

	Size nchar = halfdimprods_[ 1 ] * halfdimsizes_[ 1 ];
	grid_.resize( nchar );
	std::fill( grid_.begin(), grid_.end(), (unsigned char) 0 );

	bin_width_2x_ = bin_width_ * 2;
	half_bin_width_root_three_ = bin_width_ * 0.5 * std::sqrt( 3.0 ); // the distance from the each corner of a voxel to the center.

	//std::cout << "Bool3DGrid Constructor:\n  nbins: " << nbins[ 1 ] << " " << nbins[ 2 ] << " " << nbins[ 3 ] << "\n";
	//std::cout << "  dimprods" << dimprods_[ 1 ] << " " << dimprods_[ 2 ] << " " << dimprods_[ 3 ] << "\n";
	//std::cout << "  halfdimsizes_" << halfdimsizes_[ 1 ] << " " << halfdimsizes_[ 2 ] << " " << halfdimsizes_[ 3 ] << "\n";
	//std::cout << "  halfdimprods_" << halfdimprods_[ 1 ] << " " << halfdimprods_[ 2 ] << " " << halfdimprods_[ 3 ] << "\n";
	//std::cout << "  supervoxel_dimsizes_" << supervoxel_dimsizes_[ 1 ] << " " << supervoxel_dimsizes_[ 2 ] << " " << supervoxel_dimsizes_[ 3 ] << "\n";
	//std::cout << "  supervoxel_dimprods_" << supervoxel_dimprods_[ 1 ] << " " << supervoxel_dimprods_[ 2 ] << " " << supervoxel_dimprods_[ 3 ] << "\n";
	//std::cout << "bin width: " << bin_width_ << " bin width 2x " << bin_width_2x_ << "\n";
	//std::cout << std::endl;
}


unsigned char
Bool3DGrid::mask_from_offsets( Size xmod2, Size ymod2, Size zmod2 ) const
{
	assert( xmod2 == 0 || xmod2 == 1 );
	assert( ymod2 == 0 || ymod2 == 1 );
	assert( zmod2 == 0 || zmod2 == 1 );

	Size offset = 4 * xmod2 + 2 * ymod2 + zmod2;
	switch ( offset ) {
		case 0 : return 0x01;
		case 1 : return 0x02;
		case 2 : return 0x04;
		case 3 : return 0x08;
		case 4 : return 0x10;
		case 5 : return 0x20;
		case 6 : return 0x40;
		case 7 : return 0x80;
	}
	utility_exit_with_message( "ERROR: Bad input to Bool3DGrid mask from offsets" );
	return 0x0;
}

unsigned char
Bool3DGrid::negmask_from_offsets( Size xmod2, Size ymod2, Size zmod2 ) const
{
	assert( xmod2 == 0 || xmod2 == 1 );
	assert( ymod2 == 0 || ymod2 == 1 );
	assert( zmod2 == 0 || zmod2 == 1 );

	Size offset = 4 * xmod2 + 2 * ymod2 + zmod2;
	switch ( offset ) {
		case 0 : return 0xFE;
		case 1 : return 0xFD;
		case 2 : return 0xFB;
		case 3 : return 0xF7;
		case 4 : return 0xEF;
		case 5 : return 0xDF;
		case 6 : return 0xBF;
		case 7 : return 0x7F;
	}
	utility_exit_with_message( "ERROR: Bad input to Bool3DGrid negmask from offsets" );
	return 0x0;

}

Bool3DGrid::Size
Bool3DGrid::byte_index_from_doublebin( Bin3D const & doublebin ) const {

	assert( doublebin[ 1 ] < halfdimsizes_[ 1 ] );
	assert( doublebin[ 2 ] < halfdimsizes_[ 2 ] );
	assert( doublebin[ 3 ] < halfdimsizes_[ 3 ] );

	Size index =
		supervoxel_dimprods_[ 1 ] * ( doublebin[ 1 ] / n_doublebins_per_supervoxel ) +
		supervoxel_dimprods_[ 2 ] * ( doublebin[ 2 ] / n_doublebins_per_supervoxel ) +
		( doublebin[ 3 ] / n_doublebins_per_supervoxel ); // supervoxel index;
	/// must be converted into the starting position for that supervoxel:

	index *= n_doublebins_per_supervoxel * n_doublebins_per_supervoxel * n_doublebins_per_supervoxel;

	index += n_doublebins_per_supervoxel * n_doublebins_per_supervoxel * (doublebin[ 1 ] % n_doublebins_per_supervoxel );
	index += n_doublebins_per_supervoxel * (doublebin[ 2 ] % n_doublebins_per_supervoxel );
	index += doublebin[ 3 ] % n_doublebins_per_supervoxel;

	return index;
}


/////////////////////// BUMP GRID ////////////////////////////////////
BumpGrid::BumpGrid() :
	grids_( n_probe_radii, 0 ),
	pair_permit_overlap_( n_probe_radii ),
	general_permit_overlap_( 0.0 )
{
	for ( Size ii = 1; ii <= n_probe_radii; ++ii ) {
		grids_[ ii ] = protocols::match::Bool3DGridOP( new Bool3DGrid );
		pair_permit_overlap_[ ii ].resize( n_probe_radii );
		std::fill ( pair_permit_overlap_[ ii ].begin(), pair_permit_overlap_[ ii ].end(), 0.0 );
	}

	Real const permit_overlap[ n_probe_radii + 1 ][ n_probe_radii + 1 ] = {
		///  ZERO, H_ARO, H_ALA,   OXY,   NIT, C_CAR, C_ALA, SULPH
		{       0,     0,     0,     0,     0,     0,     0,     0 }, //  ZERO
		{       0,     0,     0,   0.5,   0.5,     0,     0,     0 }, // H_ARO
		{       0,     0,     0,     0,     0,     0,     0,     0 }, // H_ALA
		{       0,   0.5,     0,   0.3,   0.3,     0,     0,     0 }, //   OXY
		{       0,   0.5,     0,   0.3,   0.3,     0,     0,     0 }, //   NIT
		{       0,     0,     0,     0,     0,     0,     0,     0 }, // C_CAR
		{       0,     0,     0,     0,     0,     0,     0,     0 }, // C_ALA
		{       0,     0,     0,     0,     0,     0,     0,     0 }};// SULPH

	for ( Size ii = 1; ii <= n_probe_radii; ++ii ) {
		for ( Size jj = 1; jj <= n_probe_radii; ++jj ) {
			pair_permit_overlap_[ ii ][ jj ] = permit_overlap[ ii ][ jj ];
		}
	}

}

BumpGrid::BumpGrid( BumpGrid const & rhs ) :
	parent(),
	grids_( n_probe_radii, 0 ),
	pair_permit_overlap_( rhs.pair_permit_overlap_ ),
	general_permit_overlap_( rhs.general_permit_overlap_ )
{
	for ( Size ii = 1; ii <= n_probe_radii; ++ii ) {
		grids_[ ii ] = protocols::match::Bool3DGridOP( new Bool3DGrid( *rhs.grids_[ ii ] ) );
	}
}

BumpGrid::~BumpGrid() {}

BumpGrid const & BumpGrid::operator = ( BumpGrid const & rhs )
{
	for ( Size ii = 1; ii <= n_probe_radii; ++ii ) {
		grids_[ ii ] = protocols::match::Bool3DGridOP( new Bool3DGrid( *rhs.grids_[ ii ] ) );
	}
	pair_permit_overlap_ = rhs.pair_permit_overlap_;
	general_permit_overlap_ = rhs.general_permit_overlap_;

	return *this;
}

// @details  This function assumes that the bounding box is large enough
// to contain all of the (real space) spheres.  Make this bounding box large
// enough to hold all of the configuration space spheres -- extend each
// dimension by the largest probe radius.  Then extend each dimension to an
// integer.
void BumpGrid::set_bounding_box( BoundingBox const & bb )
{
	for ( Size ii = 1; ii <= n_probe_radii; ++ii ) {
		Vector lower = bb.lower() - probe_radius( (ProbeRadius) ii  );
		lower.x() = static_cast< Real > ( static_cast< int > ( lower.x() ));
		lower.y() = static_cast< Real > ( static_cast< int > ( lower.y() ));
		lower.z() = static_cast< Real > ( static_cast< int > ( lower.z() ));
		Vector upper = bb.upper() + probe_radius( (ProbeRadius) ii );
		BoundingBox iibb( lower, upper );
		grids_[ ii ] = protocols::match::Bool3DGridOP( new Bool3DGrid );
		grids_[ ii ]->set_bin_width( 0.25 ); /// HARD CODED HACK!
		grids_[ ii ]->set_bounding_box( iibb );
	}

}


/// @details The input bounding box is assumed to hold all of the (real space)
/// spheres that should be represented.  Create a larger bounding box big enough
/// to hold all of the configuration space spheres.
BumpGrid
BumpGrid::create_bump_grid_for_bb( BoundingBox const & bb ) const
{
	BumpGrid stub;

	/// Initialize the new BumpGrid with the sphere overlap parameters held in this bump grid
	stub.pair_permit_overlap_ = pair_permit_overlap_;
	stub.general_permit_overlap_ = general_permit_overlap_;

	stub.set_bounding_box( bb );
	return stub;
}

BumpGridOP
BumpGrid::create_new_bump_grid_for_bb( BoundingBox const & bb ) const {
	BumpGridOP bgop( new BumpGrid );

	/// Initialize the new BumpGrid with the sphere overlap parameters held in this bump grid
	bgop->pair_permit_overlap_ = pair_permit_overlap_;
	bgop->general_permit_overlap_ = general_permit_overlap_;

	bgop->set_bounding_box( bb );
	return bgop;
}

void BumpGrid::or_with( BumpGrid const & other )
{
	for ( Size ii = 1; ii <= n_probe_radii; ++ii ) {
		if ( grids_[ ii ]->actual_bb().intersects( other.grids_[ ii ]->actual_bb() ) ) {
			grids_[ ii ]->or_with( *other.grids_[ ii ] );
		}
	}
}
void BumpGrid::and_with( BumpGrid const & other )
{
	for ( Size ii = 1; ii <= n_probe_radii; ++ii ) {
		if ( grids_[ ii ]->actual_bb().intersects( other.grids_[ ii ]->actual_bb() ) ) {
			grids_[ ii ]->and_with( *other.grids_[ ii ] );
		}
	}
}

/// @details Input sphere is a real-space sphere -- it's radius will be increased by the
/// varying probe radii to create the varying configuration-space spheres.
void BumpGrid::or_by_sphere( Sphere const & input )
{
	for ( Size ii = 1; ii <= n_probe_radii; ++ii ) {
		grids_[ ii ]->or_by_sphere_conservative( input.first, input.second + probe_radius( (ProbeRadius) ii )  );
	}

}

void BumpGrid::or_by_sphere(
	Vector const & center, ProbeRadius radius_type
)
{
	if ( radius_type == ZERO ) return;

	for ( Size ii = 1; ii <= n_probe_radii; ++ii ) {
		grids_[ ii ]->or_by_sphere_conservative( center,
			required_separation_distance( ProbeRadius( ii ), radius_type ) );
	}
}

bool BumpGrid::overlaps( BumpGrid const & other ) const {

	for ( Size ii = 1; ii <= n_probe_radii; ++ii ) {
		if ( grids_[ ii ]->actual_bb().intersects( other.grids_[ ii ]->actual_bb() ) ) {
			return true;
		}
	}
	return false;
}

/// @brief Collision detection by a point p with a given (fixed) radius type.
/// Collision detection is performed in "configuration space", where the obstacles
/// have been convolved with a spherical probe of a given radius.  This collision
/// detection is conservative: it will not report a collision to exist that does not,
/// but may miss a collision that does exist.
bool BumpGrid::occupied( ProbeRadius radius_type, Vector const & p ) const
{
	if ( radius_type == ZERO ) return false;
	if ( grids_[ radius_type ]->actual_bb().contains( p ) ) {
		return grids_[ radius_type ]->occupied( p );
	}

	return false;
}

void BumpGrid::set_general_overlap_tolerance( Real tolerated_overlap )
{
	general_permit_overlap_ = tolerated_overlap;
}


void BumpGrid::set_pair_overlap_tolerance( ProbeRadius rad1, ProbeRadius rad2, Real tolerated_overlap )
{
	pair_permit_overlap_[ rad1 ][ rad2 ] = tolerated_overlap;
	pair_permit_overlap_[ rad2 ][ rad1 ] = tolerated_overlap;
}

BumpGrid::Real
BumpGrid::pair_overlap_tolerance( ProbeRadius rad1, ProbeRadius rad2 ) const
{
	return pair_permit_overlap_[ rad1 ][ rad2 ];
}

utility::vector1< ProbeRadius >
initialize_atomtype_2_probe_radius_map()
{
	using namespace core;
	using namespace core::chemical;

	AtomTypeSetCOP atset = core::chemical::ChemicalManager::get_instance()->atom_type_set( FA_STANDARD );
	utility::vector1< ProbeRadius > attype_2_radtype( atset->n_atomtypes(), ZERO );
	for ( Size ii = 1; ii <= atset->n_atomtypes(); ++ii ) {
		ProbeRadius ii_probe_radius( ZERO );
		AtomType const & iiat( (*atset)[ ii ] );
		if ( iiat.element() == "H" ) {
			if ( iiat.name() == "Hapo" ) ii_probe_radius = H_ALA;
			else ii_probe_radius = H_ARO; // covers Haro, Hpol, HNbb, and HOH
		} else if ( iiat.element() == "N" ) {
			ii_probe_radius = NIT;
		} else if ( iiat.element() == "O" ) {
			ii_probe_radius = OXY;
		} else if ( iiat.element() == "C" ) {
			if ( iiat.name() == "CNH2" || iiat.name()  == "COO" || iiat.name() == "CObb" ) ii_probe_radius = C_CAR;
			else ii_probe_radius = C_ALA;
		} else if ( iiat.element() == "P" || iiat.element() == "S" ) {
			ii_probe_radius = SULPH;
		} else if ( iiat.lj_radius() == 0.0 ) {
			ii_probe_radius = ZERO;
		} else {
			/// find the next smallest probe radius -- this will produce a conservative set of collisions,
			/// as some radii pairs tolerate some collision to describe hydrogen bonds.
			Real radius_shrunk = 0.89 * iiat.lj_radius();
			for ( Size jj = 1; jj <= n_probe_radii; ++jj ) {
				Real rad = BumpGrid::probe_radius( (ProbeRadius) jj );
				if ( rad <= radius_shrunk ) {
					ii_probe_radius = (ProbeRadius) jj;
				} else {
					break;
				}
			}
		}
		attype_2_radtype[ ii ] = ii_probe_radius;
		TR << "Atom type " << iiat.name() << " with ProbeRadius " << ii_probe_radius
			<< " = " << BumpGrid::probe_radius( (ProbeRadius) ii_probe_radius ) << std::endl;
	}
	return attype_2_radtype;
}

ProbeRadius
probe_radius_for_atom_type( core::Size atomtype )
{
	static const utility::vector1< ProbeRadius > attype_2_probe_radius( initialize_atomtype_2_probe_radius_map() );
	return attype_2_probe_radius[ atomtype ];
}

//// Useful non-member functions for BumpGrid.
BumpGridOP
bump_grid_to_enclose_pose( core::pose::Pose const & pose )
{
	using namespace core;

	Vector first_xyz( 0 );
	Real at1rad( 0.0 );

	if ( pose.total_residue() < 1 ) {
		return BumpGridOP( new BumpGrid );
	} else {
		bool found_an_atom( false );
		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			if ( pose.residue( ii ).natoms() > 0 ) {
				first_xyz = pose.residue( ii ).xyz( 1 );
				at1rad = probe_radius_for_atom_type( pose.residue( ii ).atom( 1 ).type() );
 				found_an_atom = true;
				break;
			}
		}
		if ( ! found_an_atom ) {
			/// weird non-empty pose where each residue has 0 atoms?
			return BumpGridOP( new BumpGrid );
		}
	}

	Vector bbl( first_xyz - at1rad ), bbu( first_xyz + at1rad );
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		for ( Size jj = 1; jj <= pose.residue( ii ).natoms(); ++jj ) {
			Real jjrad = BumpGrid::probe_radius( probe_radius_for_atom_type( pose.residue_type( ii ).atom( jj ).atom_type_index() ));
			bbl.min( pose.residue( ii ).xyz( jj ) - jjrad );
			bbu.max( pose.residue( ii ).xyz( jj ) + jjrad );
		}
	}

	BumpGridOP bgop( new BumpGrid );
	numeric::geometry::BoundingBox< Vector > bb( bbl, bbu );
	bgop->set_bounding_box( bb );

	return bgop;
}

BumpGridOP bump_grid_to_enclose_residue_backbone(
	core::conformation::Residue const & residue,
	BumpGrid const & original_grid
)
{
	using namespace core;

	Vector first_xyz( 0 );

	if ( residue.natoms() < 1 ) {
		return BumpGridOP( new BumpGrid );
	} else {
		first_xyz = residue.xyz( 1 );
	}

	Real at1rad = probe_radius_for_atom_type( residue.atom( 1 ).type() );

	core::Vector bbl( first_xyz - at1rad ), bbu( first_xyz + at1rad );
	//core::Vector bbl_strict( first_xyz ), bbu_strict( first_xyz );

	for ( Size ii = 2; ii <= residue.last_backbone_atom(); ++ii ) {
		Real iirad = BumpGrid::probe_radius( probe_radius_for_atom_type( residue.atom( ii ).type() ));
		bbl.min( residue.xyz( ii ) - iirad );
		bbu.max( residue.xyz( ii ) + iirad );
	}

	numeric::geometry::BoundingBox< Vector > bb( bbl, bbu );
	return original_grid.create_new_bump_grid_for_bb( bb );
}

BumpGridOP bump_grid_to_enclose_residue(
	core::conformation::Residue const & residue,
	BumpGrid const & original_grid
)
{
	using namespace core;

	Vector first_xyz( 0 );

	if ( residue.natoms() < 1 ) {
		return BumpGridOP( new BumpGrid );
	} else {
		first_xyz = residue.xyz( 1 );
	}

	Real at1rad = probe_radius_for_atom_type( residue.atom( 1 ).type() );

	core::Vector bbl( first_xyz - at1rad ), bbu( first_xyz + at1rad );
	//core::Vector bbl_strict( first_xyz ), bbu_strict( first_xyz );

	for ( Size ii = 2; ii <= residue.natoms(); ++ii ) {
		Real iirad = BumpGrid::probe_radius( probe_radius_for_atom_type( residue.atom( ii ).type() ));
		bbl.min( residue.xyz( ii ) - iirad );
		bbu.max( residue.xyz( ii ) + iirad );
	}

	numeric::geometry::BoundingBox< Vector > bb( bbl, bbu );

	return original_grid.create_new_bump_grid_for_bb( bb );

}

void
fill_grid_with_residue_spheres(
	core::conformation::Residue const & residue,
	BumpGrid & grid
)
{
	using namespace core;
	for ( Size ii = 1; ii <= residue.natoms(); ++ii ) {
		grid.or_by_sphere( residue.xyz( ii ), probe_radius_for_atom_type( residue.atom( ii ).type() ));
	}
}

void
fill_grid_with_residue_heavyatom_spheres(
	core::conformation::Residue const & residue,
	BumpGrid & grid
)
{
	using namespace core;
	for ( Size ii = 1; ii <= residue.nheavyatoms(); ++ii ) {
		grid.or_by_sphere( residue.xyz( ii ), probe_radius_for_atom_type( residue.atom( ii ).type() ));
	}
}


void
fill_grid_with_backbone_heavyatom_spheres(
	core::conformation::Residue const & residue,
	BumpGrid & grid
)
{
	using namespace core;
	for ( Size ii = 1; ii <= residue.last_backbone_atom(); ++ii ) {
		grid.or_by_sphere( residue.xyz( ii ), probe_radius_for_atom_type( residue.atom( ii ).type() ));
	}
}


Bool3DGridKinemageWriter::Bool3DGridKinemageWriter() :
	unselectable_( true ),
	line_color_( "red" ),
	master_( "grid" ),
	shrink_factor_( 0.75 ),
	skip_completely_buried_positions_( true ),
	write_empty_voxels_( false ),
	empty_voxel_color_( "gray" ),
	write_facets_( false ),
	facet_master_( "grid facets" ),
	facet_color_( "pinktint" ),
	transparent_facets_( false ),
	facet_alpha_( 1.0 )
{}

Bool3DGridKinemageWriter::~Bool3DGridKinemageWriter() {}

void Bool3DGridKinemageWriter::set_unselectable( bool unselectable ) { unselectable_ = unselectable; }
void Bool3DGridKinemageWriter::set_line_color( std::string const & line_color ) { line_color_ = line_color; }
void Bool3DGridKinemageWriter::set_master( std::string const & master ) { master_ = master; }
void Bool3DGridKinemageWriter::set_shrink_factor( Real shrink_factor ) { shrink_factor_ = shrink_factor; }
void Bool3DGridKinemageWriter::set_skip_completely_buried_positions( bool skip_completely_buried_positions ) { skip_completely_buried_positions_ = skip_completely_buried_positions; }

void Bool3DGridKinemageWriter::set_write_empty_voxels( bool write_empty_voxels ) { write_empty_voxels_ = write_empty_voxels; }
void Bool3DGridKinemageWriter::set_empty_voxel_color( std::string const & empty_voxel_color ) { empty_voxel_color_ = empty_voxel_color; }

void Bool3DGridKinemageWriter::set_write_facets( bool write_facets ) { write_facets_ = write_facets; }
void Bool3DGridKinemageWriter::set_facet_master( std::string const & facet_master ) { facet_master_ = facet_master; }
void Bool3DGridKinemageWriter::set_facet_color( std::string const & facet_color ) { facet_color_ = facet_color; }
void Bool3DGridKinemageWriter::set_transparent_facets( bool transparent_facets ) { transparent_facets_ = transparent_facets; }
void Bool3DGridKinemageWriter::set_facet_alpha( Real facet_alpha ) { facet_alpha_ = facet_alpha; }

void Bool3DGridKinemageWriter::write_grid_to_kinemage(
	std::ostream & ostr,
	std::string const & group_name,
	Bool3DGrid const & grid
) const
{
	Bool3DGrid::CornerPoints corner_offsets;

	/// Shrink in from the corners by an offset determined by the "shrink factor" and the bin width.
	corner_offsets[ 1 ] = shrink_factor_ * grid.bin_width() * Vector(  1,  1,  1 );
	corner_offsets[ 2 ] = shrink_factor_ * grid.bin_width() * Vector(  1,  1, -1 );
	corner_offsets[ 3 ] = shrink_factor_ * grid.bin_width() * Vector(  1, -1,  1 );
	corner_offsets[ 4 ] = shrink_factor_ * grid.bin_width() * Vector(  1, -1, -1 );
	corner_offsets[ 5 ] = shrink_factor_ * grid.bin_width() * Vector( -1,  1,  1 );
	corner_offsets[ 6 ] = shrink_factor_ * grid.bin_width() * Vector( -1,  1, -1 );
	corner_offsets[ 7 ] = shrink_factor_ * grid.bin_width() * Vector( -1, -1,  1 );
	corner_offsets[ 8 ] = shrink_factor_ * grid.bin_width() * Vector( -1, -1, -1 );

	Size const A( 1 ), B( 5 ), C( 7 ), D( 3 ), E( 2 ), F( 6 ), G( 8 ), H( 4 );

	ostr << "@group {" << group_name << "} dominant\n";
	ostr << "@vectorlist {} color= " << line_color_ << " width= 1\n";

	Bin3D ndims = grid.dimsizes();
	Bin3D bin;

	ostr.precision( 3 );

	for ( Size ii = 0; ii < ndims[ 1 ]; ++ii ) {
		bin[ 1 ] = ii;
		for ( Size jj = 0; jj < ndims[ 2 ]; ++jj ) {
			bin[ 2 ] = jj;
			for ( Size kk = 0; kk < ndims[ 3 ]; ++kk ) {
				bin[ 3 ] = kk;

				if ( ! grid.occupied( bin ) ) continue;

				bool skip = ( skip_completely_buried_positions_ );
				if ( skip ) {
					if ( bin[ 1 ] == 0 || bin[ 2 ] == 0 || bin[ 3 ] == 0 ) skip = false;
					if ( bin[ 1 ] + 1 >= ndims[ 1 ] || bin[ 2 ] + 1 >= ndims[ 2 ] || bin[ 3 ] + 1 >= ndims[ 3 ] ) skip = false;
				}
				if ( skip ) {
					Bin3D neighb;
					for ( Size ll = 0; ll <= 2; ++ll ) {
						neighb[ 1 ] = bin[ 1 ] - 1 + ll;
						for ( Size mm = 0; mm <= 2; ++mm ) {
							neighb[ 2 ] = bin[ 2 ] - 1 + mm;
							for ( Size nn = 0; nn <= 2; ++nn ) {
								neighb[ 3 ] = bin[ 3 ] - 1 + nn;
								if ( ll == 1 && mm == 1 && nn == 1 ) continue;
								if ( ! grid.occupied( neighb ) ) {
									skip = false;
									break;
								}
							}
							if ( ! skip ) break;
						}
						if ( ! skip ) break;
					}
				}
				/// All 27 of this voxel's neighbors are occupied.
				if ( skip ) continue;

				Bool3DGrid::CornerPoints c = grid.corners( bin );
				for ( Size ll = 1; ll <= 8; ++ll ) c[ ll ] += corner_offsets[ ll ];

				ostr << "{\"}P" << ( unselectable_ ? " U ": " " ) << c[A].x() << " " << c[A].y() << " " << c[A].z() << "\n";
				ostr << "{\"}"  << ( unselectable_ ? " U ": " " ) << c[B].x() << " " << c[B].y() << " " << c[B].z() << "\n";
				ostr << "{\"}"  << ( unselectable_ ? " U ": " " ) << c[C].x() << " " << c[C].y() << " " << c[C].z() << "\n";
				ostr << "{\"}"  << ( unselectable_ ? " U ": " " ) << c[D].x() << " " << c[D].y() << " " << c[D].z() << "\n";
				ostr << "{\"}"  << ( unselectable_ ? " U ": " " ) << c[A].x() << " " << c[A].y() << " " << c[A].z() << "\n";
				ostr << "{\"}"  << ( unselectable_ ? " U ": " " ) << c[E].x() << " " << c[E].y() << " " << c[E].z() << "\n";
				ostr << "{\"}"  << ( unselectable_ ? " U ": " " ) << c[F].x() << " " << c[F].y() << " " << c[F].z() << "\n";
				ostr << "{\"}"  << ( unselectable_ ? " U ": " " ) << c[G].x() << " " << c[G].y() << " " << c[G].z() << "\n";
				ostr << "{\"}"  << ( unselectable_ ? " U ": " " ) << c[H].x() << " " << c[H].y() << " " << c[H].z() << "\n";
				ostr << "{\"}"  << ( unselectable_ ? " U ": " " ) << c[E].x() << " " << c[E].y() << " " << c[E].z() << "\n";
				ostr << "{\"}P" << ( unselectable_ ? " U ": " " ) << c[B].x() << " " << c[B].y() << " " << c[B].z() << "\n";
				ostr << "{\"}"  << ( unselectable_ ? " U ": " " ) << c[F].x() << " " << c[F].y() << " " << c[F].z() << "\n";
				ostr << "{\"}P" << ( unselectable_ ? " U ": " " ) << c[C].x() << " " << c[C].y() << " " << c[C].z() << "\n";
				ostr << "{\"}"  << ( unselectable_ ? " U ": " " ) << c[G].x() << " " << c[G].y() << " " << c[G].z() << "\n";
				ostr << "{\"}P" << ( unselectable_ ? " U ": " " ) << c[D].x() << " " << c[D].y() << " " << c[D].z() << "\n";
				ostr << "{\"}"  << ( unselectable_ ? " U ": " " ) << c[H].x() << " " << c[H].y() << " " << c[H].z() << "\n";

			}
		}
	}

	if ( write_empty_voxels_ ) {
		ostr << "@vectorlist {} color= " << empty_voxel_color_ << " master= {empty voxels} width= 1\n";
		for ( Size ii = 0; ii < ndims[ 1 ]; ++ii ) {
			bin[ 1 ] = ii;
			for ( Size jj = 0; jj < ndims[ 2 ]; ++jj ) {
				bin[ 2 ] = jj;
				for ( Size kk = 0; kk < ndims[ 3 ]; ++kk ) {
					bin[ 3 ] = kk;
					if ( grid.occupied( bin ) ) continue;

					Bool3DGrid::CornerPoints c = grid.corners( bin );
					for ( Size ll = 1; ll <= 8; ++ll ) c[ ll ] += corner_offsets[ ll ];

					ostr << "{\"}P" << ( unselectable_ ? " U ": " " ) << c[A].x() << " " << c[A].y() << " " << c[A].z() << "\n";
					ostr << "{\"}"  << ( unselectable_ ? " U ": " " ) << c[B].x() << " " << c[B].y() << " " << c[B].z() << "\n";
					ostr << "{\"}"  << ( unselectable_ ? " U ": " " ) << c[C].x() << " " << c[C].y() << " " << c[C].z() << "\n";
					ostr << "{\"}"  << ( unselectable_ ? " U ": " " ) << c[D].x() << " " << c[D].y() << " " << c[D].z() << "\n";
					ostr << "{\"}"  << ( unselectable_ ? " U ": " " ) << c[A].x() << " " << c[A].y() << " " << c[A].z() << "\n";
					ostr << "{\"}"  << ( unselectable_ ? " U ": " " ) << c[E].x() << " " << c[E].y() << " " << c[E].z() << "\n";
					ostr << "{\"}"  << ( unselectable_ ? " U ": " " ) << c[F].x() << " " << c[F].y() << " " << c[F].z() << "\n";
					ostr << "{\"}"  << ( unselectable_ ? " U ": " " ) << c[G].x() << " " << c[G].y() << " " << c[G].z() << "\n";
					ostr << "{\"}"  << ( unselectable_ ? " U ": " " ) << c[H].x() << " " << c[H].y() << " " << c[H].z() << "\n";
					ostr << "{\"}"  << ( unselectable_ ? " U ": " " ) << c[E].x() << " " << c[E].y() << " " << c[E].z() << "\n";
					ostr << "{\"}P" << ( unselectable_ ? " U ": " " ) << c[B].x() << " " << c[B].y() << " " << c[B].z() << "\n";
					ostr << "{\"}"  << ( unselectable_ ? " U ": " " ) << c[F].x() << " " << c[F].y() << " " << c[F].z() << "\n";
					ostr << "{\"}P" << ( unselectable_ ? " U ": " " ) << c[C].x() << " " << c[C].y() << " " << c[C].z() << "\n";
					ostr << "{\"}"  << ( unselectable_ ? " U ": " " ) << c[G].x() << " " << c[G].y() << " " << c[G].z() << "\n";
					ostr << "{\"}P" << ( unselectable_ ? " U ": " " ) << c[D].x() << " " << c[D].y() << " " << c[D].z() << "\n";
					ostr << "{\"}"  << ( unselectable_ ? " U ": " " ) << c[H].x() << " " << c[H].y() << " " << c[H].z() << "\n";

				}
			}
		}
	}

	if ( write_facets_ ) {

		for ( Size ii = 0; ii < ndims[ 1 ]; ++ii ) {
			bin[ 1 ] = ii;
			for ( Size jj = 0; jj < ndims[ 2 ]; ++jj ) {
				bin[ 2 ] = jj;
				for ( Size kk = 0; kk < ndims[ 3 ]; ++kk ) {
					bin[ 3 ] = kk;
					if ( ! grid.occupied( bin ) ) continue;

					Bool3DGrid::CornerPoints c = grid.corners( bin );
					for ( Size ll = 1; ll <= 8; ++ll ) c[ ll ] += corner_offsets[ ll ];
						ostr << "@ribbonlist {} color= " << facet_color_ << " master= {" << facet_master_ << "} ";
						if ( transparent_facets_ ) ostr << " alpha= " << facet_alpha_ << " ";
						ostr << "\n";

						ostr << "{\"}"  << ( unselectable_ ? " U ": " " ) << c[H].x() << " " << c[H].y() << " " << c[H].z() << "\n";
						ostr << "{\"}"  << ( unselectable_ ? " U ": " " ) << c[D].x() << " " << c[D].y() << " " << c[D].z() << "\n";
						ostr << "{\"}"  << ( unselectable_ ? " U ": " " ) << c[E].x() << " " << c[E].y() << " " << c[E].z() << "\n";
						ostr << "{\"}"  << ( unselectable_ ? " U ": " " ) << c[A].x() << " " << c[A].y() << " " << c[A].z() << "\n";
						ostr << "{\"}"  << ( unselectable_ ? " U ": " " ) << c[F].x() << " " << c[F].y() << " " << c[F].z() << "\n";
						ostr << "{\"}"  << ( unselectable_ ? " U ": " " ) << c[B].x() << " " << c[B].y() << " " << c[B].z() << "\n";
						ostr << "{\"}"  << ( unselectable_ ? " U ": " " ) << c[G].x() << " " << c[G].y() << " " << c[G].z() << "\n";
						ostr << "{\"}"  << ( unselectable_ ? " U ": " " ) << c[C].x() << " " << c[C].y() << " " << c[C].z() << "\n";
						ostr << "{\"}"  << ( unselectable_ ? " U ": " " ) << c[H].x() << " " << c[H].y() << " " << c[H].z() << "\n";
						ostr << "{\"}"  << ( unselectable_ ? " U ": " " ) << c[D].x() << " " << c[D].y() << " " << c[D].z() << "\n";

						ostr << "@ribonlist {} color= " << facet_color_ << " master= {" << facet_master_ << "} ";
						if ( transparent_facets_ ) ostr << " alpha= " << facet_alpha_ << " ";
						ostr << "\n";

						ostr << "{\"}"  << ( unselectable_ ? " U ": " " ) << c[E].x() << " " << c[E].y() << " " << c[E].z() << "\n";
						ostr << "{\"} " << ( unselectable_ ? " U ": " " ) << c[H].x() << " " << c[H].y() << " " << c[H].z() << "\n";
						ostr << "{\"}"  << ( unselectable_ ? " U ": " " ) << c[F].x() << " " << c[F].y() << " " << c[F].z() << "\n";
						ostr << "{\"}"  << ( unselectable_ ? " U ": " " ) << c[G].x() << " " << c[G].y() << " " << c[G].z() << "\n";

						ostr << "@ribonlist {} color= " << facet_color_ << " master= {" << facet_master_ << "} ";
						if ( transparent_facets_ ) ostr << " alpha= " << facet_alpha_ << " ";
						ostr << "\n";

						ostr << "{\"}"  << ( unselectable_ ? " U ": " " ) << c[A].x() << " " << c[A].y() << " " << c[A].z() << "\n";
						ostr << "{\"}"  << ( unselectable_ ? " U ": " " ) << c[B].x() << " " << c[B].y() << " " << c[B].z() << "\n";
						ostr << "{\"}"  << ( unselectable_ ? " U ": " " ) << c[D].x() << " " << c[D].y() << " " << c[D].z() << "\n";
						ostr << "{\"}"  << ( unselectable_ ? " U ": " " ) << c[C].x() << " " << c[C].y() << " " << c[C].z() << "\n";

				}
			}
		}

	}

}


void Bool3DGridKinemageWriter::write_grid_to_file(
	std::string const & /*fname*/,
	std::string const & /*group_name*/,
	Bool3DGrid const & /*grid*/
) const
{

}

/*
private:
	bool unselectable_;
	std::string line_color_
	std::string master_;
	Real shrink_factor_;


	bool write_empty_voxels_;
	std::string empty_voxel_color_;

	bool write_facets_;
	std::string facet_master_;
	std::string facet_color_;
	bool transparent_facets_;
	Real facet_alpha_;
*/


}
}
