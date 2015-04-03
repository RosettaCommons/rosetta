// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/OccupiedSpaceHash.cc
/// @brief  Declaration for classes in 6D hasher
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

// Unit headers
#include <protocols/match/OccupiedSpaceHash.hh>

// Package headers
#include <protocols/match/VoxelSetIterator.hh>
#include <protocols/match/BumpGrid.hh>

//numeric headers
#include <numeric/geometry/hashing/SixDHasher.hh>

// Utility headers

// C++ headers
#include <fstream>

#include <utility/vector1.hh>


namespace protocols {
namespace match {

OccupiedSpaceHash::OccupiedSpaceHash() : initialized_( false ), revision_id_( 1 ) {}
OccupiedSpaceHash::~OccupiedSpaceHash() {}

void
OccupiedSpaceHash::set_bounding_box(
	BoundingBox const & bb
)
{
	assert( ! initialized_ );
	bb_ = bb;
}

void
OccupiedSpaceHash::set_uniform_xyz_bin_width( Real bin_width )
{
	assert( ! initialized_ );
	assert( bin_width > 0 );
	Real half_witdh = 0.5 * bin_width;
	for ( Size ii = 1; ii <= 3; ++ii ) {
		xyz_bin_widths_[ ii ] = bin_width;
		xyz_bin_halfwidths_[ ii ] = half_witdh;
	}
}

void
OccupiedSpaceHash::set_uniform_euler_angle_bin_width( Real bin_width_degrees )
{
	assert( ! initialized_ );
	assert( bin_width_degrees > 0 );
	Real half_width = 0.5 * bin_width_degrees;
	for ( Size ii = 1; ii <= 3; ++ii ) {
		euler_bin_widths_[ ii ] = bin_width_degrees;
		euler_bin_halfwidths_[ ii ] = half_width;
	}
}

void
OccupiedSpaceHash::set_xyz_bin_widths( Vector const & bin_widths )
{
	assert( ! initialized_ );
	for ( Size ii = 1; ii <= 3; ++ii ) {
		assert( bin_widths( ii ) > 0 );
		xyz_bin_widths_[ ii ] = bin_widths( ii );
		xyz_bin_halfwidths_[ ii ] = 0.5 * bin_widths( ii );
	}

}

void
OccupiedSpaceHash::set_euler_bin_widths( Vector const & euler_bin_widths )
{
	assert( ! initialized_ );
	for ( Size ii = 1; ii <= 3; ++ii ) {
		assert( euler_bin_widths( ii ) > 0 );
		euler_bin_widths_[ ii ] = euler_bin_widths( ii );
		euler_bin_halfwidths_[ ii ] = 0.5 * euler_bin_widths( ii );
	}
}

//void
//OccupiedSpaceHash::set_expected_hash_count()
//{
//}

void
OccupiedSpaceHash::initialize()
{
	assert( ! initialized_ );

	initialized_ = true;

	for ( Size ii = 1; ii <= 3; ++ii ) {
		n_xyz_bins_[ ii ] = static_cast< Size > (( bb_.upper()(ii) - bb_.lower()(ii) ) / xyz_bin_widths_[ ii ] );
		if ( bb_.lower()(ii) + n_xyz_bins_[ ii ]*xyz_bin_widths_[ ii ] < bb_.upper()( ii ) ) {
			++n_xyz_bins_[ ii ];
		}
	}

	for ( Size ii = 1; ii <= 3; ++ii ) {
		Real range = ii == 3 ? 180 : 360;
		n_euler_bins_[ ii ] = static_cast< Size > ( range / euler_bin_widths_[ ii ] );
		if ( n_euler_bins_[ ii ] * euler_bin_widths_[ ii ] < range ) {
			n_euler_bins_[ ii ];
		}
	}

	dim_prods_[ 6 ] = 1;
	dim_prods_[ 5 ] = n_euler_bins_[ 3 ];
	dim_prods_[ 4 ] = n_euler_bins_[ 2 ] * dim_prods_[ 5 ];
	dim_prods_[ 3 ] = n_euler_bins_[ 1 ] * dim_prods_[ 4 ];
	dim_prods_[ 2 ] = n_xyz_bins_[ 3 ]   * dim_prods_[ 3 ];
	dim_prods_[ 1 ] = n_xyz_bins_[ 2 ]   * dim_prods_[ 2 ];

	threeD_projection_ = Bool3DGridOP( new Bool3DGrid );
	threeD_projection_->set_bounding_box( bb_ );
	threeD_projection_->set_bin_width( xyz_bin_widths_[ 1 ] );
}

void
OccupiedSpaceHash::insert_hit_geometry( Real6 const & geom )
{
	assert( initialized_ );

	if ( revision_id_ == 1 ) { ++revision_id_; } /// round 1 is over.

	Vector point( Vector( geom[ 1 ], geom[ 2 ], geom[ 3 ] ));
	if ( ! bb_.contains( point ) ) return;

	VoxelSetIterator voxiter( bb_, n_xyz_bins_, n_euler_bins_, xyz_bin_widths_,
		euler_bin_widths_, xyz_bin_halfwidths_, euler_bin_halfwidths_, geom );

	numeric::geometry::hashing::Bin6D bin;
	Size pos;
	boost::uint64_t bin_index;
	while ( ! voxiter.at_end() ) {
		voxiter.get_bin_and_pos( bin, pos );
		bin_index = calc_bin_index( bin );
		ActiveVoxelSet::iterator iter = hash_.find( bin_index );
		if ( iter == hash_.end() ) {
			hash_.insert( std::make_pair( bin_index, bitmask_for_position( pos ) ));
		} else {
			iter->second |= bitmask_for_position( pos );
		}
		++voxiter;
	}

	project_point_to_3d( geom );
}

void
OccupiedSpaceHash::prepare_to_note_hits_for_completed_round()
{
	/// Increment the revision id.  This signals to observing ClassicMatchAlgorithms that
	/// hits generated in previous rounds might no-longer be viable; after this round
	/// any hit such that none of it's 64 voxels recieved any hits from this round is unable
	/// to yield a match.
	///
	/// Rounds that use secondary matching on upstream residues will not
	++revision_id_;

	for ( ActiveVoxelSet::iterator iter = hash_.begin(), iter_end = hash_.end();
			iter != iter_end; ++iter ) {
		iter->second = 0;
	}

	reset_3d_projection();
}

void
OccupiedSpaceHash::reset_3d_projection()
{
	threeD_projection_->clear();
}


bool
OccupiedSpaceHash::previous_round_geometry_still_matchable( Real6 const & geom )
{
	bool still_matchable( match_possible_for_hit_geometry( geom ));
	if ( still_matchable ) {
		project_point_to_3d( geom );
	}
	return still_matchable;
}

void
OccupiedSpaceHash::note_hit_geometry( Real6 const & geom )
{
	assert( initialized_ );
	Vector point( Vector( geom[ 1 ], geom[ 2 ], geom[ 3 ] ));
	if ( ! bb_.contains( point ) ) return;

	VoxelSetIterator voxiter( bb_, n_xyz_bins_, n_euler_bins_, xyz_bin_widths_,
		euler_bin_widths_, xyz_bin_halfwidths_, euler_bin_halfwidths_, geom );

	numeric::geometry::hashing::Bin6D bin;
	Size pos;
	boost::uint64_t bin_index;
	while ( ! voxiter.at_end() ) {
		voxiter.get_bin_and_pos( bin, pos );
		bin_index = calc_bin_index( bin );
		ActiveVoxelSet::iterator iter = hash_.find( bin_index );
		if ( iter != hash_.end() ) {
			boost::uint64_t mask = bitmask_for_position( pos );
			if ( ! (iter->second & mask) ) {
				iter->second |= mask;
			}
		}
		++voxiter;
	}

	project_point_to_3d( geom );
}

/// @details It must be safe for multiple threads to call function simultaneously
bool
OccupiedSpaceHash::match_possible_for_hit_geometry( Real6 const & geom ) const
{
	using namespace utility;

	assert( initialized_ );
	Vector point( Vector( geom[ 1 ], geom[ 2 ], geom[ 3 ] ));
	if ( ! bb_.contains( point ) ) return false;

	/// Note, in the first round, the occupied space hash is empty.  Hits have not yet
	/// been inserted into it.  Therefore, the purpose of the check in the first round
	/// is to ensure the hit lies inside the bounding box.
	if ( revision_id_ == 1 ) return true;

	if ( ! threeD_projection_->occupied( point ) ) return false;

	VoxelSetIterator voxiter( bb_, n_xyz_bins_, n_euler_bins_, xyz_bin_widths_,
		euler_bin_widths_, xyz_bin_halfwidths_, euler_bin_halfwidths_, geom );

	numeric::geometry::hashing::Bin6D bin;
	Size pos;
	boost::uint64_t bin_index;
	while ( ! voxiter.at_end() ) {
		voxiter.get_bin_and_pos( bin, pos );
		bin_index = calc_bin_index( bin );
		ActiveVoxelSet::const_iterator iter = hash_.find( bin_index );
		if ( iter != hash_.end() && iter->second & bitmask_for_position( pos ) ) {
			return true;
		}
		++voxiter;
	}
	return false;

/*	for ( Size ii = 1; ii <= N_HASH_MAPS; ++ii ) {
		if ( hashes_[ ii ].first->contains( geom ) ) {
			boost::uint64_t bin_index = hashes_[ ii ].first->bin_index( geom );
			ActiveVoxelSet::const_iterator iter = hashes_[ ii ].second.find( bin_index );
			if ( iter != hashes_[ ii ].second.end() ) {
				//if ( ! threeD_projection_->occupied( point ) ) {
				//	std::cout << "point is not active?" << point.x() << " " << point.y() << " " << point.z() << std::endl;
				//}
				return true;
			}
		}
	}

	/// None of the 64 bins was already occupied by this geometry; a match cannot be found with this hit!
	return false;*/
}

void
OccupiedSpaceHash::drop_unsatisfied_voxels()
{
	assert( initialized_ );
	//threeD_projection_->clear();

	for ( ActiveVoxelSet::iterator
			iter = hash_.begin(),
			iter_end = hash_.end();
			iter != iter_end; /* no increment */ ) {
		ActiveVoxelSet::iterator iter_next = iter;
		++iter_next;
		if ( iter->second == 0 ) {
			hash_.erase( iter );
		}
		iter = iter_next;
	}

}

void
OccupiedSpaceHash::write_3Dprojection_kinemage( std::string const & kin_file_name )
{
	Bool3DGridKinemageWriter writer;
	writer.set_skip_completely_buried_positions( true );
	std::ofstream fout;
	fout.open( kin_file_name.c_str() );
	writer.write_grid_to_kinemage( fout, "occspace", *threeD_projection_ );
	fout.close();
}

OccupiedSpaceHash::Size
OccupiedSpaceHash::revision_id() const
{
	return revision_id_;
}


void
OccupiedSpaceHash::project_point_to_3d( Real6 const & geom )
{
/*
	Vector point( Vector( geom[ 1 ], geom[ 2 ], geom[ 3 ] ));

	Vector lower_bound( point ), upper_bound( point );
	for ( Size ii = 1; ii <= 3; ++ii ) lower_bound( ii ) -= xyz_width_root3_[ ii ];
	for ( Size ii = 1; ii <= 3; ++ii ) if ( lower_bound( ii ) < bb_.lower()( ii ) ) lower_bound( ii ) = bb_.lower()( ii );
	for ( Size ii = 1; ii <= 3; ++ii ) upper_bound( ii ) += xyz_width_root3_[ ii ];
	for ( Size ii = 1; ii <= 3; ++ii ) if ( upper_bound( ii ) < bb_.upper()( ii ) ) upper_bound( ii ) = bb_.upper()( ii );

	threeD_projection_->or_by_box_liberal( BoundingBox( lower_bound, upper_bound ) );
*/
	Vector point( Vector( geom[ 1 ], geom[ 2 ], geom[ 3 ] ));
	assert( bb_.contains( point ));

	Vector local = point - bb_.lower();
	Vector lower = bb_.lower();
	Vector upper = bb_.lower();

	Vector local_lower;
	for ( Size ii = 1; ii <= 3; ++ii ) {
		local_lower( ii ) = static_cast< Size > ( local( ii ) / xyz_bin_widths_[ ii ] ) * xyz_bin_widths_[ ii ];
	}
	lower += local_lower;
	upper += local_lower;
	for ( Size ii = 1; ii <= 3; ++ii ) {
		upper( ii ) += xyz_bin_widths_[ ii ];
	}

	for ( Size ii = 1; ii <= 3; ++ii ) {
		if ( point( ii ) - lower( ii ) > 0.5 ) {
			upper( ii ) += 0.5 * xyz_bin_widths_[ ii ];
		} else {
			lower( ii ) -= 0.5 * xyz_bin_widths_[ ii ];
		}
	}
	for ( Size ii = 1; ii <= 3; ++ii ) if ( lower( ii ) < bb_.lower()( ii ) ) lower( ii ) = bb_.lower()( ii );
	for ( Size ii = 1; ii <= 3; ++ii ) if ( upper( ii ) > bb_.upper()( ii ) ) upper( ii ) = bb_.upper()( ii );


	threeD_projection_->or_by_box_liberal( BoundingBox( lower, upper ) );

}

/// @details return the 64-bit bitmask for a particular bit -- numbered in increasing order
/// from right to left.
boost::uint64_t
OccupiedSpaceHash::bitmask_for_position( Size pos ) const {
	assert( pos > 0 && pos <= 64 );

	switch ( pos ) {
		case  1: return 0x01;
		case  2: return 0x02;
		case  3: return 0x04;
		case  4: return 0x08;
		case  5: return 0x10;
		case  6: return 0x20;
		case  7: return 0x40;
		case  8: return 0x80;

		case  9: return 0x0100;
		case 10: return 0x0200;
		case 11: return 0x0400;
		case 12: return 0x0800;
		case 13: return 0x1000;
		case 14: return 0x2000;
		case 15: return 0x4000;
		case 16: return 0x8000;

		case 17: return 0x010000;
		case 18: return 0x020000;
		case 19: return 0x040000;
		case 20: return 0x080000;
		case 21: return 0x100000;
		case 22: return 0x200000;
		case 23: return 0x400000;
		case 24: return 0x800000;

		case 25: return 0x01000000;
		case 26: return 0x02000000;
		case 27: return 0x04000000;
		case 28: return 0x08000000;
		case 29: return 0x10000000;
		case 30: return 0x20000000;
		case 31: return 0x40000000;
		case 32: return 0x80000000;

		case 33: return 0x0100000000LL;
		case 34: return 0x0200000000LL;
		case 35: return 0x0400000000LL;
		case 36: return 0x0800000000LL;
		case 37: return 0x1000000000LL;
		case 38: return 0x2000000000LL;
		case 39: return 0x4000000000LL;
		case 40: return 0x8000000000LL;

		case 41: return 0x010000000000LL;
		case 42: return 0x020000000000LL;
		case 43: return 0x040000000000LL;
		case 44: return 0x080000000000LL;
		case 45: return 0x100000000000LL;
		case 46: return 0x200000000000LL;
		case 47: return 0x400000000000LL;
		case 48: return 0x800000000000LL;

		case 49: return 0x01000000000000LL;
		case 50: return 0x02000000000000LL;
		case 51: return 0x04000000000000LL;
		case 52: return 0x08000000000000LL;
		case 53: return 0x10000000000000LL;
		case 54: return 0x20000000000000LL;
		case 55: return 0x40000000000000LL;
		case 56: return 0x80000000000000LL;

		case 57: return 0x0100000000000000LL;
		case 58: return 0x0200000000000000LL;
		case 59: return 0x0400000000000000LL;
		case 60: return 0x0800000000000000LL;
		case 61: return 0x1000000000000000LL;
		case 62: return 0x2000000000000000LL;
		case 63: return 0x4000000000000000LL;
		case 64: return 0x8000000000000000LL;
		default: break;
	}
	utility_exit_with_message( "Critical Error in OccupiedSpaceHash::bitmask_for_position. position outside of range [1 .. 64].  Cannot continue." );
	return 0;
}

boost::uint64_t
OccupiedSpaceHash::calc_bin_index(  numeric::geometry::hashing::Bin6D const & bin ) const
{
	boost::uint64_t index( 0 );
	for ( Size ii = 1; ii <= 6; ++ii ) {
		index += bin[ ii ] * dim_prods_[ ii ];
	}
	return index;
}


}
}

