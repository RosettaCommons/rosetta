// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/OccupiedSpaceHash.hh
/// @brief  Declaration for classes in 6D hasher
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

#ifndef INCLUDED_protocols_match_OccupiedSpaceHash_hh
#define INCLUDED_protocols_match_OccupiedSpaceHash_hh

// Unit headers
#include <protocols/match/OccupiedSpaceHash.fwd.hh>

// Package headers
#include <protocols/match/BumpGrid.fwd.hh>

// project headers
#include <core/types.hh>

// Numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/geometry/BoundingBox.hh>
#include <numeric/geometry/hashing/SixDHasher.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/fixedsizearray1.hh>


/// Boost headers
#include <boost/unordered_map.hpp>

namespace protocols {
namespace match {


/// @brief This class keeps track of the active voxels "up until now" with 64 hashes.
///
/// @details
/// Each hash has a slightly shifted definition of the origin, so that points which
/// are close in 6D will end up in the same hash voxel in at least one of these 64 hashes.
/// The amount of data held in this table is intentionally small: a single boolean (32 bits).
/// The meaning is as follows:
/// In the first round of matching constraints, "false" is inserted into each of the hashes
/// for each non-colliding placement of the downstream partner.
/// The map is not altered between the end of the first round and the beginning of the second round.
/// At the end of all rounds besides the first round, the Matcher updates the OccupiedSpaceHash
/// object for each non-colliding placement of the dowstream partner found during that round:
/// for each non-colliding placement of the downstream partner and for
/// each of the 64 hashes, the 6Dofs describing the downstream partner are hashed: if the
/// hash finds an element in the hash table, that element's value is set to "true".
/// At the end of each round past the second round, the entire hash is traversed.  Any element
/// whose value is "false" represents a voxel that failed to find a match during that round. Such
/// elements are deleted from the hash.  Any element whose value is "true" DID get a hit from that
/// round and should remain in the hash.  The element's value is set to "false".
/// After the OccupiedSpaceHash has cleared out the elements which failed to find a match for the
/// previous round, the Matcher may proceed to "clean up" its older hits (from rounds previous to this one)
/// by querying them against the hash map.  Any hit in the Matcher that does not have a corresponding element
/// in any of the 64 hashes of the OccupiedSpaceHash may be deleted, since it cannot form
/// a complete match.
///
/// This class is intended to be accessed by multiple threads but in a controlled way:
/// read access function "match_possible_for_hit_geometry" is accessible to all threads
/// during the "building" stage. However, the functions note_hit_geometry and
/// drop_unsatisfied_voxels should be called by single thread only, and no other thread
/// should try to access the ActiveVoxelHashes object at that time.
class OccupiedSpaceHash : public utility::pointer::ReferenceCount {
public:
	typedef core::Real                               Real;
	typedef core::Size                               Size;
	typedef core::Vector                             Vector;
	typedef numeric::geometry::BoundingBox< Vector > BoundingBox;
	typedef numeric::geometry::hashing::Real6        Real6;
	typedef numeric::geometry::hashing::Size3        Size3;
	typedef numeric::geometry::hashing::Real3        Real3;
	typedef boost::unordered_map< boost::uint64_t, boost::uint64_t, numeric::geometry::hashing::bin_index_hasher > ActiveVoxelSet;

public:
	OccupiedSpaceHash();
	virtual ~OccupiedSpaceHash();

	void
	set_bounding_box(
		BoundingBox const & bb
	);

	void
	set_uniform_xyz_bin_width( Real bin_width );

	void
	set_uniform_euler_angle_bin_width( Real bin_width_degrees );

	void
	set_xyz_bin_widths( Vector const & bin_widths );

	void
	set_euler_bin_widths( Vector const & euler_bin_widths );

	//void
	//set_expected_hash_count();

	void
	initialize();

	void
	insert_hit_geometry( Real6 const & geom );

	void
	prepare_to_note_hits_for_completed_round();

	void
	note_hit_geometry( Real6 const & );

	void
	reset_3d_projection();

	bool
	previous_round_geometry_still_matchable( Real6 const & );

	bool
	match_possible_for_hit_geometry( Real6 const & ) const;

	void
	drop_unsatisfied_voxels();

	void write_3Dprojection_kinemage( std::string const & kin_file_name );

	Size
	revision_id() const;

private:

	void
	project_point_to_3d( Real6 const & geom );

	boost::uint64_t
	bitmask_for_position( Size pos ) const;

	boost::uint64_t
	calc_bin_index( numeric::geometry::hashing::Bin6D const & bin ) const;

private:

	BoundingBox bb_;

	bool initialized_;

	Size3 n_xyz_bins_;
	Size3 n_euler_bins_;
	utility::fixedsizearray1< boost::uint64_t, 6 > dim_prods_;

	utility::fixedsizearray1< Real, 3 > xyz_bin_widths_;
	utility::fixedsizearray1< Real, 3 > euler_bin_widths_;

	utility::fixedsizearray1< Real, 3 > xyz_bin_halfwidths_;
	utility::fixedsizearray1< Real, 3 > euler_bin_halfwidths_;

	//utility::fixedsizearray1< Real, 3 > xyz_width_root3_;

	ActiveVoxelSet hash_;

	Bool3DGridOP threeD_projection_;

	Size revision_id_;

};

}
}

#endif
