// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/OccupiedSpaceHash.hh
/// @brief  Declaration for classes in 6D hasher
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

#ifndef INCLUDED_protocols_match_MatchSet_hh
#define INCLUDED_protocols_match_MatchSet_hh

// Unit headers
#include <protocols/match/MatchSet.fwd.hh>

// Package headers
#include <protocols/match/Hit.fwd.hh>

//numeric headers
#include <numeric/geometry/hashing/SixDHasher.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/fixedsizearray1.hh>
#include <utility/FixedSizeLexicographicalIterator.fwd.hh>


/// Boost headers
#include <boost/unordered_map.hpp>

/// C++ headers
#include <list>

#include <utility/vector1_bool.hh>


namespace protocols {
namespace match {


/// @brief This object hashes hits into 6D voxels.  This hash can then be traversed
/// to retrieve the hits that hash to the same voxel (matches!).  There are 64 hashes
/// representing the 2^6 ways to perturb the bins in 6D by 1/2 of their bin width.
///
/// @details The hit hasher expects someone else to own the hits.  It takes as input
/// constant pointers to the hits that exist and uses their addresses to hash upon.
/// The hit hasher should only be used if you can guarantee that the hits it points
/// to will outlive the hasher.
class HitHasher : public utility::pointer::ReferenceCount {
public:
	typedef core::Real                               Real;
	typedef core::Size                               Size;
	typedef core::Vector                             Vector;
	typedef numeric::geometry::BoundingBox< Vector > BoundingBox;
	typedef numeric::geometry::hashing::Real3        Real3;
	typedef utility::vector1< std::list< Hit const * > >     MatchSet;
	typedef boost::unordered_map< boost::uint64_t, MatchSet, numeric::geometry::hashing::bin_index_hasher > HitHash;

public:
	HitHasher();
	virtual ~HitHasher();

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

	void
	set_nhits_per_match( Size num_geometric_constraints );

	void
	initialize();

	/// @brief Insert a hits into all 64 hash maps
	void
	insert_hit( Size geometric_constraint_id, Hit const * hit );

	/// @brief Insert a hits into a particular hash maps
	void
	insert_hit( Size which_hash_map, Size geometric_constraint_id, Hit const * hit );

	void
	clear_hash_map( Size which_hash_map );

	HitHash::const_iterator
	hit_hash_begin( Size which_hash_map ) const;

	HitHash::const_iterator
	hit_hash_end( Size which_hash_map ) const;

	numeric::geometry::hashing::SixDCoordinateBinner const &
	binner( Size which_hash_map ) const {
		return * hit_hashes_[ which_hash_map ].first;
	}

private:
	static Size const N_HASH_MAPS = 64; // 2^6 -- all combinations of offsets for 6 dimensions

	BoundingBox bb_;

	bool initialized_;

	utility::fixedsizearray1< Real, 3 > xyz_bin_widths_;
	utility::fixedsizearray1< Real, 3 > euler_bin_widths_;
	Size n_geometric_constraints_per_match_;

	utility::vector1< std::pair< numeric::geometry::hashing::SixDCoordinateBinnerOP, HitHash > > hit_hashes_;

};

/// Class for finding hit neighbors in 6D considering all 64 origin definitions (but without
/// forming all 64 hashes).
class HitNeighborFinder : public utility::pointer::ReferenceCount
{
public:
	typedef core::Real                               Real;
	typedef core::Size                               Size;
	typedef core::Vector                             Vector;
	typedef numeric::geometry::BoundingBox< Vector > BoundingBox;
	typedef numeric::geometry::hashing::Bin6D        Bin6D;

	typedef std::list< Hit const * >                                                 HitPtrList;
	typedef std::list< std::pair< core::Size, Hit const * > >                        HitIndexList;
	typedef boost::unordered_map< boost::uint64_t, HitIndexList, numeric::geometry::hashing::bin_index_hasher >  HitHash;

public:
	HitNeighborFinder();
	virtual ~HitNeighborFinder();

	/// @brief Use the same bounding box as the HitHasher / OccupiedSpaceHash
	void
	set_bounding_box(
		BoundingBox const & bb
	);

	/// @brief Give the same xyz bin witdh given to the HitHasher / OccupiedSpaceHash.
	void
	set_uniform_xyz_bin_width( Real bin_width );

	/// @brief Give the same euler-angle bin witdh given to the HitHasher / OccupiedSpaceHash.
	void
	set_uniform_euler_angle_bin_width( Real bin_width_degrees );


	/// @brief Give the same xyz bin witdhs given to the HitHasher / OccupiedSpaceHash.
	void
	set_xyz_bin_widths( Vector const & bin_widths );

	/// @brief Give the same euler-angle bin witdhs given to the HitHasher / OccupiedSpaceHash.
	void
	set_euler_bin_widths( Vector const & euler_bin_widths );

	/// @brief Call this after the bounding-box and the bin-widths have been set up.  Must
	/// be called before "add_hits" it called.  This initializes the SixDCoordinateBinner.
	void
	initialize();


	/// @brief Add all hits (using hit pointers!) from one list of hits
	void add_hits( std::list< Hit > const & hitlist );

	/// @brief Add hit pointers from the input list
	void add_hits( HitPtrList const & hitptrlist );

	/// @brief Compute the set of connected components for the set of input hits.
	/// This search iterates across the 2^6 halfbin neighbors of each hit and creates
	/// an adjacency graph to determine the connected components (CCs).  The CCs are
	/// returned in a vector where each element contains a list of hit pointers to the
	/// elements of each connected component.
	utility::vector1< HitPtrList >
	connected_components() const;


	/// @brief Find the neighbors of the given set of query hits.  This search iterates
	/// across both the upper and the lower neighbors of the query hits (3^6 neighbors).
	HitPtrList
	neighbor_hits( HitPtrList const & queryhits ) const;

private:
	/// @brief Put the hits that have been added to "all_hits" into the hash_
	void hash_hits();

	Size
	find_next_bin(
		Real6 orig_point,
		Real6 offsets,
		utility::FixedSizeLexicographicalIterator< 6 > const & halfbin_lex,
		Bin6D & next_bin
	) const;

	bool
	within_reach(
		Bin6D const & query_halfbin,
		Bin6D const & nb_halfbin,
		Bin6D const & nbbin,
		utility::FixedSizeLexicographicalIterator< 6 > const & halfbin_lex
	) const;


private:

	BoundingBox bb_;

	bool initialized_;

	utility::fixedsizearray1< Real, 3 > xyz_bin_widths_;
	utility::fixedsizearray1< Real, 3 > euler_bin_widths_;

	HitIndexList all_hits_;

	numeric::geometry::hashing::SixDCoordinateBinnerOP binner_;
	HitHash hash_;
};

/// Class for counting the number of matches given a particular discretization level.
class MatchCounter : public utility::pointer::ReferenceCount
{
public:
	typedef core::Real                               Real;
	typedef core::Size                               Size;
	typedef core::Vector                             Vector;
	typedef numeric::geometry::BoundingBox< Vector > BoundingBox;
	typedef numeric::geometry::hashing::Bin6D        Bin6D;

	typedef utility::vector1< Size >                                              HitCounts;
	typedef boost::unordered_map< boost::uint64_t, HitCounts, numeric::geometry::hashing::bin_index_hasher >  HitHash;

public:
	MatchCounter();
	virtual ~MatchCounter();

	/// @brief Use the same bounding box as the HitHasher / OccupiedSpaceHash
	void
	set_bounding_box(
		BoundingBox const & bb
	);

	void
	set_n_geometric_constraints( Size ngeomcsts );

	/// @brief Give the same xyz bin witdh given to the HitHasher / OccupiedSpaceHash.
	void
	set_uniform_xyz_bin_width( Real bin_width );

	/// @brief Give the same euler-angle bin witdh given to the HitHasher / OccupiedSpaceHash.
	void
	set_uniform_euler_angle_bin_width( Real bin_width_degrees );


	/// @brief Give the same xyz bin witdhs given to the HitHasher / OccupiedSpaceHash.
	void
	set_xyz_bin_widths( Vector const & bin_widths );

	/// @brief Give the same euler-angle bin witdhs given to the HitHasher / OccupiedSpaceHash.
	void
	set_euler_bin_widths( Vector const & euler_bin_widths );

	/// @brief Call this after the bounding-box and the bin-widths have been set up.  Must
	/// be called before "add_hits" it called.  This initializes the SixDCoordinateBinner.
	void
	initialize();

	/// @brief Add hits from a list of hits for a particular geometric constraint.
	void add_hits( Size geomcst_id, std::list< Hit > const & hitlist );

	/// @brief Add hit from the input list of hits for a particular geometric constraint.
	void add_hits( Size geomcst_id, std::list< Hit const * >  const & hitptrlist );

	/// @brief Possibly slow method to predict the total number of matches given a set of
	/// hits and a particular grid resolution.  (The main function that this class provides).
	Size
	count_n_matches() const;


private:
	/*
	/// @brief Put the hits that have been added to "all_hits" into the hash_
	void hash_hits();

	Size
	find_next_bin(
	Real6 orig_point,
	Real6 offsets,
	utility::FixedSizeLexicographicalIterator< 6 > const & halfbin_lex,
	Bin6D & next_bin
	) const;

	bool
	within_reach(
	Bin6D const & query_halfbin,
	Bin6D const & nb_halfbin,
	Bin6D const & nbbin,
	utility::FixedSizeLexicographicalIterator< 6 > const & halfbin_lex
	) const;
	*/

private:

	BoundingBox bb_;

	Size n_geom_csts_;
	bool initialized_;

	utility::fixedsizearray1< Real, 3 > xyz_bin_widths_;
	utility::fixedsizearray1< Real, 3 > euler_bin_widths_;

	numeric::geometry::hashing::SixDCoordinateBinnerOP binner_;
	HitHash hash_;
};

/// @brief Increment the euler angles and then wrap them into their appropriate ranges
numeric::geometry::hashing::Real3
advance_euler_angles(
	numeric::geometry::hashing::Real3 const & orig_angles,
	numeric::geometry::hashing::Real3 const & offsets
);


typedef utility::vector1< Hit const * > match_lite;

struct match_lite_equals {

	bool
	operator() ( match_lite const & lhs, match_lite const & rhs ) const {
		if ( lhs.size() != rhs.size() ) return false;
		for ( core::Size ii = 1; ii <= lhs.size(); ++ii ) {
			if ( lhs[ ii ] != rhs[ ii ] ) {
				return false;
			}
		}
		return true;
	}

};

struct match_lite_hasher {
public:
	typedef core::Size Size;

public:
	Size
	operator() ( match_lite const & m ) const {
		/// Crazy hash function!
		Size hash = 1;
		for ( core::Size ii = 1; ii <= m.size(); ++ii ) {
			hash += ( hash * reinterpret_cast< Size > ( m[ ii ] ) * ii ) % 5527;
		}
		return hash % 7351;
	}
};

class MatchOutputTracker
{
public:

	typedef boost::unordered_map< match_lite, bool, match_lite_hasher, match_lite_equals > MatchHash;

public:
	MatchOutputTracker();

	void
	note_output_match( match_lite const & );

	bool
	match_has_been_output( match_lite const & m ) const;

private:
	MatchHash hash_;

};


}
}

#endif
