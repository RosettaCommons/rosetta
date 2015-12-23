// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loophash/LoopHashMap.hh
/// @brief
/// @author Mike Tyka
/// @author Ken Jung

#ifndef INCLUDED_protocols_loophash_LoopHashMap_hh
#define INCLUDED_protocols_loophash_LoopHashMap_hh

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
//#include <protocols/match/Hit.fwd.hh>
//#include <protocols/match/SixDHasher.hh>
#include <boost/unordered_map.hpp>
#include <utility/pointer/ReferenceCount.hh>
#include <protocols/loophash/BackboneDB.hh>

#include <numeric/geometry/hashing/SixDHasher.hh>

#include <string>
#include <vector>
#include <map>

#include <utility/exit.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace loophash {

/// @brief This takes a pose and two residue positions and determines the rigid body transform of the Leap described by those two residues.
///        Returns true is successful or false if something went haywire and one should just ignore this loop (this can happen at the ends)
bool get_rt_over_leap( const core::pose::Pose& orig_pose, core::Size ir, core::Size jr, numeric::geometry::hashing::Real6 &rt_6 );

/// @brief This takes a pose and two residue positions and determines the rigid body transform of the Leap described by those two residues
///        Returns true is successful or false if something went haywire and one should just ignore this loop (this can happen at the ends)
///       THe difference between this and the get_rt_over_leap function is that this version doesnt make a copy of the pose which makes it faster.
///     However this means that the pose passed cannot be a const pose, even though the function restores the fold tree afterwards..
bool get_rt_over_leap_fast( core::pose::Pose& pose, core::Size ir, core::Size jr, numeric::geometry::hashing::Real6 &rt_6 );

/// @brief This takes a pose and two residue positions and determines the rigid body transform of the Leap described by those two residues.
// sheffler
bool
get_rt_over_leap_without_foldtree_bs(
	core::pose::Pose const & pose,
	core::Size ir,
	core::Size jr,
	numeric::geometry::hashing::Real6 &rt_6
);


/// @brief The LeapIndex stores information about a particular Leap. It hold the oroiginal high precision rigid body transform
///  and an Index to a Backbone Database (BackboneDB) that has the actual phi psi angles. THe storage of the precise RT takes a
/// lot of space and may be deprecated at some point, since once it is hashed, it is rarely needed and can be recomputed if it is.
/// Note that the length of the loop is not stored either, this is again done for saving memory, as huge lists of Leaps are
/// typically created all with the same length. THe length is stored and handled by the owner of LeapIndex list.
/// The LeapIndex does not store the actual backbone coordinates of the Leap. It merely contains an index (the BackboneIndex)
/// which refers to a serial store of backbone triples (phi,psi, omega) which are stored somewhere else in a BackboneDB.
/// THis is improtant to save memory storage since multiple Leaps cna share the same backbone triple and redundant storage
/// would be hugely wasteful.
struct LeapIndex{
	core::Size index;
	core::Size offset;
	boost::uint64_t key;
};
struct LegacyLeapIndex{
	float vecx;
	float vecy;
	float vecz;
	float rotx;
	float roty;
	float rotz;
	unsigned int ba;
};

/// @brief the loop hash map stores LeapIndexes and a hashmap to access those LeapIndexes quickly by their 6D coordinates.


class LoopHashMap{
public:
	/// @brief Constructor - must give loop_size
	LoopHashMap( core::Size loop_size = 10);

	LoopHashMap(LoopHashMap const & other);

	LoopHashMap & operator=(LoopHashMap const & other);

private:
	/// @brief Setup this class give a loop_size
	void setup( core::Size loop_size);

	/// @brief Add a given piece of data to the hash
	void hash_index( numeric::geometry::hashing::Real6 transform, core::Size data );
public:

	/// @brief Add a leap/loop to the HashMap
	void add_leap( const LeapIndex &leap_index, numeric::geometry::hashing::Real6 & transform );

	/// @brief Add a leap/loop to the HashMap with a key, skipping hashing
	void add_leap( const LeapIndex &leap_index, boost::uint64_t key );

	/// @brief Add an legacy leap/loop to the HashMap
	void add_legacyleap( const LegacyLeapIndex &legacyleap_index );

	/// @brief Obtain an index to a given peptide saved
	inline const LeapIndex & get_peptide( core::Size index ){
		runtime_assert( index < loopdb_.size() );
		return loopdb_[ index ];
	}

	inline core::Size n_loops() const {
		return loopdb_.size();
	}

	/// @brief Return a vector of loops with equal keys  given a key
	//  And any keys within a radius around the original key
	void radial_lookup_withkey( boost::uint64_t key, core::Size radius, std::vector < core::Size > &result );

	/// @brief Return a vector of loops with equal keys  given a key
	//  Identical to radial_lookup_withkey(radius=0), but faster
	void lookup_withkey( boost::uint64_t key, std::vector < core::Size > &result );

	/// @brief Append to a bucket of vectors in the appropriate bin, lookup by transform
	void lookup(  numeric::geometry::hashing::Real6 transform, std::vector < core::Size > &result );

	/// @brief Append to a bucket of vectors in the appropriate bin, radial lookup by transform
	void radial_lookup( core::Size radius,  numeric::geometry::hashing::Real6 transform, std::vector < core::Size > &result );

	/// @brief count hits in the appropriate bin, radial lookup by transform
	core::Size radial_count( core::Size radius, numeric::geometry::hashing::Real6 center ) const;

	/// @brief Append to a bucket of vectors in the appropriate bin, lookup by bin index
	/// Using core::Size instead of boost::uinst64_t
	void lookup( core::Size index, std::vector < core::Size > &result );

	/// @brief Returns begin() and end() of backbone_index_map_
	void  bbdb_range( std::pair< BackboneIndexMap::iterator, BackboneIndexMap::iterator > & range );

	/// @brief Returns a hashmap key given a member of a bucket
	/// Don't think boost implements this, have to manually look it up
	boost::uint64_t return_key( core::Size bb_index );

	/// @brief Query the loopsize of this LoopHashMap
	inline core::Size get_loop_size() const { return loop_size_; }

	/// @brief Reads legacy binary dbs
	void read_legacydb( std::string filename );

	/// @brief Basic IO functionality - allows reading and writing text states to/from disk
	void write_db( std::string filename );

	/// @brief Basic IO functionality - allows reading and writing text states to/from disk
	void read_db( std::string filename, std::pair< core::Size, core::Size > loopdb_range,
		std::map< core::Size, bool > & homolog_index );
	// hack to set loopdb_range default value
	inline void read_db( std::string filename ) {
		std::pair< core::Size, core::Size > range( 0, 0 );
		std::map< core::Size, bool > homolog_index;
		read_db( filename, range, homolog_index );
	}

	/// @brief Return the memory usage of this class
	void mem_foot_print();

	/// @brief Sorts the loopdb_ by leap_index.index
	void sort();

private:  // Private data

	/// @brief  A class that will take a 6D rigid body transform and turn it into a serial hashbin
	/// number for hashing. THis is the actual hash so to speak.
	numeric::geometry::hashing::SixDCoordinateBinnerOP  hash_;

	/// @brief The actual Boost-based hashmap
	BackboneIndexMap                          backbone_index_map_;

	/// @brief List of LeadIndexes
	std::vector< LeapIndex>                   loopdb_;

	/// @brief The length of the the loops in number of residues
	core::Size                                loop_size_;

	/// @brief A functor for sort()
	struct by_index {
		bool operator()( LeapIndex const &a, LeapIndex const &b ) const {
			return a.index < b.index;
		}
	};


};


} // namespace loops
} // namespace protocols


#endif


