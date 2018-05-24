// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/hashing/hasher_data.hh
/// @brief Data storage for the Hasher
/// @author Minnie Langlois (minnie@email.unc.edu)

#ifndef INCLUDED_protocols_sewing_hashing_hasher_data_hh
#define INCLUDED_protocols_sewing_hashing_hasher_data_hh

//Sewing headers
#include <protocols/sewing/data_storage/SmartSegment.hh>
//External headers
#include <boost/unordered_map.hpp>
#include <boost/functional/hash.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <map>
#include <set>
// Utility headers
#include <utility/fixedsizearray1.hh>
#include <utility/vector1.hh>

#include <core/id/AtomID.hh> //forward header isn't working for some reason?
#include <core/types.hh>
namespace protocols {
namespace sewing {
namespace hashing {

/////////////////////////// Main HashMap for model coordinates ///////////////////////////////

//The HashKey is the 3 indices that describe the location the feature in discretized space, and the
//atomno of the point being described. The hasher
//is responsible for creating a (hopefully) uniform distribution of keys in the table.
typedef utility::fixedsizearray1<int, 3> HashKey;

struct coord_hash : std::unary_function<HashKey, core::Size> {
	core::Size operator()(HashKey const & key) const {
		size_t seed = 0;
		boost::hash_combine(seed, key[1]);
		boost::hash_combine(seed, key[2]);
		boost::hash_combine(seed, key[3]);
		return seed;
	}
};

struct coord_equal_to : std::binary_function<HashKey, HashKey, bool> {
	typedef utility::fixedsizearray1const_iterator<int, 3> const_iterator;
	bool operator()(HashKey const & key1, HashKey const & key2) const {
		return key1[1] == key2[1] && key1[2] == key2[2] && key1[3] == key2[3];
	}
};

//The HashValue is the BasisSet (3 points), and the Model (Pose) used to define the coordinate
//frame
struct HashValue {
	core::Size basis_resnum;//the residue used to make the basis set(N,CA,C).
	core::Size basis_segment_id; //The segment used to make the basis set
	core::Size segment_id;//The segment from which the atom for this key was found
	core::Size resnum;//the residues from which the atom for this key was found
	core::Size atomno;//the atomno for the atom

	friend
	bool
	operator< ( HashValue const & a, HashValue const & b ) {
		if ( a.basis_resnum < b.basis_resnum ) return true;
		if ( a.basis_resnum > b.basis_resnum ) return false;

		if ( a.segment_id < b.segment_id ) return true;
		if ( a.segment_id > b.segment_id ) return false;

		if ( a.resnum < b.resnum ) return true;
		if ( a.resnum > b.resnum ) return false;

		if ( a.atomno < b.atomno ) return true;
		if ( a.atomno > b.atomno ) return false;

		return false;
	}
};

typedef boost::unordered_map<HashKey,HashValue,coord_hash,coord_equal_to> HashMap;

/////////////////////////// End Main HashMap for model coordinates ///////////////////////////////

typedef std::map< core::id::AtomID, core::id::AtomID > AtomMap;


/////////////////////////// End HashMap for score results ///////////////////////////////

////////////////////////// Begin HasherSettings ////////////////////////////////////////
struct HasherSettings {

	HasherSettings():
		max_clash_score(0),
		min_hash_score(0),
		boxes_per_dimension(0),
		hash_between_termini(false)
	{}

	HasherSettings(
		core::Size max_clash_score_temp,
		core::Size min_hash_score_temp,
		core::Size boxes_per_dimension_temp,
		bool hash_between_termini_temp
	) {
		max_clash_score = max_clash_score_temp;
		min_hash_score = min_hash_score_temp;
		boxes_per_dimension = boxes_per_dimension_temp;
		hash_between_termini = hash_between_termini_temp;
	}
	core::Size max_clash_score;
	core::Size min_hash_score;
	core::Size boxes_per_dimension;
	bool hash_between_termini;

};

//SegmentVector
class SegmentVector: public utility::vector1< data_storage::SmartSegmentOP > {
	//just a vector1 of SmartSegmentOPs with a fancy destructor
public:
	virtual ~SegmentVector() {
		for ( data_storage::SmartSegmentOP cur_seg: *this ) {
			cur_seg->set_n_terminal_neighbor(nullptr);
			cur_seg->set_c_terminal_neighbor(nullptr);
		}
	}
	core::Size get_version() const{ return version_; }
	void set_version( core::Size ver ){ version_ = ver; }

private:
	core::Size version_;
};

typedef utility::pointer::shared_ptr< SegmentVector > SegmentVectorOP;
typedef utility::pointer::shared_ptr< SegmentVector const > SegmentVectorCOP;

//EdgeMap
typedef std::map< core::Size, std::set< core::Size > > EdgeMap;
typedef utility::pointer::shared_ptr< EdgeMap > EdgeMapOP;

//BasisMap
typedef std::map< core::Size, std::map < core::Size, utility::vector1< std::pair< core::Size, core::Size > > > > BasisMap;
typedef utility::pointer::shared_ptr< BasisMap > BasisMapOP;
//Map of segment id to a map of second segment id to vector of possible basis pairs between them (resnums only)

///////////////// Begin AlignmentSettings ///////////////////
struct AlignmentSettings {
	AlignmentSettings():
		recursive_depth_(2)
	{
		match_segments_.clear();
		segments_.clear();
		//pose_segment_starts_.clear();
		//pose_segment_ends_.clear();
	}
	AlignmentSettings(utility::vector1< core::Size > match_segments, utility::vector1< std::pair< core::Size, core::Size > > segments, core::Size recursive_depth):
		match_segments_( match_segments ),// which segments to hash - commandline option
		segments_( segments ), // residues that define segment (start & end) - commandline option
		recursive_depth_( recursive_depth )
	{}
	utility::vector1< core::Size > match_segments_;
	utility::vector1< std::pair< core::Size, core::Size > > segments_;
	core::Size recursive_depth_;

	//How many edges in do you want to go for each input segment? E.g. recursive_depth_ = 2 generates edges from the model and
	//edges from segments that have edges to the model.
};

struct IdealContact{
	core::Size ligand_atom = 1;
	//Required contacts indicate 1) the atom number within the ligand and 2) the number of contacts required for that atom
	//core::Size min_contacts;
	core::Size max_contacts = 4;
	core::Real preferred_distance = 2.1;
	core::Real preferred_contact_angle = 109.5; //Angle between contacts
	core::Real preferred_dihedral_1 = 30; // Dihedral about ligand atom: contact_base - contact - ligand_atom - other_contact
	core::Real preferred_dihedral_2 = 30; // Dihedral about ligand atom: contact - ligand_atom - other_contact - other_base
};


} //hashing
} //protocols
} //sewing


#endif //protocols/sewing/hashing/hasher_data_hh

