// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/sewing/hashing/Hasher.hh
///
/// @brief A geometric hashing class used by the SEWING protocol. The hasher functions by hashing the coordinates
/// of a given Model after transformation into a local coordinate from defined by the N,CA,C atoms of each residues.
/// Once a hash table is generate, new models (or models that have already been inserted into the table) can be "scored"
/// against the table by following the same transformation procedure. Any two residues from different models that hash many
/// of their models to the same quarter angstrom voxels indicate a good alignment between two models.
///
/// @author Tim Jacobs

#ifndef INCLUDED_protocols_sewing_hashing_Hasher_hh
#define INCLUDED_protocols_sewing_hashing_Hasher_hh

//Unit headers
#include <protocols/sewing/hashing/Hasher.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

//Package headers
#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <protocols/sewing/conformation/Model.hh>

//Utility headers
#include <core/types.hh>
#include <utility/vector1.fwd.hh>
#include <utility/fixedsizearray1.hh>
#include <numeric/xyzVector.hh>
#include <numeric/HomogeneousTransform.fwd.hh>
#include <boost/unordered_map.hpp>
#include <boost/functional/hash.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

//C++ headers
#include <map>
#include <set>


#include <basic/options/option.hh>
#include <basic/options/keys/sewing.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>


namespace protocols {
namespace sewing  {


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
	int model_id;//the identifier for the entire self-contained unit
	core::Size basis_resnum;//the residue used to make the basis set(N,CA,C).

	core::Size segment_id;//The segment from which the atom for this key was found
	core::Size resnum;//the residues from which the atom for this key was found
	core::Size atomno;//the atomno for the atom

	friend
	bool
	operator< ( HashValue const & a, HashValue const & b ) {
		if ( a.model_id < b.model_id ) return true;
		if ( a.model_id > b.model_id ) return false;

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

typedef boost::unordered_map<HashKey,utility::vector1< HashValue >,coord_hash,coord_equal_to> HashMap;

/////////////////////////// End Main HashMap for model coordinates ///////////////////////////////

/////////////////////////// HashMap for score results ///////////////////////////////

struct Basis {
	Basis():
		model_id(0),
		resnum(0)
	{}

	Basis(int model_id_temp, core::Size resnum_temp){
		model_id = model_id_temp;
		resnum = resnum_temp;
	}

	friend
	bool
	operator< (Basis const & a, Basis const & b) {
		//if ( a.model_id < b.model_id ) return true;
		//if ( a.model_id > b.model_id ) return false;

		//if ( a.resnum < b.resnum ) return true;
		//if ( a.resnum > b.resnum ) return false;

		//return false;

		return boost::tie(a.model_id, a.resnum) < boost::tie(b.model_id, b.resnum);
	}

	int model_id;
	core::Size resnum;
};

typedef std::pair< Basis, Basis > BasisPair;
typedef std::pair< core::Size, core::Size> SegmentPair;
typedef std::map< core::id::AtomID, core::id::AtomID > AtomMap;

struct basis_pair_hash : std::unary_function<BasisPair, core::Size> {
	core::Size operator()(BasisPair const & basis_pair) const {
		size_t seed = 0;
		boost::hash_combine(seed, basis_pair.first.model_id);
		boost::hash_combine(seed, basis_pair.first.resnum);
		boost::hash_combine(seed, basis_pair.second.model_id);
		boost::hash_combine(seed, basis_pair.second.resnum);
		return seed;
	}
};

struct basis_pair_equal_to : std::binary_function<BasisPair, BasisPair, bool> {
	bool operator()(BasisPair const & bp1, BasisPair const & bp2) const {
		return
			bp1.first.model_id == bp2.first.model_id &&
			bp1.first.resnum == bp2.first.resnum &&
			bp1.second.model_id == bp2.second.model_id &&
			bp1.second.resnum == bp2.second.resnum;
	}
};

struct HashResult {

	HashResult():
		clash_count(0)
	{}

	std::map< SegmentPair, AtomMap > segment_matches;
	std::map< SegmentPair, core::Size > segment_match_counts;
	core::Size clash_count;
};

typedef std::map< BasisPair, HashResult > ScoreResults;
//typedef boost::unordered_map< BasisPair, HashResult, basis_pair_hash, basis_pair_equal_to > ScoreResults;
typedef std::pair< BasisPair, HashResult > ScoreResult;

/////////////////////////// End HashMap for score results ///////////////////////////////



class Hasher : public utility::pointer::ReferenceCount
{

public:

	///@brief default constructor
	Hasher();

	ScoreResult
	score_one(
		Model const & m1,
		SewResidue const & m1_basis,
		Model const & m2,
		SewResidue const & m2_basis,
		core::Size box_length
	);

	///@brief Insert this model into the hash table
	void
	insert(
		Model const & model
	);

	// The 2nd score function with 8 arguments in Hasher.cc
	ScoreResults
	score(
		Model const & model,
		core::Size num_segment_matches,
		core::Size min_segment_score,
		core::Size max_clash_score,
		std::set<core::Size> const & score_segments,
		bool store_atoms,
		core::Size box_length
	) const;


	// The 1st score function with 7 arguments in Hasher.cc
	///@brief Score the given model against the models in the hash table
	ScoreResults
	score(
		Model const & model,
		core::Size num_segment_matches,
		core::Size min_segment_score,
		core::Size max_clash_score,
		bool store_atoms,
		core::Size box_length
	) const;


	///@brief
	void
	score_basis( // with "regular" 27 neighborhood lookup boxes
		ScoreResults & alignment_scores,
		Model const & transformed_model,
		SewResidue const & basis_residue,
		bool store_atoms
	) const;

	// very similar as in score_basis
	///@brief
	void
	score_basis_125( // with "larger" 125 neighborhood lookup boxes
		ScoreResults & alignment_scores,
		Model const & transformed_model,
		SewResidue const & basis_residue,
		bool store_atoms
	) const;

	///@brief Trim the ScoreResults to remove weak hits
	void
	trim_scores(
		ScoreResults & scores,
		core::Size num_segment_matches,
		core::Size min_segment_score,
		core::Size max_clash_score
	) const;

	///@brief keep only best scoring alignment between two models
	ScoreResults
	remove_duplicates(
		ScoreResults const & scores
	) const;

	///@brief remove edges between segments that both have 'next' or 'previous' segments
	void
	remove_connection_inconsistencies(
		std::map< int, Model > const & models,
		ScoreResults & scores
	) const;

	///@brief const accessor to the underlying hash map
	HashMap const &
	hash_map() const;

	///@brief retrive hits from the bin corresponding to the key (3D voxel) and all neighboring
	///quarter-angstrom bins
	//utility::vector1<HashValue>
	void
	neighborhood_lookup(
		HashKey const & key,
		utility::fixedsizearray1<HashMap::const_iterator, 27> & hit_its
	) const;

	// new with box_length
	///@brief retrive hits from the bin corresponding to the key (3D voxel) and all neighboring
	///quarter-angstrom bins
	//utility::vector1<HashValue>
	void
	neighborhood_lookup_125(
		HashKey const & key,
		utility::fixedsizearray1<HashMap::const_iterator, 125> & hit_its
	) const;

	///@brief Transform all features to the local coordinate frame of the basis set
	Model
	transform_model(
		Model const & model,
		SewResidue const & basis_residue
	) const;

	///@brief Hash the transformed residues into the HashMap
	void
	hash_model(
		Model const & transformed_model,
		SewResidue const & basis_residue
	);

	///@brief Create the HashKey from a SewAtom
	HashKey
	generate_key(
		SewAtom const & atom
	) const;

	///@brief Serialize hash table to disk
	void
	write_to_disk(std::string filename) const;

	///@brief Populate the hash table from one on disk.
	void
	read_from_disk(std::string filename);

private:

	HashMap hash_map_;

};


} //sewing namespace
} //protocols namespace

#endif
