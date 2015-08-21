// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/downstream/SecondaryMatcherToUpstreamResidue.hh
/// @brief  Class declaration for secondary matcher that generates upstream-only hits
///         matching the geometry of one upstream residue with another upstream residue
///         generated in a previous round.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Florian Richter (flosopher@gmail.com)

#ifndef INCLUDED_protocols_match_downstream_SecondaryMatcherToUpstreamResidue_hh
#define INCLUDED_protocols_match_downstream_SecondaryMatcherToUpstreamResidue_hh

// Unit headers
#include <protocols/match/downstream/SecondaryMatcherToUpstreamResidue.fwd.hh>

// Package headers
#include <protocols/match/BumpGrid.fwd.hh>
#include <protocols/match/Matcher.hh>
#include <protocols/match/Hit.hh>
#include <protocols/match/downstream/DownstreamAlgorithm.hh>
#include <protocols/match/downstream/SecMatchResiduePairEvaluator.fwd.hh>
#include <protocols/match/upstream/UpstreamBuilder.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.fwd.hh>

// ObjexxFCL headers
#include <ObjexxFCL/FArray2D.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <list>
#include <map>
#include <string>

#include <core/id/AtomID.hh>
#include <utility/OrderedTuple.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace match {
namespace downstream {

/// @brief A class for an algorithm.  Given a conformation of the upstream partner,
/// the algorithm is responsible for producing a set of hits.
class SecondaryMatcherToUpstreamResidue : public DownstreamAlgorithm
{
public:
	typedef DownstreamAlgorithm                   parent;
	typedef utility::fixedsizearray1< Size, 2 >    Size2;
	typedef utility::OrderedTuple< Size2 >    Size2Tuple;
	typedef std::pair< SecMatchResiduePairEvaluatorCOP, Size > Evaluator_MCFI_ID_Pair;
	typedef std::list< Evaluator_MCFI_ID_Pair > EvaluatorSet;

public:
	SecondaryMatcherToUpstreamResidue( Size geom_cst_id );

	virtual ~SecondaryMatcherToUpstreamResidue();

	virtual
	DownstreamAlgorithmOP
	clone() const;

	/// @brief Main driver function for hit generation.  This DownstreamAlgorithm
	/// structures it's iteration over the hits from previous rounds as follows:
	/// for i = 1:n_build_positions
	///    recover_rotamer_coordinates_from_previous_round( hits_for_build_point_i );
	///    initialize rotcoord data for all downstream algorithms with the same geom_cst_id
	///    #omp parallel for /// All class access below this point is const and parallelizable
	///    for j = 1:n_build_positions
	///       /// call this function to start k loop: matcher.upstream_builder[ geom_cst_id() ]->build( j )
	///       for k = 1:n_rotamers_j
	///          /// call to start l loop: downstream_algorithm->build( j, k, rotamer_k ) )
	///          for l = 1:n_rotamers_from_build_point_i
	///             if ( respair_evaluator_->evaluate_residues( rotamer_k, rotamer_l )
	///                hit_list.append( Hit( j, k, ... ));
	///             return hit_list
	/// There are two important consequences to this hit-generation layout.
	/// 1. The coordinates for rotamer_k are computed n_build_position times.
	/// 2. Only a single build-position i has it's hit coordinates in memory at any point in time.
	/// This is a clear trade-off between performance and memory with a few caveats:
	/// A. It is very easy to bound the volume of space where build-point i's rotamers lie,
	/// so it should be easy to prune rotamer building, so rotamer k will be build many fewer than
	/// n_build_position times.
	/// B. By not trying to store all rotamers in memory at once, we do not impose any undue restrictions
	/// on the number of rotamers that can be enumerated.  This is especially important if we're
	/// using backbone flexibility to search a wider region of conformation space.
	virtual
	std::list< Hit >
	build_hits_at_all_positions(
		Matcher & matcher
	);


	/// @brief Prune hits away from the target_geomcst's hit list following a change to the
	/// hits for my geom_cst_id().  Pruning hits from the target_geomcst's hit list will
	/// trigger a round of peripheral-hitlist-change responses.
	virtual
	void
	respond_to_primary_hitlist_change( Matcher & matcher, Size round_just_completed );

	/// @brief Remove my hits if my target_geomcst's hit list has been shortened.  This
	/// will not trigger a round of peripheral-hitlist-change responses.
	virtual
	void
	respond_to_peripheral_hitlist_change( Matcher & matcher );


	/// @brief Iterate across the hits from a particular upstream build point i
	/// that were generated in a previous round, and see if the geometry of the
	/// input upstream_residue has "satisfactory interactions" with the
	/// hits from upstream-build-point i; if so, it appends a Hit to the hitlist
	/// returned at the end of the method.  (Also, see comments for the
	/// build_at_all_positions method.)
	virtual
	std::list< Hit >
	build(
		Size const scaffold_build_point_id,
		Size const upstream_conf_id,
		core::conformation::Residue const & upstream_residue
	) const;

	/// @brief returns true; this secondary matcher does not describe the location
	/// of the downstream partner
	virtual
	bool
	upstream_only() const;

	/// @brief This method returns 'false' since this matcher does not describe
	/// the coordinates of the downstream partner at all.
	virtual
	bool
	generates_primary_hits() const;


	/// @brief Prepare a map between upstream hits of the target-geomcst and
	/// a list of Hit const *'s of this geom_cst_id(). This map will be used
	/// in the function hits_to_include_with_partial_match.
	virtual
	void
	prepare_for_match_enumeration( Matcher const & matcher );

	/// @brief Return the set of hits to be iterated across
	virtual
	HitPtrListCOP
	hits_to_include_with_partial_match( match_dspos1 const & m ) const;

	virtual
	Size
	n_possible_hits_per_upstream_conformation() const;

	//void
	//set_match_restype( core::chemical::ResidueTypeCOP match_restype );

	void
	set_target_geomcst_id( Size target_geomcst_id );

	void
	add_target_restype( core::chemical::ResidueTypeCOP target_restype );

	void
	add_evaluator_for_target_restype(
		core::chemical::ResidueTypeCOP  target_restype,
		SecMatchResiduePairEvaluatorCOP evaluator,
		Size                            mcfi_id_for_evaluator
	);

	/// @brief Invoked by SecondaryMatchUpstreamResProcessor; avoids multiple inherritance,
	/// while letting the SecondaryMatcherToUpstreamResidue
	void
	process_hit(
		Hit const & hit,
		core::conformation::Residue const & upstream_conformation
	);

private:

	void
	prepare_for_hit_generation(
		Matcher & matcher
	);

	bool
	prepare_for_hit_generation_at_target_build_point(
		Matcher & matcher,
		upstream::ScaffoldBuildPoint const & target_build_point
	);

	/// @brief Allow another SecondaryMatcherToUpstreamResidue to set my
	/// TargetRotamerCoords object so that we can share this data.
	void
	set_target_rotamer_coords( TargetRotamerCoordsOP target_geomcst_coords );

	void count_rotamer(
		core::conformation::Residue const & upstream_conformation
	);

	void store_rotamer_coords(
		Hit const & hit,
		core::conformation::Residue const & upstream_conformation
	);

	void reorder_restypes(
		upstream::UpstreamBuilder const & builder
	);

private:

	utility::vector1< core::Size > smUR_pose_build_resids_;

	/// The id for the geom-cst that this SecondaryMatcher is dependent upon;
	/// the upstream geometries for residues generated by the target-geometric
	/// constraint are examined in the build() method.
	Size target_geomcst_id_;

	std::map< core::chemical::ResidueTypeCOP, Size >    target_restype_index_map_;
	utility::vector1< core::chemical::ResidueTypeCOP >  target_restypes_;
	utility::vector1< EvaluatorSet > respair_evaluators_;

	//  std::map < std::string, DownstreamBuilderCOP > dsbuilders_;

	bool count_rotamers_per_target_restype_;
	core::chemical::ResidueTypeCOP last_seen_restype_;
	Size last_seen_restype_index_;
	Size count_rotamer_for_lastseen_restype_;
	Matcher::HitListConstIterator target_hits_for_focused_build_point_begin_;
	Matcher::HitListConstIterator target_hits_for_focused_build_point_end_;
	Matcher::HitListConstIterator target_hits_end_;

	utility::vector1< Size > n_rotamers_per_target_restype_;
	TargetRotamerCoordsOP target_geomcst_coords_;

	std::map< Size2Tuple, HitPtrListOP > my_hits_for_target_hit_map_;

	Size n_viable_target_hits_since_last_pruning_;
};

class TargetRotamerCoords : public utility::pointer::ReferenceCount
{
public:
	typedef core::Size   Size;
	typedef core::Vector Vector;

public:

	TargetRotamerCoords();
	virtual ~TargetRotamerCoords();

	void set_num_restypes( Size n_restypes );
	void set_restype( Size restype_index, core::chemical::ResidueTypeCOP restype );
	void set_required_atoms( Size restype_index, utility::vector1< bool > const & atom_required );

	void set_num_target_rotamers(
		utility::vector1< Size > const & n_rotamers_per_target_restype
	);

	void set_num_target_rotamers(
		Size target_restype_id,
		Size n_rotamers
	);

	void set_coordinates_for_rotamer(
		Size restype_index,
		Size rotamer_index,
		Hit const & hit,
		core::conformation::Residue const & rescoords
	);

	Size
	n_restypes() const {
		return target_restypes_.size();
	}

	core::chemical::ResidueTypeCOP
	restype( Size restype_index ) const {
		return target_restypes_[ restype_index ];
	}

	Size
	n_rots_total() const {
		return n_rots_total_;
	}

	Size
	n_rotamers_for_restype( Size restype_id ) const {
		return coords_[ restype_id ].size2();
	}

	Size
	n_atoms_for_restype( Size restype_id ) const {
		return atom_ids_for_coordinates_[ restype_id ].size();
	}

	Vector const &
	coord(
		Size restype_index,
		Size rotamer_index,
		Size which_atom
	) const {
		return coords_[ restype_index ]( which_atom, rotamer_index );
	}

	Size
	restype_atomno(
		Size restype_index,
		Size which_atom
	) const {
		return atom_ids_for_coordinates_[ restype_index ][ which_atom ];
	}

	Hit const &
	hit (
		Size restype_index,
		Size rotamer_index
	) const {
		return hit_data_[ restype_index ][ rotamer_index ];
	}

	void
	set_clash_checking( Size rotamer_index) {
		build_coords_for_clash_checking_[ rotamer_index ] = true;
	}

	bool
	get_clash_checking( Size rotamer_index) {
		return build_coords_for_clash_checking_[ rotamer_index ];
	}

	void
	set_coords_for_clash_check(
		Size rotamer_index,
		utility::vector1< Vector > & coords
	) {
		coords_for_clash_checking_[ rotamer_index ].resize( coords.size() );
		coords_for_clash_checking_[ rotamer_index ] = coords;
	}

	utility::vector1< Vector >
	get_coords_for_clash_check(
		Size rotamer_index
	) {
		return coords_for_clash_checking_[ rotamer_index ];
	}

	void
	set_clash_check_types(
		Size n_rotamers
	){
		build_coords_for_clash_checking_.resize( n_rotamers );
		coords_for_clash_checking_.resize( n_rotamers );
	}

	void
	set_ds_atom_ids_needed(
		utility::vector1 < core::id::AtomID  > atom_ids
	){
		ds_atom_ids_needed_ = atom_ids;
	}

	utility::vector1 < core::id::AtomID  >
	get_ds_atom_ids_needed(
	){
		return ds_atom_ids_needed_;
	}


private:
	utility::vector1< core::chemical::ResidueTypeCOP > target_restypes_;
	utility::vector1< utility::vector1< Size > > atom_ids_for_coordinates_;
	utility::vector1< ObjexxFCL::FArray2D< Vector > > coords_;
	Size n_rots_total_;
	// In parallel to the coordinates for each of the target hits, keep the
	// hit data that we need for later generating hits for the match
	// first = scaffold_build_id, second = upstream_conf_id;
	utility::vector1< utility::vector1< Hit > > hit_data_;

	utility::vector1< bool > build_coords_for_clash_checking_;
	utility::vector1< core::id::AtomID  > ds_atom_ids_needed_;
	utility::vector1< utility::vector1< Vector > > coords_for_clash_checking_;
};

/// @brief A simple class to respond to the UpstreamBuilder's
/// process_hit method and pass on the coordinates to its "owning"
/// SecondaryMatcherToUpstreamResidue object.
class SecondaryMatchUpstreamResProcessor : public upstream::UpstreamResidueProcessor
{
public:
	typedef upstream::UpstreamResidueProcessor parent;

public:

	SecondaryMatchUpstreamResProcessor(
		SecondaryMatcherToUpstreamResidue & sec_matcher
	);

	virtual
	void
	process_hit(
		Hit const & hit,
		core::conformation::Residue const & upstream_conformation
	);

private:

	SecondaryMatcherToUpstreamResidue & sec_matcher_;

};


/// @brief A simple struct to use in list.sort() to ensure that the
/// hits returned by a secondary matcher which has possibly generated upstream
/// hits out-of-order, will return an ordered-hit-list in its
/// build_hits_at_all_positions() method.
///
/// @details This struct compares the upstream portion of the hits it's returning
/// ensuring that the rotamer indices (the upstream_conf_ids()) are in ascending
/// order for each scaffold build point.
struct us_secmatch_hit_compare
{
	bool operator () ( Hit const & lhs, Hit const & rhs ) const {
		if ( lhs.scaffold_build_id() == rhs.scaffold_build_id() ) {
			return lhs.upstream_conf_id() < rhs.upstream_conf_id();
		} else {
			return lhs.scaffold_build_id() < rhs.scaffold_build_id();
		}
	}
};


}
}
}

#endif
