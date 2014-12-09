// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/downstream/SecondaryMatcherToDownstreamResidue.hh
/// @brief  Class declaration for secondary matcher that generates upstream-only hits
///         matching the geometry of one upstream residue with another upstream residue
///         generated in a previous round.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Florian Richter (flosopher@gmail.com)

#ifndef INCLUDED_protocols_match_downstream_SecondaryMatcherToDownstreamResidue_hh
#define INCLUDED_protocols_match_downstream_SecondaryMatcherToDownstreamResidue_hh

// Unit headers
#include <protocols/match/downstream/SecondaryMatcherToDownstreamResidue.fwd.hh>

// Package headers
#include <protocols/match/Matcher.hh>
// AUTO-REMOVED #include <protocols/match/Hit.hh>
#include <protocols/match/downstream/DownstreamAlgorithm.hh>
#include <protocols/match/downstream/SecMatchResiduePairEvaluator.fwd.hh>
// AUTO-REMOVED #include <protocols/match/upstream/UpstreamBuilder.hh>
#include <protocols/match/downstream/SecondaryMatcherToUpstreamResidue.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
// AUTO-REMOVED #include <core/id/AtomID.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
// C++ headers
#include <list>

#include <core/id/AtomID.hh>


namespace protocols {
namespace match {
namespace downstream {

/// @brief A class for an algorithm.  Given a conformation of the downstream partner,
/// the algorithm is responsible for producing a set of hits.
class SecondaryMatcherToDownstreamResidue : public DownstreamAlgorithm
{
public:
	typedef DownstreamAlgorithm                   parent;
	typedef std::pair< SecMatchResiduePairEvaluatorCOP, Size > Evaluator_MCFI_ID_Pair;
	typedef std::list< Evaluator_MCFI_ID_Pair > EvaluatorSet;

public:
	SecondaryMatcherToDownstreamResidue( core::pose::PoseCOP upstream_pose, Size geom_cst_id );

	virtual ~SecondaryMatcherToDownstreamResidue();

	virtual
	DownstreamAlgorithmOP
	clone() const;

	/// @brief Main driver function for hit generation.  This DownstreamAlgorithm
	/// structures it's iteration over the hits from previous rounds as follows:
	/// for i = 1:this->geom_cst_id() - 1
	///    if ( ! matcher->representative_downstream_algorithm( i )->generates_primary_hits() ) continue;
	///    for j = 1:n_build_points_for_geomcst( i )
	///       recover_downstream_coordinates_from_previous_round( hit_subset_j );
	///       initialize TaretRotamerCoords data for all downstream algorithms with the same geom_cst_id
	///       #omp parallel for /// All class access below this point is const and parallelizable
	///       for k = 1:n_build_positions
	///          /// call this function to start l loop: matcher.upstream_builder[ geom_cst_id() ]->build( k )
	///             for l = 1:n_rotamers_k
	///             /// call to start m loop: downstream_algorithm->build( k, l, rotamer_l ) )
	///             for m = 1:n_hits_in_block_j
	///                if ( respair_evaluator_->evaluate_residues( rotamer_l, rotamer_m )
	///                   hit_list.append( Hit( k, l, i, 1, hit[ m ].second() ));
	///                return hit_list
	/// There are two important consequences to this hit-generation layout.
	/// 1. The coordinates for rotamer_l are computed sum( i, n_build_points_for_geomcst( i )) times.
	/// 2. The number of downstream target coordinates that live in memory at the same time is bound by some constant (10K).
	/// This is a clear trade-off between performance and memory
	/// NOTE: the most time consuming portion will likely be the m loop, and not the repeated construction
	/// of coordinates in the j loop.  Reguardless of how many times we rebuild coordinates for rotamer j, the expense
	/// will primarily lie in the same place: the call to evaluate_residues( rotamer_l, rotamer_m ).
	/// NOTE: if there are ways to iterate across the j loop differently, it likely possible to prune m/l combinations early
	/// and thereby improve running time.
	virtual
	std::list< Hit >
	build_hits_at_all_positions(
		Matcher & matcher
	);


	/// @brief mimic the classic matcher's reset of the Occupied space hash.
	virtual
	void
	respond_to_primary_hitlist_change( Matcher & matcher, Size round_just_completed );

	/// @brief Remove my hits if they fall into a volume of the occupied space hash
	/// that is no longer occupied.
	virtual
	void
	respond_to_peripheral_hitlist_change( Matcher & matcher );


	/// @brief Iterate across the conformations of the downstream residue coming from hits
	/// generated in previous rounds, and add hits for each upstream residue that
	/// (Also, see comments for the
	/// build_at_all_positions method.)
	virtual
	std::list< Hit >
	build(
		Size const scaffold_build_point_id,
		Size const upstream_conf_id,
		core::conformation::Residue const & upstream_residue
	) const;

	/// @brief returns false; this secondary matcher describes the location
	/// of the downstream partner even though it does not generate that location
	/// itself.  Matches may be found by hashing the 6D coordinate of the
	/// downstream partner.
	virtual
	bool
	upstream_only() const;

	/// @brief This method returns 'false' since this matcher does not describe
	/// the coordinates of the downstream partner at all.
	virtual
	bool
	generates_primary_hits() const;

	HitPtrListCOP
	hits_to_include_with_partial_match( match_dspos1 const & ) const;

	virtual
	Size
	n_possible_hits_per_upstream_conformation() const;

	void
	set_downstream_restype( core::chemical::ResidueTypeCOP downstream_restype );

	void
	set_focused_geomcst_id( Size focused_geomcst_id );

	void
	add_evaluator(
		SecMatchResiduePairEvaluatorCOP evaluator,
		Size mcfi_id
	);

  void
  set_catalytic_atoms(
    utility::vector1< core::Size > catalytic_atoms
  ){
    catalytic_atoms_ = catalytic_atoms;
  }

//	void
//	set_dsbuilders(
//  	std::map<std::string, DownstreamBuilderCOP> dsbuilders
//	);

	//assign minimum separation distance to downstream-upstrean pair
	//min_sep_d2_from_upstream_atoms_[ downstream ][ upstream ].second distance
	//min_sep_d2_from_upstream_atoms_[ downstream ][ upstream ].first
//	void
//	initialize_upstream_residue(
//  	core::chemical::ResidueTypeCOP  us_res /*upstream residue*/
//	);

private:

	void
	prepare_for_hit_generation(
		Matcher & matcher
	);

	void
	prepare_for_hit_generation_for_geomcst(
		Matcher & matcher,
		Size target_geomcst_id
	);

	bool
	prepare_for_hit_generation_at_target_build_point(
		Matcher & matcher,
		Size target_geomcst_id,
		upstream::ScaffoldBuildPoint const & target_build_point
	);

	/// @brief Allow another SecondaryMatcherToDownstreamResidue to set my
	/// TargetRotamerCoords object so that we can share this data.
	void
	set_target_rotamer_coords( TargetRotamerCoordsOP target_geomcst_coords );

	void
	set_ds_coords_needed(
		utility::vector1< core::id::AtomID > ds_coords
	){
		downstream_atom_coordinates_needed_ = ds_coords;
	}

	bool
	ds_atom_present(
  	Size index
	) const;

  utility::vector1< utility::vector1< utility::vector1< std::pair< core::Size, core::Real > > > >
  get_min_sep_d2_from_upstream_atoms() const;

private:

	core::chemical::ResidueTypeCOP downstream_restype_;
	core::pose::PoseCOP upstream_pose_;
	EvaluatorSet respair_evaluators_;
//  std::map < std::string, DownstreamBuilderCOP > dsbuilders_;

	Size focused_geomcst_id_;
	Matcher::HitListConstIterator hits_for_focused_geomcst_and_build_point_begin_;
	Matcher::HitListConstIterator hits_for_focused_geomcst_and_build_point_end_;
	Matcher::HitListConstIterator hits_for_focused_geomcst_end_;

	mutable TargetRotamerCoordsOP target_downstream_coords_;
	utility::vector1< core::id::AtomID > downstream_atom_coordinates_needed_;
  //catalytic_atoms_[2] and catalytic_atoms_[1] upstream two atoms
  //(built from distance and angle) from the constraint
  //catalytic_atoms_[3] and catalytic_atoms_[4] downstream two atoms
  //(built from distance and angle) from the constraint
	utility::vector1< core::Size > catalytic_atoms_;

	//utility::vector1< core::id::AtomID > downstream_atoms_for_clash_checking_ ;
	Size occspace_rev_id_at_last_update_;
};

}
}
}

#endif
