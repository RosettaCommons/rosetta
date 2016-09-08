// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/downstream/SecondaryMatcherToDownstreamResidue.cc
/// @brief  Class implementation for secondary matcher that generates upstream-only hits
///         matching the geometry of one upstream residue with another upstream residue
///         generated in a previous round.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Florian Richter (flosopher@gmail.com)


// Unit headers
#include <protocols/match/downstream/SecondaryMatcherToDownstreamResidue.hh>

// Package headers
#include <protocols/match/Matcher.hh>
#include <protocols/match/Hit.hh>
#include <protocols/match/OccupiedSpaceHash.hh>
#include <protocols/match/downstream/DownstreamAlgorithm.hh>
#include <protocols/match/downstream/DownstreamBuilder.hh>
#include <protocols/toolbox/match_enzdes_util/LigandConformer.hh>
#include <protocols/match/downstream/SecMatchResiduePairEvaluator.hh>
#include <protocols/match/downstream/SecondaryMatcherToUpstreamResidue.hh>
#include <protocols/match/upstream/ScaffoldBuildPoint.hh>
#include <protocols/match/upstream/UpstreamBuilder.hh>
#include <core/pose/Pose.hh>

// Project headers
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>

#include <basic/Tracer.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <numeric/xyzVector.hh>

// C++ headers
#include <list>
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace match {
namespace downstream {

static THREAD_LOCAL basic::Tracer TR( "protocols.match.downstream.SecondaryMatcherToDownstreamResidue" );

SecondaryMatcherToDownstreamResidue::SecondaryMatcherToDownstreamResidue(
	core::pose::PoseCOP upstream_pose,
	Size geom_cst_id
) :
	parent( geom_cst_id ),
	upstream_pose_(upstream_pose)
{
	for ( core::Size ii=1; ii<=4; ++ii ) { catalytic_atoms_.push_back(0); }
}

SecondaryMatcherToDownstreamResidue::~SecondaryMatcherToDownstreamResidue() {}

DownstreamAlgorithmOP
SecondaryMatcherToDownstreamResidue::clone() const
{
	return DownstreamAlgorithmOP( new SecondaryMatcherToDownstreamResidue( *this ) );
}


std::list< Hit >
SecondaryMatcherToDownstreamResidue::build_hits_at_all_positions(
	Matcher & matcher
)
{
	prepare_for_hit_generation( matcher );

	utility::vector1< upstream::ScaffoldBuildPointCOP > const & my_build_points
		( matcher.per_constraint_build_points( geom_cst_id() ) );
	utility::vector1< std::list< Hit > > hits( my_build_points.size() );

	for ( Size ii = 1; ii <= geom_cst_id() - 1; ++ii ) {

		if ( ! matcher.representative_downstream_algorithm( ii )->generates_primary_hits() ) continue;

		prepare_for_hit_generation_for_geomcst( matcher, ii );

		/// nab the launch points for my target
		utility::vector1< upstream::ScaffoldBuildPointCOP > const & target_build_points
			( matcher.per_constraint_build_points( focused_geomcst_id_ ) );

		for ( Size jj = 1; jj <= target_build_points.size(); ++jj ) {
			if ( ! prepare_for_hit_generation_at_target_build_point( matcher, ii, *target_build_points[ jj ] ) ) {
				continue;
			}
			TR << "Secondary matching against geomcst " << ii << " hits from build point " << target_build_points[ jj ]->original_insertion_point() << std::endl;
			// Parallelize here.
#ifdef USE_OPENMP
			#pragma omp parallel for
#endif
			for ( Size kk = 1; kk <= my_build_points.size(); ++kk ) {
				//Most of the time we don't want to build residue at the same point
				//as the target build point, but in the case of backbone matching
				//if ( target_build_points[ jj ]->index() != my_build_points[ kk ]->index() ) {
				//std::cout << "APL DEBUG build points SecondaryMatcherToDownstreamResidue::build_hits_at_all_positions " << kk << " matcher.upstream_builder: " << matcher.upstream_builder( geom_cst_id() )() << std::endl;
				std::list< Hit > kk_hits = matcher.upstream_builder( geom_cst_id() )->build( * my_build_points[ kk ] );
				hits[ kk ].splice( hits[ kk ].end(), kk_hits );
				//}
			}
		}
	}

	std::list< Hit > all_hits;
	us_secmatch_hit_compare compare;
	for ( Size ii = 1; ii <= my_build_points.size(); ++ii ) {
		hits[ ii ].sort( compare );
		all_hits.splice( all_hits.end(), hits[ ii ] );
	}

	return all_hits;

}

void
SecondaryMatcherToDownstreamResidue::respond_to_primary_hitlist_change(
	Matcher & matcher,
	Size round_just_completed
)
{
	OccupiedSpaceHashOP occspace = matcher.occ_space_hash();

	if ( geom_cst_id() == round_just_completed ) {
		TR << "Finished round " << geom_cst_id() << " with " << matcher.hits( geom_cst_id() ).size();
		TR << " hits." << std::endl;
	}


	if ( geom_cst_id() != round_just_completed ) {
		TR << "Updating the occupied space grid with " << matcher.hits( geom_cst_id() ).size();
		TR << " hits from geometric constraint # " << geom_cst_id() << std::endl;
	}

	occspace->prepare_to_note_hits_for_completed_round();

	/// Prepare to clear the cobwebs.  Note which voxels in the occ_space_hash_
	/// could lead to matches, and which voxels could not possibly lead to matches.
	for ( Matcher::HitListConstIterator
			iter = matcher.hits( geom_cst_id() ).begin(),
			iter_end = matcher.hits( geom_cst_id() ).end();
			iter != iter_end; ++iter ) {
		occspace->note_hit_geometry( iter->second() );
	}

	occspace->drop_unsatisfied_voxels();

	occspace_rev_id_at_last_update_ = occspace->revision_id();

}


/// @details Drop hits that had previously seemed viable after another round completed;
/// during that round, certain previously occupied regions of 6D were not filled
/// with new hits.  Any previously-generated hit that falls into a region of 6D which is
/// no longer occupied should be elminated since it could not ever result in a match;
/// it is inviable.
void
SecondaryMatcherToDownstreamResidue::respond_to_peripheral_hitlist_change(
	Matcher & matcher
)
{
	OccupiedSpaceHashOP occspace = matcher.occ_space_hash();
	if ( occspace_rev_id_at_last_update_ == occspace->revision_id() ) {
		/// Nothing about the occupied space hash has changed since I last pruned
		/// my non-viable hits.  There are no more non-viable hits that need pruning.
		return;
	}

	Matcher::HitListIterator iter = matcher.hit_list_begin(   geom_cst_id() );
	Matcher::HitListIterator iter_end = matcher.hit_list_end( geom_cst_id() );

	Size drop_count( 0 );

	while ( iter != iter_end  ) {
		Matcher::HitListIterator iter_next = iter;
		++iter_next;
		if ( ! occspace->previous_round_geometry_still_matchable( iter->second() ) ) {
			matcher.erase_hit( *this, geom_cst_id(), iter );
			++drop_count;
		}
		iter = iter_next;
	}

	TR << "Erased " << drop_count << " round " << geom_cst_id();
	TR << " hits with ";
	TR << matcher.hits( geom_cst_id() ).size() << " hits remaining." << std::endl;

	occspace_rev_id_at_last_update_ = occspace->revision_id();

}


/// @brief Iterate across the hits from a particular upstream build point i
/// that were generated in a previous round, and see if the geometry of the
/// input upstream_residue has "satisfactory interactions" with the
/// hits from upstream-build-point i; if so, it appends a Hit to the hitlist
/// returned at the end of the method.  (Also, see comments for the
/// build_at_all_positions method.)
std::list< Hit >
SecondaryMatcherToDownstreamResidue::build(
	Size const scaffold_build_point_id,
	Size const upstream_conf_id,
	core::conformation::Residue const & upstream_residue
) const
{
	std::list< Hit > hits;
	Size const downstream_natoms = target_downstream_coords_->n_atoms_for_restype( 1 );
	core::conformation::Residue target_residue( *downstream_restype_, false );
	core::conformation::Residue target_residue_clash_check( *downstream_restype_, false );
	utility::vector1< core::Vector > ds_coords( ( target_downstream_coords_->get_ds_atom_ids_needed() ).size()  );

	// flo oct 09: moderate hack assuming that the ligand will always be last
	target_residue.seqpos( upstream_pose_->size() + 1);

	for ( Size ii = 1; ii <= target_downstream_coords_->n_rotamers_for_restype( 1 ); ++ii ) {
		/// Set the coordinates for only the atoms that I need to give to the evaluator
		for ( Size kk = 1; kk <= downstream_natoms; ++kk ) {
			target_residue.set_xyz(
				target_downstream_coords_->restype_atomno( 1, kk ),
				target_downstream_coords_->coord( 1, ii, kk ));
		}

		for ( EvaluatorSet::const_iterator eval_iter = respair_evaluators_.begin(),
				eval_iter_end = respair_evaluators_.end();
				eval_iter != eval_iter_end; ++eval_iter ) {

			if ( eval_iter->first->evaluate_residues( upstream_residue, target_residue ) ) {
				//detect collision between downstream evaluator atoms
				//and heavy sidechain atoms
				if ( ! are_colliding( upstream_residue, target_residue, downstream_atom_coordinates_needed_, catalytic_atoms_ ) ) {

					//get coordinates from the downstream otherwise build and store them
					if ( target_downstream_coords_->get_clash_checking( ii ) ) {
						ds_coords = target_downstream_coords_->get_coords_for_clash_check( ii );
					} else {
						get_dsbuilder()->coordinates_from_hit( target_downstream_coords_->hit(1,ii), target_downstream_coords_->get_ds_atom_ids_needed(), ds_coords );
						target_downstream_coords_->set_coords_for_clash_check( ii, ds_coords );
						target_downstream_coords_->set_clash_checking( ii );
					}
					//fills downstream residue with coordinates for clash checking atoms
					for ( Size ll = 1; ll <= ( target_downstream_coords_->get_ds_atom_ids_needed() ).size(); ++ll )  {
						target_residue_clash_check.set_xyz( target_downstream_coords_->get_ds_atom_ids_needed()[ ll ].atomno(), ds_coords[ ll ] );
					}

					if ( ! are_colliding( upstream_residue, target_residue_clash_check, target_downstream_coords_->get_ds_atom_ids_needed(), catalytic_atoms_ ) ) {
						Hit hit;
						hit.first()[ 1 ] = scaffold_build_point_id;
						hit.first()[ 2 ] = upstream_conf_id;
						hit.first()[ 3 ] = eval_iter->second;
						hit.first()[ 4 ] = (target_downstream_coords_->hit(1,ii)).downstream_conf_id();
						//flo sep '13 the below line is probably a bug, was old line, above is better
						//hit.first()[ 4 ] = 1; // do not store the downstream rotamer id or the focused_geomcst_id_; we never try to recover this placement
						hit.second()     = target_downstream_coords_->hit( 1, ii ).second(); // downstream coordinate from previous round
						hits.push_back( hit );
					}
				}
			}
		}
	}
	return hits;
}

bool
SecondaryMatcherToDownstreamResidue::upstream_only() const
{
	return false;
}

bool
SecondaryMatcherToDownstreamResidue::generates_primary_hits() const
{
	return false;
}


HitPtrListCOP
SecondaryMatcherToDownstreamResidue::hits_to_include_with_partial_match( match_dspos1 const & ) const
{
	HitPtrListCOP empty;
	utility_exit_with_message( "Cannot invoke SecondaryMatcherToDownstreamResidue::hits_to_include_with_partial_match()" );
	return empty;
}

SecondaryMatcherToDownstreamResidue::Size
SecondaryMatcherToDownstreamResidue::n_possible_hits_per_upstream_conformation() const
{
	return target_downstream_coords_->n_rots_total();
}


void
SecondaryMatcherToDownstreamResidue::set_downstream_restype( core::chemical::ResidueTypeCOP downstream_restype )
{
	downstream_restype_ = downstream_restype;
}


void
SecondaryMatcherToDownstreamResidue::set_focused_geomcst_id( Size focused_geomcst_id )
{
	focused_geomcst_id_ = focused_geomcst_id;
}

void
SecondaryMatcherToDownstreamResidue::add_evaluator(
	SecMatchResiduePairEvaluatorCOP evaluator,
	Size mcfi_id
)
{
	respair_evaluators_.push_back( std::make_pair( evaluator, mcfi_id ) );
}


void
SecondaryMatcherToDownstreamResidue::prepare_for_hit_generation(
	Matcher & matcher
)
{
	TR << "Preparing for hit generation" << std::endl;

	/// Initialize the target_geomcst_coords_ object;
	target_downstream_coords_ = TargetRotamerCoordsOP( new TargetRotamerCoords );
	Size n_target_restypes = 1;
	target_downstream_coords_->set_num_restypes( n_target_restypes );
	target_downstream_coords_->set_restype( n_target_restypes, downstream_restype_ );

	/// Initialize the other downstream algorithms for this geometric constraint
	/// so that they point all point to the same target_downstream_coords_ object
	std::list< DownstreamAlgorithmOP > const & dsalgs = matcher.nonconst_downstream_algorithms( geom_cst_id() );
	std::list< SecondaryMatcherToDownstreamResidueOP > secmatchers;
	for ( std::list< DownstreamAlgorithmOP >::const_iterator iter = dsalgs.begin(),
			iter_end = dsalgs.end(); iter != iter_end; ++iter ) {
		SecondaryMatcherToDownstreamResidueOP secmatcher =
			utility::pointer::dynamic_pointer_cast< SecondaryMatcherToDownstreamResidue> ( *iter );
		runtime_assert( secmatcher != 0 );
		if ( secmatcher.get() != this ) {
			secmatcher->set_target_rotamer_coords( target_downstream_coords_ );
		}
		secmatchers.push_back( secmatcher );
	}

	utility::vector1< bool > atom_required( downstream_restype_->natoms(), false );
	downstream_atom_coordinates_needed_.clear();

	bool any_evaluators_require_all_coords( false );
	for ( std::list< SecondaryMatcherToDownstreamResidueOP >::const_iterator
			secmatch_iter = secmatchers.begin(), secmatch_iter_end = secmatchers.end();
			secmatch_iter != secmatch_iter_end; ++secmatch_iter ) {

		for ( EvaluatorSet::const_iterator eval_iter = (*secmatch_iter)->respair_evaluators_.begin(),
				eval_iter_end = (*secmatch_iter)->respair_evaluators_.end();
				eval_iter != eval_iter_end; ++eval_iter ) {

			if ( eval_iter->first->require_all_target_residue_atom_coordinates() ) {
				any_evaluators_require_all_coords = true;
				break;
			}
		}
		if ( any_evaluators_require_all_coords ) break;
	}

	if ( any_evaluators_require_all_coords ) {
		std::fill( atom_required.begin(), atom_required.end(), true );
	} else {
		for ( Size ii = 1; ii <= downstream_restype_->natoms(); ++ii ) {

			for ( std::list< SecondaryMatcherToDownstreamResidueOP >::const_iterator
					secmatch_iter = secmatchers.begin(), secmatch_iter_end = secmatchers.end();
					secmatch_iter != secmatch_iter_end; ++secmatch_iter ) {

				for ( EvaluatorSet::const_iterator eval_iter = (*secmatch_iter)->respair_evaluators_.begin(),
						eval_iter_end = (*secmatch_iter)->respair_evaluators_.end();
						eval_iter != eval_iter_end; ++eval_iter ) {

					if ( eval_iter->first->require_target_atom_coordinate( ii ) ) {
						atom_required[ ii ] = true;
						break;
					}
				}
				if ( atom_required[ ii ] ) break;
			}
		} // for ii
	}

	target_downstream_coords_->set_required_atoms( 1, atom_required );
	for ( Size ii = 1; ii <= downstream_restype_->natoms(); ++ii ) {
		if ( atom_required[ ii ] ) downstream_atom_coordinates_needed_.push_back( core::id::AtomID( ii, 1 ) );
	}

	focused_geomcst_id_ = 0;
}

void
SecondaryMatcherToDownstreamResidue::prepare_for_hit_generation_for_geomcst(
	Matcher & matcher,
	Size target_geomcst_id
)
{
	focused_geomcst_id_ = target_geomcst_id;
	set_dsbuilder( matcher.downstream_builder( focused_geomcst_id_ ) );

	hits_for_focused_geomcst_and_build_point_begin_ = matcher.hits( focused_geomcst_id_ ).begin();
	hits_for_focused_geomcst_and_build_point_end_ = hits_for_focused_geomcst_and_build_point_begin_;
	hits_for_focused_geomcst_end_ = matcher.hits( focused_geomcst_id_ ).end();

	/// Make sure all the downstream algorithms for this geom_cst_id()
	/// know which geomcst is the source for this  next batch of hits.
	std::list< DownstreamAlgorithmOP > const & dsalgs = matcher.nonconst_downstream_algorithms( geom_cst_id() );
	for ( std::list< DownstreamAlgorithmOP >::const_iterator
			iter = dsalgs.begin(), iter_end = dsalgs.end();
			iter != iter_end; ++iter ) {
		if ( iter->get() != this ) {
			SecondaryMatcherToDownstreamResidueOP other = utility::pointer::dynamic_pointer_cast< SecondaryMatcherToDownstreamResidue > ( *iter );
			runtime_assert( other != 0 );
			other->set_focused_geomcst_id( focused_geomcst_id_ );
			//set downstram atom coords needed for all downstream algorithms
			//and also downstreambuilders that are needed to get coordinates
			//clash checking
			other->set_ds_coords_needed( downstream_atom_coordinates_needed_);
			other->set_dsbuilder( matcher.downstream_builder( focused_geomcst_id_ ) );
		}
	}

	toolbox::match_enzdes_util::LigandConformer lig_confs =  toolbox::match_enzdes_util::LigandConformer( * get_dsbuilder()->get_lig_conformers( 1 ) );
	utility::vector1< core::id::AtomID > downstream_atoms_for_clash_checking;

	for ( Size natoms_ds = 1; natoms_ds <= lig_confs.n_collision_check_atoms(); ++natoms_ds ) {
		Size id = lig_confs.collision_check_id_2_restype_id( natoms_ds );
		if ( ! ds_atom_present(id) ) {
			downstream_atoms_for_clash_checking.push_back( core::id::AtomID(id,1) );
		}
	}

	target_downstream_coords_->set_ds_atom_ids_needed(downstream_atoms_for_clash_checking);

}

bool
SecondaryMatcherToDownstreamResidue::ds_atom_present(
	Size index
) const {
	for ( Size ii = 1; ii <= downstream_atom_coordinates_needed_.size(); ++ii ) {
		if ( index == downstream_atom_coordinates_needed_[ii].atomno() ) return true;
	}
	return false;
}

bool
SecondaryMatcherToDownstreamResidue::prepare_for_hit_generation_at_target_build_point(
	Matcher & /*matcher*/,
	Size target_geomcst_id,
	upstream::ScaffoldBuildPoint const & target_build_point
)
{
	TR << "Preparing to examine geomcst-" << target_geomcst_id << "'s hits built from scaffold build point " << target_build_point.original_insertion_point() << std::endl;


	Size target_build_id = target_build_point.index();

	/// Find the range of iterators that cover the hits for the upstream matcher
	while ( hits_for_focused_geomcst_and_build_point_begin_ != hits_for_focused_geomcst_end_ ) {
		if ( hits_for_focused_geomcst_and_build_point_begin_->scaffold_build_id() == target_build_id ) {
			break;
		} else if ( hits_for_focused_geomcst_and_build_point_begin_->scaffold_build_id() > target_build_id ) {
			break;
		}
		++hits_for_focused_geomcst_and_build_point_begin_;
	}
	hits_for_focused_geomcst_and_build_point_end_ = hits_for_focused_geomcst_and_build_point_begin_;
	Size count_target_hits( 0 );
	while ( hits_for_focused_geomcst_and_build_point_end_ != hits_for_focused_geomcst_end_ ) {
		if ( hits_for_focused_geomcst_and_build_point_end_->scaffold_build_id() != target_build_id ) {
			break;
		}
		++hits_for_focused_geomcst_and_build_point_end_;
		++count_target_hits;
	}

	if ( hits_for_focused_geomcst_and_build_point_begin_ == hits_for_focused_geomcst_and_build_point_end_ ) {
		TR << "No geomcst-" << focused_geomcst_id_ << " hits built from scaffold build point " << target_build_point.original_insertion_point() << std::endl;
		return false;
	}

	/// Now rebuild the coordinates for the hits that have been selected.
	target_downstream_coords_->set_num_target_rotamers( 1, count_target_hits );
	target_downstream_coords_->set_clash_check_types( count_target_hits );

	// DownstreamBuilderOP dsbuilder = matcher.downstream_builder( focused_geomcst_id_ );
	//set_dsbuilder( matcher.downstream_builder( focused_geomcst_id_ ) );
	runtime_assert( get_dsbuilder() != 0 );

	core::conformation::Residue dsrescoords( *downstream_restype_, false );
	utility::vector1< core::Vector > coords( downstream_atom_coordinates_needed_.size() );

	Size count_ds_hits( 0 );
	for ( Matcher::HitListConstIterator iter = hits_for_focused_geomcst_and_build_point_begin_,
			iter_end = hits_for_focused_geomcst_and_build_point_end_; iter != iter_end; ++iter ) {
		++count_ds_hits;
		get_dsbuilder()->coordinates_from_hit( *iter, downstream_atom_coordinates_needed_, coords );
		for ( Size ii = 1; ii <= downstream_atom_coordinates_needed_.size(); ++ii ) {
			dsrescoords.set_xyz( downstream_atom_coordinates_needed_[ ii ].atomno(), coords[ ii ] );
		}
		target_downstream_coords_->set_coordinates_for_rotamer(
			1, count_ds_hits,
			*iter, dsrescoords );
	}

	/// This would also be the appropriate time to initialize a pruning object to
	/// quickly reject rotamers in the UpstreamBuilder.
	TR << "Examining " << count_target_hits << " hits built from scaffold build point "
		<< target_build_point.original_insertion_point() << " representing " << target_downstream_coords_->n_rots_total()
		<< " unique upstream rotamers" << std::endl;

	return true;
}

void
SecondaryMatcherToDownstreamResidue::set_target_rotamer_coords(
	TargetRotamerCoordsOP target_downstream_coords
)
{
	target_downstream_coords_ = target_downstream_coords;
}

}
}
}
