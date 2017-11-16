// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/downstream/SecondaryMatcherToUpstreamResidue.cc
/// @brief  Class implementation for secondary matcher that generates upstream-only hits
///         matching the geometry of one upstream residue with another upstream residue
///         generated in a previous round.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Florian Richter (flosopher@gmail.com)


// Unit headers
#include <protocols/match/downstream/SecondaryMatcherToUpstreamResidue.hh>

// Package headers
#include <protocols/match/Matcher.hh>
#include <protocols/match/Hit.hh>
#include <protocols/match/downstream/DownstreamAlgorithm.hh>
#include <protocols/match/downstream/SecMatchResiduePairEvaluator.hh>
#include <protocols/match/upstream/ScaffoldBuildPoint.hh>
#include <protocols/match/upstream/UpstreamBuilder.hh>

// Project headers
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>

#include <basic/Tracer.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <list>
#include <string>

// Boost headers
#include <boost/unordered_map.hpp>

#include <utility/OrderedTuple.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace match {
namespace downstream {

static basic::Tracer TR( "protocols.match.downstream.SecondaryMatcherToUpstreamResidue" );

struct hash_upstream_hit
{
	core::Size operator() ( SecondaryMatcherToUpstreamResidue::Size2Tuple const & uphit ) const {
		return ( 271 * uphit.data()[ 1 ] + uphit.data()[ 2 ] ) % 1699;
	}
};


SecondaryMatcherToUpstreamResidue::SecondaryMatcherToUpstreamResidue(
	Size geom_cst_id
) :
	parent( geom_cst_id ),
	target_geomcst_id_( 0 ),
	count_rotamers_per_target_restype_( false ),
	last_seen_restype_index_( 0 ),
	count_rotamer_for_lastseen_restype_( 0 ),
	n_viable_target_hits_since_last_pruning_( 0 )
{
}

SecondaryMatcherToUpstreamResidue::~SecondaryMatcherToUpstreamResidue() {}

DownstreamAlgorithmOP
SecondaryMatcherToUpstreamResidue::clone() const
{
	return DownstreamAlgorithmOP( new SecondaryMatcherToUpstreamResidue( *this ) );
}

/// @details This DownstreamAlgorithm structures it's iteration over the target hits
/// from previous rounds as follows:
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
std::list< Hit >
SecondaryMatcherToUpstreamResidue::build_hits_at_all_positions(
	Matcher & matcher
)
{
	prepare_for_hit_generation( matcher );

	//kui 102109
	//smUR_pose_build_resids_ = matcher.get_pose_build_resids();

	/// nab the launch points for both my target and myself
	utility::vector1< upstream::ScaffoldBuildPointCOP > const & target_build_points
		( matcher.per_constraint_build_points( target_geomcst_id_ ) );
	utility::vector1< upstream::ScaffoldBuildPointCOP > const & my_build_points
		( matcher.per_constraint_build_points( geom_cst_id() ) );

	utility::vector1< std::list< Hit > > hits( my_build_points.size() );

	for ( Size ii = 1; ii <= target_build_points.size(); ++ii ) {
		if ( ! prepare_for_hit_generation_at_target_build_point( matcher, *target_build_points[ ii ] ) ) {
			continue;
		}

		//TR << "Secondary matching against geomcst " << target_geomcst_id_ << " hits from build point " << target_build_points[ ii ]->index() << std::endl;
		//TR << "Secondary matching against geomcst " << target_geomcst_id_ << " hits from protein build point " << target_build_points[ ii ]->original_insertion_point() << std::endl;

		Size ii_total_found_matches( 0 );
#ifdef USE_OPENMP
		#pragma omp parallel for reduction( + : ii_total_found_matches )
#endif
		for ( Size jj = 1; jj <= my_build_points.size(); ++jj ) {
			//Most of the time we don't want to build residue at the same point
			//as the target build point, but in the case of backbone matching
			//if ( target_build_points[ ii ]->index() != my_build_points[ jj ]->index() ) {
			///TR << " Building sec. match res on protein build point:  " << my_build_points[ jj ]->original_insertion_point() << std::endl;
			std::list< Hit > jj_hits = matcher.upstream_builder( geom_cst_id() )->build( * my_build_points[ jj ] );
			ii_total_found_matches += jj_hits.size();
			hits[ jj ].splice( hits[ jj ].end(), jj_hits );
			//}
		}
		TR << "Secondary matching against geomcst " << target_geomcst_id_ << " hits from protein build point " << target_build_points[ ii ]->original_insertion_point() << " produced " << ii_total_found_matches << " hits " << std::endl;
	}

	std::list< Hit > all_hits;
	us_secmatch_hit_compare compare;
	for ( Size ii = 1; ii <= my_build_points.size(); ++ii ) {
		hits[ ii ].sort( compare );
		all_hits.splice( all_hits.end(), hits[ ii ] );
	}

	return all_hits;

}


/// @brief Prune hits away from the target_geomcst's hit list following a change to the
/// hits for my geom_cst_id().  Pruning hits from the target_geomcst's hit list will
/// trigger a round of peripheral-hitlist-change responses.
void
SecondaryMatcherToUpstreamResidue::respond_to_primary_hitlist_change(
	Matcher & matcher,
	Size /*round_just_completed*/
)
{

	boost::unordered_map< Size2Tuple, bool, hash_upstream_hit > targets_hits_matched;
	Matcher::HitList const & my_hits = matcher.hits( geom_cst_id() );
	for ( Matcher::HitListConstIterator iter = my_hits.begin(), iter_end = my_hits.end();
			iter != iter_end; ++iter ) {
		Size2 target_hit;
		target_hit[ 1 ] = static_cast< Size > ( iter->second()[ 1 ] );
		target_hit[ 2 ] = static_cast< Size > ( iter->second()[ 2 ] );
		boost::unordered_map< Size2Tuple, bool, hash_upstream_hit >::const_iterator
			hash_iter = targets_hits_matched.find( target_hit );
		if ( hash_iter == targets_hits_matched.end() ) {
			targets_hits_matched[ target_hit ] = true;
		}
	}

	Matcher::HitListIterator target_hit_iter = matcher.hit_list_begin( target_geomcst_id_ );
	Matcher::HitListIterator target_hit_iter_end = matcher.hit_list_end( target_geomcst_id_ );

	Size drop_count( 0 );
	while ( target_hit_iter != target_hit_iter_end ) {
		Matcher::HitListIterator target_hit_iter_next = target_hit_iter;
		++target_hit_iter_next;

		Size2 target_hit;
		target_hit[ 1 ] = static_cast< Size > ( target_hit_iter->first()[ 1 ] );
		target_hit[ 2 ] = static_cast< Size > ( target_hit_iter->first()[ 2 ] );
		boost::unordered_map< Size2Tuple, bool, hash_upstream_hit >::const_iterator
			hash_iter = targets_hits_matched.find( target_hit );
		if ( hash_iter == targets_hits_matched.end() ) {
			matcher.erase_hit( *this, target_geomcst_id_, target_hit_iter );
			++drop_count;
		}

		target_hit_iter = target_hit_iter_next;
	}

	n_viable_target_hits_since_last_pruning_ = matcher.hits( target_geomcst_id_ ).size();

	TR << "Erased " << drop_count << " round " << target_geomcst_id_;
	TR << " hits with ";
	TR << n_viable_target_hits_since_last_pruning_ << " hits remaining." << std::endl;
}


/// @brief Remove my hits if my target_geomcst's hit list has been shortened.  This
/// will not trigger a round of peripheral-hitlist-change responses.
void
SecondaryMatcherToUpstreamResidue::respond_to_peripheral_hitlist_change(
	Matcher & matcher
)
{
	if ( n_viable_target_hits_since_last_pruning_
			== matcher.hits( target_geomcst_id_ ).size() ) {
		return;
	}

	boost::unordered_map< Size2Tuple, bool, hash_upstream_hit > target_hit_hash;
	Matcher::HitList const & target_hits( matcher.hits( target_geomcst_id_ ));
	for ( Matcher::HitListConstIterator iter = target_hits.begin(), iter_end = target_hits.end();
			iter != iter_end; ++iter ) {
		Size2 target_hit;
		target_hit[ 1 ] = iter->first()[ 1 ];
		target_hit[ 2 ] = iter->first()[ 2 ];
		boost::unordered_map< Size2Tuple, bool, hash_upstream_hit >::const_iterator
			hash_iter = target_hit_hash.find( target_hit );
		if ( hash_iter == target_hit_hash.end() ) {
			target_hit_hash[ target_hit ] = true;
		}
	}

	Matcher::HitListIterator hit_iter = matcher.hit_list_begin( geom_cst_id() );
	Matcher::HitListIterator hit_iter_end = matcher.hit_list_end( geom_cst_id() );

	Size drop_count( 0 );
	while ( hit_iter != hit_iter_end ) {
		Matcher::HitListIterator hit_iter_next = hit_iter;
		++hit_iter_next;

		Size2 target_hit;
		target_hit[ 1 ] = static_cast< Size > ( hit_iter->second()[ 1 ] );
		target_hit[ 2 ] = static_cast< Size > ( hit_iter->second()[ 2 ] );
		boost::unordered_map< Size2Tuple, bool, hash_upstream_hit >::const_iterator
			hash_iter = target_hit_hash.find( target_hit );
		if ( hash_iter == target_hit_hash.end() ) {
			matcher.erase_hit( *this, geom_cst_id(), hit_iter );
			++drop_count;
		}

		hit_iter = hit_iter_next;
	}

	n_viable_target_hits_since_last_pruning_ = matcher.hits( target_geomcst_id_ ).size();

	TR << "Erased " << drop_count << " round " << geom_cst_id();
	TR << " hits with ";
	TR << matcher.hits( geom_cst_id() ).size() << " hits remaining." << std::endl;

}


/// @brief Iterate across the hits from a particular upstream build point i
/// that were generated in a previous round, and see if the geometry of the
/// input upstream_residue has "satisfactory interactions" with the
/// hits from upstream-build-point i; if so, it appends a Hit to the hitlist
/// returned at the end of the method.  (Also, see comments for the
/// build_at_all_positions method.)
std::list< Hit >
SecondaryMatcherToUpstreamResidue::build(
	Size const scaffold_build_point_id,
	Size const upstream_conf_id,
	core::conformation::Residue const & upstream_residue
) const
{
	std::list< Hit > hits;
	for ( Size ii = 1; ii <= target_geomcst_coords_->n_restypes(); ++ii ) {
		core::conformation::Residue target_residue( *target_geomcst_coords_->restype( ii ), false );
		Size const ii_natoms = target_geomcst_coords_->n_atoms_for_restype( ii );
		core::chemical::ResidueTypeCOP us_res_type = upstream_residue.type_ptr();
		//TR << " sec. matched  residue type:  " << us_res_type->name() << std::endl;

		for ( Size jj = 1; jj <= target_geomcst_coords_->n_rotamers_for_restype( ii ); ++jj ) {
			//TR << "number of rotamers for restype " << target_geomcst_coords_->n_rotamers_for_restype( ii ) << std::endl;


			target_residue.seqpos(smUR_pose_build_resids_[target_geomcst_coords_->hit( ii, jj ).scaffold_build_id()]);
			/// Set the coordinates for only the atoms that I need to give to the evaluator
			for ( Size kk = 1; kk <= ii_natoms; ++kk ) {
				target_residue.set_xyz(
					target_geomcst_coords_->restype_atomno( ii, kk ),
					target_geomcst_coords_->coord( ii, jj, kk ));
			}

			for ( EvaluatorSet::const_iterator eval_iter = respair_evaluators_[ ii ].begin(),
					eval_iter_end = respair_evaluators_[ ii ].end();
					eval_iter != eval_iter_end; ++eval_iter ) {

				if ( eval_iter->first->evaluate_residues( upstream_residue, target_residue ) ) {
					//Do a clash check by passing the current hit, upstream residue
					//and a dsbuilder that contains the clash info between the
					//upstream restype and the downstream residue
					//if ( ! are_colliding( target_geomcst_coords_->hit(ii,jj), upstream_residue, dsbuilders_.find(  us_res_type->name() )->second ) ){
					Hit hit;
					hit.first()[ 1 ] = scaffold_build_point_id;
					hit.first()[ 2 ] = upstream_conf_id;
					hit.first()[ 3 ] = eval_iter->second; // mcfi ID
					hit.first()[ 4 ] = 1;
					hit.second()[ 1 ] = target_geomcst_coords_->hit( ii, jj ).scaffold_build_id();
					hit.second()[ 2 ] = target_geomcst_coords_->hit( ii, jj ).upstream_conf_id();
					hit.second()[ 3 ] = target_geomcst_coords_->hit( ii, jj ).external_geom_id();
					hit.second()[ 4 ] = 0.0;
					hit.second()[ 5 ] = 0.0;
					hit.second()[ 6 ] = 0.0;
					hits.push_back( hit );
					//TR << " ADDING A HIT! " << "rotamers for restype " << jj << std::endl;
					//}
				}
			}
		}
	}
	return hits;
}

bool
SecondaryMatcherToUpstreamResidue::upstream_only() const
{
	return true;
}

bool
SecondaryMatcherToUpstreamResidue::generates_primary_hits() const
{
	return false;
}


/// @brief Prepare a map between upstream hits of the target-geomcst and
/// a list of Hit const *'s of this geom_cst_id(). This map will be used
/// in the function hits_to_include_with_partial_match.
void
SecondaryMatcherToUpstreamResidue::prepare_for_match_enumeration( Matcher const & matcher )
{
	my_hits_for_target_hit_map_.clear();
	Matcher::HitList const & hits( matcher.hits( geom_cst_id() ));
	for ( Matcher::HitListConstIterator iter = hits.begin(), iter_end = hits.end();
			iter != iter_end; ++iter ) {
		Size2 target_hit;
		target_hit[ 1 ] = static_cast< Size > ( iter->second()[ 1 ] );
		target_hit[ 2 ] = static_cast< Size > ( iter->second()[ 2 ] );
		std::map< Size2Tuple, HitPtrListOP >::const_iterator list_iter =
			my_hits_for_target_hit_map_.find( target_hit );
		if ( list_iter == my_hits_for_target_hit_map_.end() ) {
			my_hits_for_target_hit_map_[ target_hit ] = utility::pointer::shared_ptr<class numeric::kdtree::WrappedPrimitive<class std::list<const class protocols::match::Hit *, class std::allocator<const class protocols::match::Hit *> > > >( new HitPtrList );
			list_iter = my_hits_for_target_hit_map_.find( target_hit );
			//std::cout << "add to map " << target_hit[ 1 ] << " " << target_hit[ 2 ] << std::endl;
		}
		list_iter->second->val().push_back( & (*iter ) ); /// add a pointer to the matcher's hit
	}
}

/// @brief Return the set of hits to be iterated across
HitPtrListCOP
SecondaryMatcherToUpstreamResidue::hits_to_include_with_partial_match(
	match_dspos1 const & m
) const
{
	HitPtrListCOP hitptrlist( 0 );
	Size2 target_hit;
	target_hit[ 1 ] = static_cast< Size > ( m.upstream_hits[ target_geomcst_id_ ].scaffold_build_id() );
	target_hit[ 2 ] = static_cast< Size > ( m.upstream_hits[ target_geomcst_id_ ].upstream_conf_id() );
	std::map< Size2Tuple, HitPtrListOP >::const_iterator list_iter =
		my_hits_for_target_hit_map_.find( target_hit );
	if ( list_iter != my_hits_for_target_hit_map_.end() ) {
		hitptrlist = list_iter->second;
	} else {
		std::cerr << "fetch failed! " << target_hit[ 1 ] << " " << target_hit[ 2 ] << std::endl;
		std::cerr << "something weird happened!" << std::endl;
	}
	return hitptrlist;
}

SecondaryMatcherToUpstreamResidue::Size
SecondaryMatcherToUpstreamResidue::n_possible_hits_per_upstream_conformation() const
{
	return target_geomcst_coords_->n_rots_total();
}

//void
//SecondaryMatcherToUpstreamResidue::set_match_restype( core::chemical::ResidueTypeCOP match_restype )
//{
//
//}

void
SecondaryMatcherToUpstreamResidue::set_target_geomcst_id( Size target_geomcst_id )
{
	target_geomcst_id_ = target_geomcst_id;
}

void
SecondaryMatcherToUpstreamResidue::add_target_restype( core::chemical::ResidueTypeCOP target_restype )
{
	if ( target_restype_index_map_.find( target_restype ) == target_restype_index_map_.end() ) {
		target_restypes_.push_back( target_restype );
		target_restype_index_map_[ target_restype ] = target_restypes_.size();
	}
}

void
SecondaryMatcherToUpstreamResidue::add_evaluator_for_target_restype(
	core::chemical::ResidueTypeCOP  target_restype,
	SecMatchResiduePairEvaluatorCOP evaluator,
	Size                            mcfi_id_for_evaluator
)
{
	runtime_assert( target_restype_index_map_.find( target_restype ) != target_restype_index_map_.end() );
	if ( target_restypes_.size() > respair_evaluators_.size() ) {
		respair_evaluators_.resize( target_restypes_.size() );
	}
	respair_evaluators_[ target_restype_index_map_[ target_restype ] ].push_back(
		std::make_pair( evaluator, mcfi_id_for_evaluator ));
}

void
SecondaryMatcherToUpstreamResidue::process_hit(
	Hit const & hit,
	core::conformation::Residue const & upstream_conformation
)
{
	if ( count_rotamers_per_target_restype_ ) {
		count_rotamer( upstream_conformation );
	} else {
		store_rotamer_coords( hit, upstream_conformation );
	}
}

void
SecondaryMatcherToUpstreamResidue::prepare_for_hit_generation(
	Matcher & matcher
)
{
	TR << "Preparing for hit generation" << std::endl;
	target_hits_for_focused_build_point_begin_ = matcher.hits( target_geomcst_id_ ).begin();
	target_hits_for_focused_build_point_end_ = target_hits_for_focused_build_point_begin_;
	target_hits_end_ = matcher.hits( target_geomcst_id_ ).end();

	/// Initialize the target_geomcst_coords_ object;
	target_geomcst_coords_ = TargetRotamerCoordsOP( new TargetRotamerCoords );
	upstream::UpstreamBuilderCOP usbuilder = matcher.upstream_builder( target_geomcst_id_ );

	Size n_target_restypes = usbuilder->n_restypes_to_build();
	target_geomcst_coords_->set_num_restypes( n_target_restypes );
	n_rotamers_per_target_restype_.resize( n_target_restypes );

	std::list< DownstreamAlgorithmOP > const & dsalgs = matcher.nonconst_downstream_algorithms( geom_cst_id() );
	std::list< SecondaryMatcherToUpstreamResidueOP > secmatch_algs;
	for ( std::list< DownstreamAlgorithmOP >::const_iterator iter = dsalgs.begin(),
			iter_end = dsalgs.end(); iter != iter_end; ++iter ) {
		SecondaryMatcherToUpstreamResidueOP secmatcher =
			utility::pointer::dynamic_pointer_cast< SecondaryMatcherToUpstreamResidue > ( *iter );
		runtime_assert( secmatcher != 0 );
		secmatcher->smUR_pose_build_resids_ = matcher.get_pose_build_resids();
		secmatch_algs.push_back( secmatcher );
		secmatcher->reorder_restypes( *usbuilder );
	}

	for ( Size ii = 1; ii <= n_target_restypes; ++ii ) {
		core::chemical::ResidueTypeCOP iirestype = usbuilder->restype( ii );
		target_geomcst_coords_->set_restype( ii, iirestype );
		utility::vector1< bool > atom_required( iirestype->natoms(), false );

		bool any_evaluator_requires_all_atoms( false );
		for ( std::list< SecondaryMatcherToUpstreamResidueOP >::const_iterator iter = secmatch_algs.begin(),
				iter_end = secmatch_algs.end(); iter != iter_end; ++iter ) {

			for ( EvaluatorSet::const_iterator eval_iter = (*iter)->respair_evaluators_[ ii ].begin(),
					eval_iter_end = (*iter)->respair_evaluators_[ ii ].end();
					eval_iter != eval_iter_end; ++eval_iter ) {
				if ( eval_iter->first->require_all_target_residue_atom_coordinates() ) {
					any_evaluator_requires_all_atoms = true;
					break;
				}
			}
			if ( any_evaluator_requires_all_atoms ) break;
		}
		if ( any_evaluator_requires_all_atoms ) {
			std::fill( atom_required.begin(), atom_required.end(), true );
		} else {
			for ( Size jj = 1; jj <= iirestype->natoms(); ++jj ) {
				for ( std::list< SecondaryMatcherToUpstreamResidueOP >::const_iterator iter = secmatch_algs.begin(),
						iter_end = secmatch_algs.end(); iter != iter_end; ++iter ) {

					for ( EvaluatorSet::const_iterator eval_iter = (*iter)->respair_evaluators_[ ii ].begin(),
							eval_iter_end = (*iter)->respair_evaluators_[ ii ].end();
							eval_iter != eval_iter_end; ++eval_iter ) {

						if ( eval_iter->first->require_target_atom_coordinate( jj ) ) {
							atom_required[ jj ] = true;
							break;
						}
					}
					if ( atom_required[ jj ] ) break;
				}
			}
		}
		target_geomcst_coords_->set_required_atoms( ii, atom_required );
	}
}

bool
SecondaryMatcherToUpstreamResidue::prepare_for_hit_generation_at_target_build_point(
	Matcher & matcher,
	upstream::ScaffoldBuildPoint const & target_build_point
)
{
	TR << "Preparing to examine geomcst-" << target_geomcst_id_ << "'s hits built from scaffold build point " << target_build_point.original_insertion_point() << std::endl;

	Size target_build_id = target_build_point.index();

	/// Find the range of iterators that cover the hits for the upstream matcher
	while ( target_hits_for_focused_build_point_begin_ != target_hits_end_ ) {
		if ( target_hits_for_focused_build_point_begin_->scaffold_build_id() == target_build_id ) {
			break;
		} else if ( target_hits_for_focused_build_point_begin_->scaffold_build_id() > target_build_id ) {
			break;
		}
		++target_hits_for_focused_build_point_begin_;
	}
	target_hits_for_focused_build_point_end_ = target_hits_for_focused_build_point_begin_;
	Size count_target_hits( 0 );
	while ( target_hits_for_focused_build_point_end_ != target_hits_end_ ) {
		if ( target_hits_for_focused_build_point_end_->scaffold_build_id() != target_build_id ) {
			break;
		}
		++target_hits_for_focused_build_point_end_;
		++count_target_hits;
	}

	if ( target_hits_for_focused_build_point_begin_ == target_hits_for_focused_build_point_end_ ) {
		TR << "No geomcst-" << target_geomcst_id_ << " hits built from scaffold build point " << target_build_point.original_insertion_point() << std::endl;
		return false;
	}

	std::list< Hit > unique_upstream_hits;
	Size last_conf_id( 0 );
	for ( Matcher::HitListConstIterator iter = target_hits_for_focused_build_point_begin_;
			iter != target_hits_for_focused_build_point_end_; ++iter ) {
		//get_all_hits_for_clash_cheicking
		if ( last_conf_id != iter->upstream_conf_id() ) {
			unique_upstream_hits.push_back( *iter );
			last_conf_id = iter->upstream_conf_id();
		}
	}

	/// first pass over the hits; count the number of hits per residue type.
	count_rotamers_per_target_restype_ = true;
	last_seen_restype_index_ = 1;
	last_seen_restype_ = target_restypes_[ 1 ];
	std::fill( n_rotamers_per_target_restype_.begin(), n_rotamers_per_target_restype_.end(), 0 );
	SecondaryMatchUpstreamResProcessor proc( *this );
	matcher.upstream_builder( target_geomcst_id_ )->recover_hits(
		unique_upstream_hits.begin(), unique_upstream_hits.end(),
		target_build_point, proc );

	target_geomcst_coords_->set_num_target_rotamers( n_rotamers_per_target_restype_ );

	/// second pass over the hits; store the rotamers in the TargetRotamerCoords object
	count_rotamers_per_target_restype_ = false;
	last_seen_restype_index_ = 1;
	last_seen_restype_ = target_restypes_[ 1 ];
	count_rotamer_for_lastseen_restype_ = 0;
	matcher.upstream_builder( target_geomcst_id_ )->recover_hits(
		unique_upstream_hits.begin(), unique_upstream_hits.end(),
		target_build_point, proc );

	/// Now update all of the downstream-algorithms for this geom_cst_id() so that they're all pointing
	/// to the same TargetRotamerCoords object.

	std::list< DownstreamAlgorithmOP > const & dsalgs = matcher.nonconst_downstream_algorithms( geom_cst_id() );
	for ( std::list< DownstreamAlgorithmOP >::const_iterator
			iter = dsalgs.begin(), iter_end = dsalgs.end();
			iter != iter_end; ++iter ) {
		if ( iter->get() != this ) {
			SecondaryMatcherToUpstreamResidueOP other =
				utility::pointer::dynamic_pointer_cast< SecondaryMatcherToUpstreamResidue > ( *iter );
			runtime_assert( other != 0 );
			//TR << "SecondaryMatcherToUpstreamResidue * other" << other << std::endl;
			other->set_target_rotamer_coords( target_geomcst_coords_ );
		}
	}

	/// This would also be the appropriate time to initialize a pruning object to
	/// quickly reject rotamers in the UpstreamBuilder.
	TR << "Examining " << count_target_hits << " hits built from scaffold build point "
		<< target_build_point.original_insertion_point() << " representing " << target_geomcst_coords_->n_rots_total()
		<< " unique upstream rotamers" << std::endl;

	return true;
}

void
SecondaryMatcherToUpstreamResidue::set_target_rotamer_coords(
	TargetRotamerCoordsOP target_geomcst_coords
)
{
	target_geomcst_coords_ = target_geomcst_coords;
}

void
SecondaryMatcherToUpstreamResidue::count_rotamer(
	core::conformation::Residue const & upstream_conformation
)
{
	//std::cout << "count_rotamer: " << upstream_conformation.name() << std::endl;
	if ( & upstream_conformation.type() != last_seen_restype_.get() ) {
		++last_seen_restype_index_;
		while ( last_seen_restype_index_ <= target_restypes_.size() ) {
			if ( target_restypes_[ last_seen_restype_index_ ].get() == & upstream_conformation.type() ) break;
			++last_seen_restype_index_;
		}
		runtime_assert( last_seen_restype_index_ <= target_restypes_.size() );
		last_seen_restype_ = upstream_conformation.type_ptr();
	}
	++n_rotamers_per_target_restype_[ last_seen_restype_index_ ];
}

void
SecondaryMatcherToUpstreamResidue::store_rotamer_coords(
	Hit const & hit,
	core::conformation::Residue const & upstream_conformation
)
{
	//std::cout << "store_rotamer_coords: " << upstream_conformation.name() << std::endl;
	if ( & upstream_conformation.type() == last_seen_restype_.get() ) {
		++count_rotamer_for_lastseen_restype_;
		runtime_assert( count_rotamer_for_lastseen_restype_ <= n_rotamers_per_target_restype_[ last_seen_restype_index_ ] );
	} else {
		++last_seen_restype_index_;
		while ( last_seen_restype_index_ <= target_restypes_.size() ) {
			if ( target_restypes_[ last_seen_restype_index_ ].get() == & upstream_conformation.type() ) break;
			++last_seen_restype_index_;
		}
		runtime_assert( last_seen_restype_index_ <= target_restypes_.size() );
		last_seen_restype_ = upstream_conformation.type_ptr();
		count_rotamer_for_lastseen_restype_ = 1;

	}

	target_geomcst_coords_->set_coordinates_for_rotamer(
		last_seen_restype_index_, count_rotamer_for_lastseen_restype_,
		hit, upstream_conformation );
}

void SecondaryMatcherToUpstreamResidue::reorder_restypes(
	upstream::UpstreamBuilder const & usbuilder
)
{
	utility::vector1< Size > old_2_new( target_restypes_.size(), 0 );
	//TR << "TARGET RESTYPES: "<< target_restypes_.size() << " UPSTREAMBUILDER RESTYPES: " << usbuilder.n_restypes_to_build() << std::endl;
	runtime_assert( target_restypes_.size() == usbuilder.n_restypes_to_build() );
	for ( Size ii = 1; ii <= usbuilder.n_restypes_to_build(); ++ii ) {
		std::map< core::chemical::ResidueTypeCOP, Size >::const_iterator
			olditer = target_restype_index_map_.find( usbuilder.restype( ii ) );
		//TR << "TARGET RESTYPES: " << usbuilder.restype( ii )->name() << std::endl;
		runtime_assert( olditer != target_restype_index_map_.end() );
		old_2_new[ olditer->second ] = ii;
	}
	for ( Size ii = 1; ii <= usbuilder.n_restypes_to_build(); ++ii ) {
		runtime_assert( old_2_new[ ii ] != 0 );
	}

	utility::vector1< core::chemical::ResidueTypeCOP >  old_target_restypes( target_restypes_ );
	utility::vector1< EvaluatorSet > old_respair_evaluators( respair_evaluators_ );

	for ( Size ii = 1; ii <= usbuilder.n_restypes_to_build(); ++ii ) {
		target_restypes_[    old_2_new[ ii ] ] = old_target_restypes[    ii ];
		respair_evaluators_[ old_2_new[ ii ] ] = old_respair_evaluators[ ii ];
		target_restype_index_map_[ target_restypes_[ old_2_new[ ii ] ] ] = old_2_new[ ii ];
	}

}

// TargetRotamerCoords class
TargetRotamerCoords::TargetRotamerCoords() :
	n_rots_total_( 0 )
{
}

TargetRotamerCoords::~TargetRotamerCoords() {}

void TargetRotamerCoords::set_num_restypes( Size n_restypes )
{
	target_restypes_.resize( n_restypes );
	atom_ids_for_coordinates_.resize( n_restypes );
	coords_.resize( n_restypes );
	hit_data_.resize( n_restypes );
}

void TargetRotamerCoords::set_restype(
	Size restype_index,
	core::chemical::ResidueTypeCOP restype )
{
	target_restypes_[ restype_index ] = restype;
}

void TargetRotamerCoords::set_required_atoms(
	Size restype_index,
	utility::vector1< bool > const & atom_required
)
{
	Size count_required( 0 );
	runtime_assert( atom_required.size() == target_restypes_[ restype_index ]->natoms() );

	for ( Size ii = 1; ii <= atom_required.size(); ++ii ) if ( atom_required[ ii ] ) ++count_required;
	atom_ids_for_coordinates_[ restype_index ].resize( count_required );
	count_required = 0;
	for ( Size ii = 1; ii <= atom_required.size(); ++ii ) {
		if ( atom_required[ ii ] ) {
			++count_required;
			atom_ids_for_coordinates_[ restype_index ][ count_required ] = ii;
			//std::cout << "Requiring atom " << target_restypes_[ restype_index ]->atom_name( ii ) << std::endl;
		}
	}
}

void TargetRotamerCoords::set_num_target_rotamers(
	utility::vector1< Size > const & n_rotamers_per_target_restype
)
{
	runtime_assert( n_rotamers_per_target_restype.size() == target_restypes_.size() );
	n_rots_total_ = 0;
	for ( Size ii = 1; ii <= n_rotamers_per_target_restype.size(); ++ii ) {
		n_rots_total_ += n_rotamers_per_target_restype[ ii ];
		if ( n_rotamers_per_target_restype[ ii ] != 0 ) {
			//std::cout << "Setting n_rotamers for " << target_restypes_[ ii ]->name() << " " << atom_ids_for_coordinates_[ ii ].size() << " " << n_rotamers_per_target_restype[ ii ] << std::endl;
			coords_[ ii ].dimension( atom_ids_for_coordinates_[ ii ].size(), n_rotamers_per_target_restype[ ii ] );
			hit_data_[ ii ].resize( n_rotamers_per_target_restype[ ii ] );
		} else {
			coords_[ ii ].clear();
			hit_data_[ ii ].clear();
		}
	}
}


void TargetRotamerCoords::set_num_target_rotamers(
	Size target_restype_id,
	Size n_rotamers
)
{

	if ( coords_.size() == 1 && target_restype_id == 1 ) {
		n_rots_total_ = n_rotamers;
	}
	if ( n_rotamers != 0 ) {
		//std::cout << "Setting n_rotamers for " << target_restypes_[ ii ]->name() << " " << atom_ids_for_coordinates_[ ii ].size() << " " << n_rotamers_per_target_restype[ ii ] << std::endl;
		coords_[ target_restype_id ].dimension( atom_ids_for_coordinates_[ target_restype_id ].size(), n_rotamers );
		hit_data_[ target_restype_id ].resize( n_rotamers );
	} else {
		coords_[ target_restype_id ].clear();
		hit_data_[ target_restype_id ].clear();
	}

}


void TargetRotamerCoords::set_coordinates_for_rotamer(
	Size restype_index,
	Size rotamer_index,
	Hit const & hit,
	core::conformation::Residue const & rescoords
)
{
	for ( Size ii = 1, ii_end = n_atoms_for_restype( restype_index ); ii <= ii_end; ++ii ) {
		coords_[ restype_index ]( ii, rotamer_index ) = rescoords.xyz( atom_ids_for_coordinates_[ restype_index ][ ii ] );
		/*std::cout << "Coord: " << rescoords.name() << " " << ii << " " << rescoords.atom_name( atom_ids_for_coordinates_[ restype_index ][ ii ] );
		std::cout << " " << coords_[ restype_index ]( ii, rotamer_index ).x();
		std::cout << " " << coords_[ restype_index ]( ii, rotamer_index ).y();
		std::cout << " " << coords_[ restype_index ]( ii, rotamer_index ).z() << std::endl;*/
	}
	hit_data_[ restype_index ][ rotamer_index ] = hit;
}

SecondaryMatchUpstreamResProcessor::SecondaryMatchUpstreamResProcessor(
	SecondaryMatcherToUpstreamResidue & sec_matcher
) :
	sec_matcher_( sec_matcher )
{}

void
SecondaryMatchUpstreamResProcessor::process_hit(
	Hit const & hit,
	core::conformation::Residue const & upstream_conformation
)
{
	sec_matcher_.process_hit( hit, upstream_conformation );
}


}
}
}
