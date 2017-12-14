// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/output/UpstreamHitCacher.hh
/// @brief  Declaration for class to cache and recall the conformations of upstream hits.
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

// Unit headers
#include <protocols/match/output/UpstreamHitCacher.hh>

// Package headers
#include <protocols/match/Hit.hh>
#include <protocols/match/Matcher.hh>
#include <protocols/match/upstream/ScaffoldBuildPoint.hh>

// Project headers
#include <core/conformation/Residue.hh>

// Utility headers
#include <utility>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/exit.hh>

// C++ headers
#include <map>

#include <utility/OrderedTuple.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace match {
namespace output {


UpstreamHitCacher::UpstreamHitCacher( MatcherCOP matcher ) :
	matcher_(std::move( matcher )),
	n_geometric_constraints_( matcher_ ? matcher_->n_geometric_constraints() : 0 ),
	n_confs_to_cache_( 100 ),
	index_for_rotamer_( n_geometric_constraints_ ),
	which_cst_being_processed_( 0 ),
	queue_head_( n_geometric_constraints_, 0 ),
	scafrot_pair_for_conf_( n_geometric_constraints_ ),
	upstream_confs_( n_geometric_constraints_ )
{
	resize_arrays();
}

UpstreamHitCacher::~UpstreamHitCacher() = default;

void
UpstreamHitCacher::set_cache_size( Size n_rotamers_to_cache )
{
	n_confs_to_cache_ = n_rotamers_to_cache;
	resize_arrays();
}


core::conformation::ResidueCOP
UpstreamHitCacher::upstream_conformation_for_hit( Size cst_id, Hit const & hit )
{
	core::conformation::ResidueCOP residue;
#ifdef USE_OPENMP
	#pragma omp critical ( upstream_hit_cacher_upstream_conformation_for_hit )
#endif
	{
		ScaffoldRotamerPair srp; srp[ 1 ] = hit.scaffold_build_id(); srp[ 2 ] = hit.upstream_conf_id();
		ScaffoldRotamerTuple srt( srp );
		Size index = already_in_queue( cst_id, srt );
		if ( index == 0 ) {
			index = fetch( cst_id, srt );
		}
		residue = upstream_confs_[ cst_id ][ index ];
	}
	return residue;
}


void
UpstreamHitCacher::process_hit(
	Hit const & hit,
	core::conformation::Residue const & upstream_conformation
)
{
	runtime_assert( which_cst_being_processed_ != 0 );

	Size const cst_id = which_cst_being_processed_;

	if ( queue_head_[ cst_id ] == 0 ) {
		/// empty queue condition
		++queue_head_[ cst_id ];
	} else {
		Size next = queue_head_[ cst_id ] + 1;
		if ( next == n_confs_to_cache_ + 1 ) next = 1;

		if ( upstream_confs_[ cst_id ][ next ] ) {
			// erase the entry for this rotamer in the index_for_rotamer_ map;
			auto iter =
				index_for_rotamer_[ cst_id ].find( scafrot_pair_for_conf_[ cst_id ][ next ] );
			runtime_assert( iter != index_for_rotamer_[ cst_id ].end() );
			index_for_rotamer_[ cst_id ].erase( iter );
		}
		queue_head_[ cst_id ] = next;
	}

	Size index = queue_head_[ cst_id ];

	upstream_confs_[ cst_id ][ index ] = core::conformation::ResidueCOP( core::conformation::ResidueOP( new core::conformation::Residue( upstream_conformation ) ) );

	ScaffoldRotamerPair srp; srp[ 1 ] = hit.scaffold_build_id(); srp[ 2 ] = hit.upstream_conf_id();
	ScaffoldRotamerTuple srt( srp );

	scafrot_pair_for_conf_[ cst_id ][ index ] =  srt;
	index_for_rotamer_[ cst_id ][ srt ] = index;

}


void
UpstreamHitCacher::resize_arrays()
{
	debug_assert( index_for_rotamer_.size() == n_geometric_constraints_ );
	debug_assert( scafrot_pair_for_conf_.size() == n_geometric_constraints_ );
	debug_assert( upstream_confs_.size() == n_geometric_constraints_ );
	debug_assert( queue_head_.size() == n_geometric_constraints_ );

	std::fill( queue_head_.begin(), queue_head_.end(), 0 );
	for ( Size ii = 1; ii <= n_geometric_constraints_; ++ii ) {
		index_for_rotamer_[ ii ].clear();
		scafrot_pair_for_conf_[ ii ].clear();
		scafrot_pair_for_conf_[ ii ].resize( n_confs_to_cache_ );
		upstream_confs_[ ii ].clear();
		upstream_confs_[ ii ].resize( n_confs_to_cache_ );
	}
}


/// @brief Returns 0 if the scaffold/rotamer pair is not already in the queue, and
/// the non-zero index in the queue if it is.
UpstreamHitCacher::Size
UpstreamHitCacher::already_in_queue( Size cst_id, ScaffoldRotamerTuple const & rotid ) const
{
	auto iter = index_for_rotamer_[ cst_id ].find( rotid );
	if ( iter == index_for_rotamer_[ cst_id ].end() ) {
		return 0;
	} else {
		return iter->second;
	}
}

/// @brief Construct the rotamer for the requested scaffold/rotamer pair and
/// put it into the queue, evicting the previous queue resident if necessary.
UpstreamHitCacher::Size
UpstreamHitCacher::fetch( Size cst_id, ScaffoldRotamerTuple const & rotid )
{
	which_cst_being_processed_ = cst_id;

	Hit hit;
	hit.first()[ 1 ] = rotid.data()[ 1 ];
	hit.first()[ 2 ] = rotid.data()[ 2 ];

	matcher_->upstream_builder( cst_id )->recover_hit(
		hit,
		* matcher_->build_point( hit.scaffold_build_id() ),
		* this
	);

	which_cst_being_processed_ = 0;

	return queue_head_[ cst_id ];
}

/*
private:
MatcherCOP matcher_;

Size n_geometric_constraints_;
Size n_confs_to_cache_;

utility::vector1< std::map< ScaffRotamerTuple, Size > > index_for_rotamer_;

utility::vector1< Size > queue_head_;
utility::vector1< utility::vector1< ScaffoldRotamerTuple > >           scafrot_pair_for_conf_;
utility::vector1< utility::vector1< core::conformation::ResidueCOP > > upstream_confs_;

};
*/

}
}
}

