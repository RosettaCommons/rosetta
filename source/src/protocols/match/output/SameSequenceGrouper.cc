// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/output/SameSequenceGrouper.cc
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

// Unit headers
#include <protocols/match/output/SameSequenceGrouper.hh>

// Package headers
#include <protocols/match/Hit.hh>
#include <protocols/match/downstream/DownstreamBuilder.hh>
#include <protocols/match/output/UpstreamHitCacher.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>

// Utility headers
#include <utility/exit.hh>

#include <utility/OrderedTuple.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace match {
namespace output {

SameSequenceGrouper::SameSequenceGrouper() : n_geometric_constraints_( 0 ) {}
SameSequenceGrouper::SameSequenceGrouper( Size ncst ) : n_geometric_constraints_( ncst ) {}

SameSequenceGrouper::~SameSequenceGrouper() {}

SameSequenceGrouper::Size
SameSequenceGrouper::assign_group_for_match(
	match const & m
)
{
	return assign_group_for_match( match_dspos1( m, 1 ) );
}

SameSequenceGrouper::Size
SameSequenceGrouper::assign_group_for_match(
	match_dspos1 const & m
)
{
	runtime_assert( m.upstream_hits.size() == n_geometric_constraints_ );

	utility::vector1< Size > seq_vector( n_geometric_constraints_ * 2, 0 );
	for ( Size ii = 1; ii <= n_geometric_constraints_; ++ii ) {
		seq_vector[ ii ] = m.upstream_hits[ ii ].scaffold_build_id();
		seq_vector[ n_geometric_constraints_ + ii ] =
			hit_cacher_->upstream_conformation_for_hit( ii, fake_hit( m.upstream_hits[ ii ] ) )->aa();
	}

	SequenceMap::const_iterator iter = sequence_indexer_.find( seq_vector );
	if ( iter == sequence_indexer_.end() ) {
		Size next_index = sequence_indexer_.size() + 1;
		sequence_indexer_[ seq_vector ] = next_index;
		return next_index;
	} else {
		return iter->second;
	}
}

void
SameSequenceGrouper::reset()
{
	sequence_indexer_.clear();
}

void
SameSequenceGrouper::set_n_geometric_constraints( Size n_csts )
{
	n_geometric_constraints_ = n_csts;
}

void
SameSequenceGrouper::set_hit_cacher( UpstreamHitCacherOP cacher )
{
	hit_cacher_ = cacher;
}


SameSequenceAndDSPositionGrouper::SameSequenceAndDSPositionGrouper() :
	SameSequenceGrouper(), rms_group_cutoff_(1.0) {}
SameSequenceAndDSPositionGrouper::SameSequenceAndDSPositionGrouper( Size ncst ) :
	SameSequenceGrouper( ncst ), rms_group_cutoff_(1.0)
{
	dsbuilders_.resize( ncst, NULL );
}

SameSequenceAndDSPositionGrouper::~SameSequenceAndDSPositionGrouper() {}

void
SameSequenceAndDSPositionGrouper::set_n_geometric_constraints( Size n_csts )
{
	parent::set_n_geometric_constraints( n_csts );
	dsbuilders_.resize( n_csts );
}

void
SameSequenceAndDSPositionGrouper::set_rms_group_cutoff( Real cutoff )
{
	rms_group_cutoff_ = cutoff;
}


SameSequenceGrouper::Size
SameSequenceAndDSPositionGrouper::assign_group_for_match(
	match_dspos1 const & m
){
	//std::cout << "seq ds assigning group for match " << std::endl;
	Size sequence_group( parent::assign_group_for_match( m ) );
	Size ds_group( assign_downstream_position_group_for_match( m ) );
	std::pair< Size, Size > seq_dspos_pair( sequence_group, ds_group );

	std::map< std::pair< Size, Size >, Size >::iterator seqpos_it = sequence_pos_map_.find( seq_dspos_pair );
	if( seqpos_it == sequence_pos_map_.end() ){
		Size newgroup( sequence_pos_map_.size() + 1 );
		sequence_pos_map_.insert( std::make_pair( seq_dspos_pair, newgroup ) );
		return newgroup;
	}
	return seqpos_it->second;
}


void
SameSequenceAndDSPositionGrouper::reset()
{
	parent::reset();
	sequence_pos_map_.clear();
	representative_dspos_.clear();
}

void
SameSequenceAndDSPositionGrouper::set_relevant_atom_ids(
	utility::vector1< core::id::AtomID > const & relevant_atom_ids
)
{
	relevant_atom_ids_ = relevant_atom_ids;
}

void
SameSequenceAndDSPositionGrouper::set_downstream_builder(
	Size geomcst_id,
	downstream::DownstreamBuilderCOP dsbuilder
)
{
	runtime_assert( dsbuilders_.size() >= geomcst_id );
	dsbuilders_[geomcst_id] = dsbuilder;
}


SameSequenceGrouper::Size
SameSequenceAndDSPositionGrouper::assign_downstream_position_group_for_match(
	match_dspos1 const & m
){
	utility::vector1< Vector > dspos_coords( relevant_atom_ids_.size() );
	Hit fullhit = full_hit( m );
	dsbuilders_[ m.originating_geom_cst_for_dspos ]->coordinates_from_hit( fullhit, relevant_atom_ids_, dspos_coords );
	runtime_assert( dspos_coords.size() == relevant_atom_ids_.size() );
	Real lowest_rms( 1000000.0 );

	for( core::Size ii = 1; ii <= representative_dspos_.size(); ++ii ){
		//calculate rms without superposition
		Real cur_rms(0.0);
		for( core::Size jj = 1; jj <= relevant_atom_ids_.size(); ++jj ){
			Vector diff = dspos_coords[jj] - representative_dspos_[ii][jj];
			cur_rms += diff.length_squared();
		}
		Real rms = std::sqrt(cur_rms / relevant_atom_ids_.size() );
		if( rms < lowest_rms ) lowest_rms = rms;
		if( rms  <= rms_group_cutoff_ ) return ii;
	}
	//if we've made it till here, that means this ds position is new
	//std::cout << "found new group with lowest rms of " << lowest_rms << " to previous one." << std::endl;
	representative_dspos_.push_back( dspos_coords );
	return representative_dspos_.size();
}

} //namespace output
} //namespace match
} //namespace protocols
