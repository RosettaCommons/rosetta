// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/output/MatchConsolidator.cc
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

// Unit headers
#include <protocols/match/output/MatchConsolidator.hh>

// Package headers
#include <protocols/match/Hit.hh>

#include <protocols/match/output/MatchEvaluator.hh>
#include <protocols/match/output/MatchGrouper.hh>
#include <protocols/match/output/OutputWriter.hh>
#include <protocols/match/output/MatchScoreWriter.hh>

// Utility headers
#include <utility/exit.hh>

// C++ headers
#include <algorithm>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace match {
namespace output {

MatchConsolidator::MatchConsolidator() : n_to_output_per_group_( 5 ) {}

MatchConsolidator::~MatchConsolidator() {}

void
MatchConsolidator::begin_processing()
{
	runtime_assert( grouper_ != 0 );
	runtime_assert( evaluator_ != 0 );
	MatchProcessor::begin_processing();
	match_groups_.clear();
	grouper_->reset();

}

void
MatchConsolidator::process_match(
	match const & m
)
{
	runtime_assert( grouper_ != 0 );
	runtime_assert( evaluator_ != 0 );
	runtime_assert( writer_ != 0 );
	//runtime_assert( match_dspos1_groups_.size() == 0 );
	note_match_processed();

	if ( !this->passes_filters( m ) ) return;

	Size group = grouper_->assign_group_for_match( m );

	if ( group > match_groups_.size() ) {
		match_groups_.resize( group );
		match_groups_[ group ] = protocols::match::output::BestMatchesCollectionOP( new BestMatchesCollection( n_to_output_per_group_, false ) );
	}
	Real score = evaluator_->score( m );
	match_groups_[ group ]->add_match( m, score );
	//std::cout << ".";
}

void
MatchConsolidator::process_match(
	match_dspos1 const & m
)
{
	runtime_assert( grouper_ != 0 );
	runtime_assert( evaluator_ != 0 );
	runtime_assert( writer_ != 0 );
	//runtime_assert( match_groups_.size() == 0 );
	note_match_processed();

	if ( !this->passes_filters( m ) ) return;

	Size group = grouper_->assign_group_for_match( m );

	if ( group > match_groups_.size() ) {
		match_groups_.resize( group );
		match_groups_[ group ] = protocols::match::output::BestMatchesCollectionOP( new BestMatchesCollection( n_to_output_per_group_, true ) );
	}
	Real score = evaluator_->score( m );
	match_groups_[ group ]->add_match_dspos1( m, score );
	//std::cout << ".";
}


void
MatchConsolidator::end_processing()
{
	/// Can only handle regular match objects or match_dspos1's but not both.
	//runtime_assert( match_dspos1_groups_.size() == 0 || match_groups_.size() == 0 );

	runtime_assert( grouper_ != 0 );
	runtime_assert( evaluator_ != 0 );
	runtime_assert( match_score_writer_ != 0 );

	//std::cout << std::endl;
	//std::cout << "match groups: " << match_groups_.size() << std::endl;

	for ( Size ii = 1; ii <= match_groups_.size(); ++ii ) {
		//std::cout << ":";
		for ( Size jj = 1; jj <= match_groups_[ ii ]->n_kept_matches(); ++jj ) {
			if ( ! match_groups_[ ii ]->dspos1_mode() ) {
				writer_->record_match( match_groups_[ ii ]->kept_match( jj ) , evaluator_ , match_score_writer_ );
			} else {
				writer_->record_match( match_groups_[ ii ]->kept_match_dspos1( jj ) );
			}
		}
	}
	match_score_writer_->write_match_scores();
	MatchProcessor::end_processing();
}

void
MatchConsolidator::end_processing_of_regular_match_groups()
{

	/*for ( Size ii = 1; ii <= match_groups_.size(); ++ii ) {
	Size ii_n_matches = match_groups_[ ii ].size();
	utility::vector1< match > ii_matches( ii_n_matches );
	utility::vector1< Real  > ii_scores(  ii_n_matches, 0.0 );
	std::copy( match_groups_[ ii ].begin(), match_groups_[ ii ].end(), ii_matches.begin() );
	for ( Size jj = 1; jj <= ii_n_matches; ++jj ) {
	ii_scores[ jj ] = evaluator_->score( ii_matches[ jj ] );
	}
	utility::vector1< Size > top_score_indices( n_to_output_per_group_ );

	// either sort or, if it's going to be fast, use utility::arg_least_several
	if ( n_to_output_per_group_ > 50 ) {
	utility::vector1< std::pair< Real, Size > > score_index_pairs( ii_n_matches );
	for ( Size jj = 1; jj <= ii_n_matches; ++jj ) {
	score_index_pairs[ jj ] = std::make_pair( ii_scores[ jj ], jj );
	}
	std::sort( score_index_pairs.begin(), score_index_pairs.end(), utility::SortFirst< Real, Size >() );
	for ( Size jj = 1; jj <= n_to_output_per_group_; ++jj ) {
	if ( jj > ii_n_matches ) {
	top_score_indices.resize( ii_n_matches );
	break;
	}
	top_score_indices[ jj ] = score_index_pairs[ jj ].second;
	}
	} else {
	utility::arg_least_several( ii_scores, top_score_indices );
	}

	for ( Size jj = 1; jj <= top_score_indices.size(); ++jj ) {
	writer_->record_match( ii_matches[ top_score_indices[ jj ]] );
	}
	}*/
}

void
MatchConsolidator::end_processing_of_match_dspos1_groups()
{

	/* for ( Size ii = 1; ii <= match_dspos1_groups_.size(); ++ii ) {
	Size ii_n_matches = match_dspos1_groups_[ ii ].size();
	utility::vector1< match_dspos1 > ii_matches( ii_n_matches );
	utility::vector1< Real >         ii_scores(  ii_n_matches, 0.0 );
	std::copy( match_dspos1_groups_[ ii ].begin(), match_dspos1_groups_[ ii ].end(), ii_matches.begin() );
	for ( Size jj = 1; jj <= ii_n_matches; ++jj ) {
	ii_scores[ jj ] = evaluator_->score( ii_matches[ jj ] );
	}
	utility::vector1< Size > top_score_indices( n_to_output_per_group_ );

	// either sort or, if it's going to be fast, use utility::arg_least_several
	if ( n_to_output_per_group_ > 50 ) {
	utility::vector1< std::pair< Real, Size > > score_index_pairs( ii_n_matches );
	for ( Size jj = 1; jj <= ii_n_matches; ++jj ) {
	score_index_pairs[ jj ] = std::make_pair( ii_scores[ jj ], jj );
	}
	std::sort( score_index_pairs.begin(), score_index_pairs.end(), utility::SortFirst< Real, Size >() );
	for ( Size jj = 1; jj <= n_to_output_per_group_; ++jj ) {
	if ( jj > ii_n_matches ) {
	top_score_indices.resize( ii_n_matches );
	break;
	}
	top_score_indices[ jj ] = score_index_pairs[ jj ].second;
	}
	} else {
	utility::arg_least_several( ii_scores, top_score_indices );
	}

	for ( Size jj = 1; jj <= top_score_indices.size(); ++jj ) {
	writer_->record_match( ii_matches[ top_score_indices[ jj ]] );
	}
	}*/
}


void
MatchConsolidator::set_n_to_output_per_group( Size setting )
{
	n_to_output_per_group_ = setting;
}

void
MatchConsolidator::set_grouper( MatchGrouperOP grouper )
{
	grouper_ = grouper;
}

void
MatchConsolidator::reset_grouper()
{
	if ( grouper_ ) grouper_->reset();
}

BestMatchesCollection::~BestMatchesCollection() {}

BestMatchesCollection::BestMatchesCollection(
	Size n_to_keep,
	bool dspos1_mode
) :
	// KAB - below line commented out when removing -Wunused-private-field on 2014-09-11
	// n_top_matches_to_keep_( n_to_keep ),
	dspos1_mode_( dspos1_mode ),
	scores_heap_( n_to_keep ),
	best_matches_( dspos1_mode ? 0 : n_to_keep ),
	best_match_dspos1s_( dspos1_mode ? n_to_keep : 0 )
{
}

void
BestMatchesCollection::add_match( match const & m, Real score )
{
	debug_assert( ! dspos1_mode_ );
	Size new_pos = index_for_new_match( score );
	if ( new_pos != 0 ) {
		best_matches_[ new_pos ] = m;
	}
}

void
BestMatchesCollection::add_match_dspos1( match_dspos1 const & m, Real score )
{
	debug_assert( dspos1_mode_ );
	Size new_pos = index_for_new_match( score );
	if ( new_pos != 0 ) {
		best_match_dspos1s_[ new_pos ] = m;
	}
}

BestMatchesCollection::Size
BestMatchesCollection::n_kept_matches() const
{
	return scores_heap_.size();
}

match const &
BestMatchesCollection::kept_match( Size which_match ) const
{
	debug_assert( ! dspos1_mode_ );
	debug_assert( (int) which_match <= scores_heap_.size() );
	return best_matches_[ which_match ];
}

match_dspos1 const &
BestMatchesCollection::kept_match_dspos1( Size which_match ) const
{
	debug_assert( ! dspos1_mode_ );
	debug_assert( (int) which_match <= scores_heap_.size() );
	return best_match_dspos1s_[ which_match ];
}


BestMatchesCollection::Size
BestMatchesCollection::index_for_new_match( Real score )
{
	Real negscore = -1 * score;
	if ( scores_heap_.size() < scores_heap_.capacity() ) {
		int new_pos( scores_heap_.size() + 1 );
		bool err;
		scores_heap_.heap_insert( new_pos, negscore, err );
		debug_assert( ! err );
		return new_pos;
	} else {
		if ( scores_heap_.heap_head() < negscore ) {
			int old_pos; float oldscore; bool err;
			scores_heap_.heap_extract( old_pos, oldscore, err );
			debug_assert( ! err );
			debug_assert( negscore > oldscore );
			scores_heap_.heap_insert( old_pos, negscore, err );
			debug_assert( ! err );
			return old_pos; // overwrite this position
		} // else return 0 .. moved below to appease the compiler
	}
	return 0;
}

}
}
}
