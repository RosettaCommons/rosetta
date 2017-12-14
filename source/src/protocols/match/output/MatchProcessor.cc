// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/output/MatchProcessor.cc
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini
/// @author Florian Richter( floric@u.washington.edu)

// Unit headers
#include <protocols/match/output/MatchProcessor.hh>

// package headers
#include <protocols/match/output/MatchFilter.hh>
#include <protocols/match/output/UpstreamDownstreamCollisionFilter.hh>
#include <protocols/match/output/OutputWriter.hh>
#include <protocols/match/output/MatchEvaluator.hh>
#include <protocols/match/output/MatchScoreWriter.hh>

//Project headers
#include <basic/Tracer.hh>


// Utility headers
#include <utility/exit.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace match {
namespace output {

static basic::Tracer TR( "protocols.match.output.MatchProcessor" );

MatchProcessor::MatchProcessor()
: match_processing_successful_(false), writer_(/* NULL */),
	num_matches_processed_(0), up_down_filt_(/* NULL */), up_coll_filt_(nullptr)
{
	filter_failcounts_.clear();
}


MatchProcessor::~MatchProcessor() = default;

void
MatchProcessor::begin_processing()
{
	runtime_assert( writer_ != nullptr );
	num_matches_processed_ = 0;
	match_processing_successful_ = false;
	filter_failcounts_.clear();
	writer_->prepare_for_output_writing();
}

void
MatchProcessor::end_processing()
{

	TR << "A total of " << num_matches_processed_ << " matches were processed." << std::endl;
	core::Size total_filter_fails(0);
	if ( filter_failcounts_.size() != 0 ) {
		for ( std::map< std::string, core::Size >::const_iterator map_it( filter_failcounts_.begin() ), map_end( filter_failcounts_.end() ); map_it != map_end; ++map_it ) {
			TR << map_it->second << " matches did not pass filter " << map_it->first << ".\n" ;
			total_filter_fails += map_it->second;
		}
		if ( num_matches_processed_ == total_filter_fails ) TR << "Apparently no matches passed all filters." << std::endl;
	} else TR << "No matches failed any filter.\n";

	if ( num_matches_processed_ > total_filter_fails ) {
		match_processing_successful_ = true;
	}
	TR.flush();
	writer_->end_output_writing();
}

void
MatchProcessor::set_output_writer( OutputWriterOP writer )
{
	writer_ = writer;
}

OutputWriterOP
MatchProcessor::output_writer()
{
	return writer_;
}

void
MatchProcessor::note_filter_fail( std::string filter_name )
{
	auto map_it( filter_failcounts_.find( filter_name ) );
	if ( map_it == filter_failcounts_.end() ) {
		filter_failcounts_.insert( std::pair< std::string, core::Size >( filter_name, 1 ) );
	} else map_it->second++;
}

void
MatchProcessor::note_match_processed()
{
	num_matches_processed_++;
}

bool
MatchProcessor::passes_filters(
	match const & m
)
{

	for ( std::list< MatchFilterOP >::const_iterator
			iter = filters_.begin(), iter_end = filters_.end();
			iter != iter_end; ++iter ) {
		if (  ! (*iter)->passes_filter( m ) ) {
			note_filter_fail( (*iter)->filter_name() );
			return false;
		}
	}

	for ( std::list< StateAccumulatingMatchFilterOP >::const_iterator
			iter = filters_with_state_.begin(), iter_end = filters_with_state_.end();
			iter != iter_end; ++iter ) {
		(*iter)->note_match_accepted( m );
	}
	return true;
}

bool
MatchProcessor::passes_filters(
	match_dspos1 const & m
)
{
	for ( std::list< MatchFilterOP >::const_iterator
			iter = filters_.begin(), iter_end = filters_.end();
			iter != iter_end; ++iter ) {
		if (  ! (*iter)->passes_filter( m ) ) {
			note_filter_fail( (*iter)->filter_name() );
			return false;
		}
	}

	for ( std::list< StateAccumulatingMatchFilterOP >::const_iterator
			iter = filters_with_state_.begin(), iter_end = filters_with_state_.end();
			iter != iter_end; ++iter ) {
		(*iter)->note_match_accepted( m );
	}
	return true;
}


void
MatchProcessor::add_filter( MatchFilterOP filter )
{

	UpstreamCollisionFilterCOP up_coll_filt(
		utility::pointer::dynamic_pointer_cast< UpstreamCollisionFilter > ( filter ) );
	if ( up_coll_filt ) {
		up_coll_filt_ = up_coll_filt;
		return;  //upstream collision filter not happening at match enumeration stage anymore
	}

	filters_.push_back( filter );

	/// Is this a state-accululating filter?
	StateAccumulatingMatchFilterOP filter_with_state(
		utility::pointer::dynamic_pointer_cast< StateAccumulatingMatchFilter > ( filter ) );
	if ( filter_with_state ) {
		filters_with_state_.push_back( filter_with_state );
	}

	UpstreamDownstreamCollisionFilterCOP up_down_filt(
		utility::pointer::dynamic_pointer_cast< UpstreamDownstreamCollisionFilter > ( filter ) );
	if ( up_down_filt ) up_down_filt_ = up_down_filt;

}


void
MatchProcessor::reset_filters()
{
	for ( std::list< StateAccumulatingMatchFilterOP >::const_iterator
			iter = filters_with_state_.begin(), iter_end = filters_with_state_.end();
			iter != iter_end; ++iter ) {
		(*iter)->reset();
	}
}

void
MatchProcessor::clear_filters()
{
	filters_.clear();
}

UpstreamDownstreamCollisionFilterCOP
MatchProcessor::up_down_filt() const
{
	return up_down_filt_;
}

UpstreamCollisionFilterCOP
MatchProcessor::up_coll_filt() const
{
	return up_coll_filt_;
}

void
MatchProcessor::set_evaluator( MatchEvaluatorOP evaluator )
{
	evaluator_ = evaluator;
}

void
MatchProcessor::set_match_score_writer( MatchScoreWriterOP scorewriter)
{
	match_score_writer_ = scorewriter;
}

}
}
}
