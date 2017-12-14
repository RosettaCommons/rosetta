// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/output/MatchOutputter.cc
/// @brief  Implementation of class to write output matches that pass filters
/// This class does not "look ahead" to future matches to decide whether the current match
/// to process should be output, however, filters are able to keep a history of what they've
/// output so far.
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

// Unit headers
#include <protocols/match/output/MatchOutputter.hh>

// Package headers
#include <protocols/match/output/MatchProcessor.hh>
#include <protocols/match/output/OutputWriter.hh>
#include <protocols/match/output/MatchEvaluator.hh>
#include <protocols/match/output/MatchScoreWriter.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/exit.hh>

// C++ headers

namespace protocols {
namespace match {
namespace output {

MatchOutputter::MatchOutputter()
: MatchProcessor()
{}

MatchOutputter::~MatchOutputter() = default;

void
MatchOutputter::begin_processing()
{
	MatchProcessor::begin_processing();
}

void
MatchOutputter::end_processing()
{
	match_score_writer_->write_match_scores();
	MatchProcessor::end_processing();
}

void
MatchOutputter::process_match(
	match const & m
)
{
	note_match_processed();
	if ( !this->passes_filters( m ) ) return;

	if ( writer_ ) {
		runtime_assert( evaluator_ != nullptr );
		runtime_assert( match_score_writer_ != nullptr );
		writer_->record_match( m , evaluator_ , match_score_writer_ );
	}
}


void
MatchOutputter::process_match(
	match_dspos1 const & m
)
{
	note_match_processed();
	if ( !this->passes_filters( m ) ) return;

	if ( writer_ ) {
		writer_->record_match( m );
	}
}

}
}
}
