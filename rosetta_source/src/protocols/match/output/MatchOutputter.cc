// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <iostream>

namespace protocols {
namespace match {
namespace output {

MatchOutputter::MatchOutputter()
	: MatchProcessor()
{}

MatchOutputter::~MatchOutputter() {}

void
MatchOutputter::begin_processing()
{
	MatchProcessor::begin_processing();
}

void
MatchOutputter::end_processing()
{
	MatchProcessor::end_processing();
}

void
MatchOutputter::process_match(
	match const & m
)
{
	note_match_processed();
	if( !this->passes_filters( m ) ) return;

	if ( writer_ ) {
		writer_->record_match( m );
	}
}


void
MatchOutputter::process_match(
	match_dspos1 const & m
)
{
	note_match_processed();
	if( !this->passes_filters( m ) ) return;

	if ( writer_ ) {
		writer_->record_match( m );
	}
}

}
}
}
