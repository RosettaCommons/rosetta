// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/output/MatchProcessor.hh
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), Florian Richter (floric@u.washington.edu), porting to mini

#ifndef INCLUDED_protocols_match_output_MatchProcessor_hh
#define INCLUDED_protocols_match_output_MatchProcessor_hh

// Unit headers
#include <protocols/match/output/MatchProcessor.fwd.hh>

// Package headers
#include <protocols/match/Hit.fwd.hh>
#include <protocols/match/output/MatchFilter.fwd.hh>
#include <protocols/match/output/OutputWriter.fwd.hh>
#include <protocols/match/output/UpstreamCollisionFilter.fwd.hh>
#include <protocols/match/output/UpstreamDownstreamCollisionFilter.fwd.hh>
#include <protocols/match/output/MatchEvaluator.fwd.hh>
#include <protocols/match/output/MatchScoreWriter.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// c++ headers
#include <map>

#ifdef WIN32
#include <protocols/match/Hit.hh>
#endif


namespace protocols {
namespace match {
namespace output {

class MatchProcessor : public utility::pointer::ReferenceCount {
public:

	MatchProcessor();

	virtual
	~MatchProcessor();

	/// @brief Invoked by the Matcher before it begins feeding matches to the processor.
	/// Base-class has a no-op implementation.
	virtual
	void
	begin_processing();

	virtual
	void
	process_match(
		match const & m
	) = 0;

	virtual
	void
	process_match(
		match_dspos1 const & m
	) = 0;

	/// @brief Invoked by the Matcher after it finishes feeding matches to the processor.
	/// Base-class has a no-op implementation.
	virtual
	void
	end_processing();

	bool
	match_processing_successful() const {
		return match_processing_successful_; }

	void
	set_output_writer( OutputWriterOP writer );

	void
	set_evaluator( MatchEvaluatorOP evaluator );

	void
	set_match_score_writer( MatchScoreWriterOP scorewriter );

	/// @brief const access to the output writer, e.g.
	/// for the MatcherMover when incorporating matches
	/// into the pose
	OutputWriterOP
	output_writer();

	void
	add_filter( MatchFilterOP filter );

	void
	reset_filters();

	void
	clear_filters();

	UpstreamDownstreamCollisionFilterCOP
	up_down_filt() const;

	UpstreamCollisionFilterCOP
	up_coll_filt() const;


protected:

	bool
	passes_filters(
		match const & m
	);

	bool
	passes_filters(
		match_dspos1 const & m
	);

	void
	note_match_processed();

	bool match_processing_successful_;

	OutputWriterOP writer_;

	MatchEvaluatorOP evaluator_;
	MatchScoreWriterOP match_score_writer_;

private:

	void
	note_filter_fail( std::string filter_name );

	std::list< MatchFilterOP > filters_;/// NOTE: these filters are not given the opportunity to accumulate state.
	std::list< StateAccumulatingMatchFilterOP > filters_with_state_;
	std::map< std::string, core::Size > filter_failcounts_;
	core::Size num_matches_processed_;

	/// NOTE: this is used for clash checking between upstream only hits and the ligand
	UpstreamDownstreamCollisionFilterCOP up_down_filt_;

	UpstreamCollisionFilterCOP up_coll_filt_;

};

}
}
}

#endif
