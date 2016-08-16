// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/output/ProcessorFactory.hh
/// @brief  Declaration for a factory class responsible for instantiating
///         MatchProcessor classes.
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

#ifndef INCLUDED_protocols_match_output_ProcessorFactory_hh
#define INCLUDED_protocols_match_output_ProcessorFactory_hh

// Unit headers
#include <protocols/match/output/ProcessorFactory.fwd.hh>

// Package headers
#include <protocols/match/Matcher.fwd.hh>
#include <protocols/match/MatcherTask.fwd.hh>
#include <protocols/match/output/MatchProcessor.fwd.hh>
#include <protocols/match/output/MatchEvaluator.fwd.hh>
#include <protocols/match/output/MatchFilter.fwd.hh>
#include <protocols/match/output/MatchGrouper.fwd.hh>
#include <protocols/match/output/OutputWriter.fwd.hh>
#include <protocols/match/output/UpstreamHitCacher.fwd.hh>
#include <protocols/match/output/MatchScoreWriter.fwd.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <list>

namespace protocols {
namespace match {
namespace output {

class ProcessorFactory : public utility::pointer::ReferenceCount {
public:
	typedef core::Real Real;
	typedef core::Size Size;
public:
	ProcessorFactory();

	virtual
	~ProcessorFactory();

	static
	MatchProcessorOP
	create_processor(
		MatcherCOP matcher,
		MatcherTaskCOP task
	);

private:
	static
	MatchEvaluatorOP
	create_evaluator(
		MatcherCOP matcher,
		MatcherTaskCOP task,
		UpstreamHitCacherOP cacher
	);

	static
	std::list< MatchFilterOP >
	create_filters(
		MatcherCOP matcher,
		MatcherTaskCOP task,
		UpstreamHitCacherOP cacher
	);

	static
	MatchGrouperOP
	create_grouper(
		MatcherCOP matcher,
		MatcherTaskCOP task,
		UpstreamHitCacherOP cacher
	);


	static
	OutputWriterOP
	create_output_writer(
		MatcherCOP matcher,
		MatcherTaskCOP task,
		UpstreamHitCacherOP cacher
	);

	static
	MatchScoreWriterOP
	create_match_score_writer();

};

}
}
}

#endif
