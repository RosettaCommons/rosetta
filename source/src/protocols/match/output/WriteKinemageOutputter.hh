// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/output/WriteKinemageOutputter.hh
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

#ifndef INCLUDED_protocols_match_output_WriteKinemageOutputter_hh
#define INCLUDED_protocols_match_output_WriteKinemageOutputter_hh

// Unit headers
#include <protocols/match/output/WriteKinemageOutputter.fwd.hh>

// Package Headers
#include <protocols/match/Hit.fwd.hh>
#include <protocols/match/output/OutputWriter.hh>
#include <protocols/match/output/UpstreamHitCacher.fwd.hh>
#include <protocols/match/output/WriteUpstreamCoordinateKineamge.fwd.hh>
#include <protocols/match/output/MatchEvaluator.fwd.hh>
#include <protocols/match/output/MatchScoreWriter.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1_bool.hh>


namespace protocols {
namespace match {
namespace output {

class WriteKinemageOutputter : public OutputWriter
{
public:
	virtual ~WriteKinemageOutputter();

	virtual
	void
	record_match( match const & m , MatchEvaluatorOP evaluator, MatchScoreWriterOP match_score_writer );

	/// @brief evaluator and score writer are not passed in because single-downstream-position match
	///currently have no way of being evaluated
	virtual
	void
	record_match( match_dspos1 const & m );

	void
	set_coordinate_cacher( UpstreamHitCacherOP );

	void
	set_n_geomcst( Size n_geomcst );

	void
	set_kin_writer( WriteUpstreamHitKinemageOP us_writer );

	void
	set_downstream_writer( Size geomcst_id, DownstreamCoordinateKinemageWriterOP ds_writer );

private:

	UpstreamHitCacherOP coordinate_cacher_;
	WriteUpstreamHitKinemageOP upstream_kin_writer_;
	utility::vector1< DownstreamCoordinateKinemageWriterOP > dswriters_;

};

}
}
}

#endif
