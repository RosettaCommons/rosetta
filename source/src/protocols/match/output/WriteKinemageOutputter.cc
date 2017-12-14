// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/output/WriteKinemageOutputter.cc
/// @brief  Forward declaration of class to write output matches that have (presumably) passed the output filters.
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

// Unit headers
#include <protocols/match/output/WriteKinemageOutputter.hh>

// Package Headers
#include <protocols/match/Hit.hh>
#include <protocols/match/output/UpstreamHitCacher.hh>
#include <protocols/match/output/WriteUpstreamCoordinateKineamge.hh>
#include <protocols/match/output/MatchEvaluator.hh>
#include <protocols/match/output/MatchScoreWriter.hh>

// Project Headers
#include <core/conformation/Residue.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace match {
namespace output {

WriteKinemageOutputter::~WriteKinemageOutputter() = default;


void
WriteKinemageOutputter::record_match( match const & m , MatchEvaluatorOP evaluator, MatchScoreWriterOP match_score_writer )
{
	upstream_kin_writer_->start_new_match();

	for ( Size ii = 1; ii <= m.size(); ++ii ) {
		core::conformation::ResidueCOP conf =
			coordinate_cacher_->upstream_conformation_for_hit( ii, m[ ii ] );

		upstream_kin_writer_->set_dswriter( dswriters_[ ii ] );
		upstream_kin_writer_->geom_id( ii );
		upstream_kin_writer_->output_hit( m[ ii ], *conf );
		match_score_writer->add_match( evaluator->score(m) );
	}
}

void
WriteKinemageOutputter::record_match( match_dspos1 const & m )
{
	upstream_kin_writer_->start_new_match();
	Size const ds_hit_id = m.originating_geom_cst_for_dspos;
	upstream_kin_writer_->set_dswriter( dswriters_[ ds_hit_id ] );

	for ( Size ii = 1; ii <= m.upstream_hits.size(); ++ii ) {
		core::conformation::ResidueCOP conf =
			coordinate_cacher_->upstream_conformation_for_hit( ii, fake_hit( m.upstream_hits[ ii ] ) );

		upstream_kin_writer_->geom_id( ii );
		if ( ii == ds_hit_id ) {
			upstream_kin_writer_->output_hit( full_hit( m ), *conf );
		} else {
			upstream_kin_writer_->output_upstream_coordinates( m.upstream_hits[ ii ], *conf );
		}
	}
}

void
WriteKinemageOutputter::set_coordinate_cacher( UpstreamHitCacherOP cacher )
{
	coordinate_cacher_ = cacher;
}

void
WriteKinemageOutputter::set_n_geomcst( Size n_geomcst )
{
	debug_assert( dswriters_.size() == 0 ); // only set once!
	dswriters_.resize( n_geomcst );
}

void
WriteKinemageOutputter::set_kin_writer(
	WriteUpstreamHitKinemageOP writer
)
{
	upstream_kin_writer_ = writer;
}

void
WriteKinemageOutputter::set_downstream_writer(
	Size geomcst_id,
	DownstreamCoordinateKinemageWriterOP ds_writer
)
{
	dswriters_[ geomcst_id ] = ds_writer;
}

}
}
}
