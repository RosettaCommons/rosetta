// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/upstream/UpstreamBuilder.cc
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

// Unit header
#include <protocols/match/upstream/UpstreamBuilder.hh>

// Package header
#include <protocols/match/BumpGrid.hh>
#include <protocols/match/upstream/ProteinUpstreamBuilder.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace match {
namespace upstream {

UpstreamBuilder::~UpstreamBuilder() {}

bool UpstreamBuilder::compatible(
	Hit const & my_hit,
	ScaffoldBuildPoint const & build_point_mine,
	UpstreamBuilder const & other,
	Hit const & other_hit,
	ScaffoldBuildPoint const & build_point_other,
	bool first_dispatch
) const
{
	if ( ! first_dispatch ) {
		utility_exit_with_message( "CRITICAL ERROR: UpstreamBuilder::compatible() double-dispatch failure" );
	}
	return other.compatible( other_hit, build_point_other, *this, my_hit, build_point_mine, false );

}

bool UpstreamBuilder::compatible(
	Hit const & my_hit,
	ScaffoldBuildPoint const & build_point_mine,
	ProteinUpstreamBuilder const & other,
	Hit const & other_hit,
	ScaffoldBuildPoint const & build_point_other,
	bool first_dispatch
) const
{
	if ( ! first_dispatch ) {
		utility_exit_with_message( "CRITICAL ERROR: UpstreamBuilder::compatible() double-dispatch failure" );
	}
	return other.compatible( other_hit, build_point_other, *this, my_hit, build_point_mine, false );

}


void
UpstreamBuilder::set_bb_grid( BumpGridCOP bbgrid )
{
	bbgrid_ = bbgrid;
}

UpstreamResidueProcessor::~UpstreamResidueProcessor() {}


}
}
}
