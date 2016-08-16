// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/downstream/DownstreamBuilder.cc
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini
/// @author Florian Richter (floric@u.washington.edu)

/// Unit headers
#include <protocols/match/downstream/DownstreamBuilder.hh>

/// Package headers
#include <protocols/match/BumpGrid.hh>
#include <protocols/match/OccupiedSpaceHash.hh>
#include <protocols/match/downstream/ActiveSiteGrid.hh>
#include <protocols/match/downstream/RigidLigandBuilder.hh>
#include <protocols/match/downstream/LigandConformerBuilder.hh>

//Project header
#include <basic/Tracer.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace match {
namespace downstream {

static THREAD_LOCAL basic::Tracer TR( "protocols.match.downstream.SecondaryMatcherToDownstreamResidue" );

DownstreamBuilder::DownstreamBuilder() : parent(), bb_grid_( /* 0 */ ) {}

DownstreamBuilder::DownstreamBuilder( DownstreamBuilder const & other ) :
	parent(),
	bb_grid_( other.bb_grid_ )
{}

DownstreamBuilder::~DownstreamBuilder()
{}


void
DownstreamBuilder::set_bb_grid(
	BumpGridCOP bbgrid
)
{
	bb_grid_ = bbgrid;
}

void
DownstreamBuilder::set_occupied_space_hash(
	OccupiedSpaceHashCOP occ_space
)
{
	space_ = occ_space;
}

void
DownstreamBuilder::set_active_site_grid(
	ActiveSiteGridCOP active_site_grid
)
{
	active_site_grid_ = active_site_grid;
}

bool
DownstreamBuilder::compatible(
	Hit const & my_hit,
	DownstreamBuilder const & other,
	Hit const & other_hit,
	bool first_dispatch
) const
{
	if ( ! first_dispatch ) {
		utility_exit_with_message( "CRITICAL ERROR: DownstreamBuilder::compatible() double-dispatch failure" );
	}
	return other.compatible( other_hit, *this, my_hit, false );
}

bool
DownstreamBuilder::compatible(
	Hit const & my_hit,
	RigidLigandBuilder const & other,
	Hit const & other_hit,
	bool first_dispatch
) const
{
	if ( ! first_dispatch ) {
		utility_exit_with_message( "CRITICAL ERROR: DownstreamBuilder::compatible() double-dispatch failure" );
	}
	return other.compatible( other_hit, *this, my_hit, false );
}

bool
DownstreamBuilder::compatible(
	Hit const & my_hit,
	LigandConformerBuilder const & other,
	Hit const & other_hit,
	bool first_dispatch
) const
{
	if ( ! first_dispatch ) {
		utility_exit_with_message( "CRITICAL ERROR: DownstreamBuilder::compatible() double-dispatch failure" );
	}
	return other.compatible( other_hit, *this, my_hit, false );
}

}
}
}
