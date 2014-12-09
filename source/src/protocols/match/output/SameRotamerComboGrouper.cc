// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/output/SameRotamerComboGrouper.cc
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

// Unit headers
#include <protocols/match/output/SameRotamerComboGrouper.hh>

// Package headers
#include <protocols/match/output/MatchGrouper.hh>
#include <protocols/match/Hit.hh>

// Utility headers
#include <utility/exit.hh>

#include <utility/OrderedTuple.hh>


namespace protocols {
namespace match {
namespace output {

SameRotamerComboGrouper::SameRotamerComboGrouper() : n_geometric_constraints_( 0 ) {}
SameRotamerComboGrouper::SameRotamerComboGrouper( Size ncst ) : n_geometric_constraints_( ncst ) {}

SameRotamerComboGrouper::~SameRotamerComboGrouper() {}

SameRotamerComboGrouper::Size
SameRotamerComboGrouper::assign_group_for_match(
	match const & m
)
{
	return assign_group_for_match( match_dspos1( m, 1 ) );
}

SameRotamerComboGrouper::Size
SameRotamerComboGrouper::assign_group_for_match(
	match_dspos1 const & m
)
{
	runtime_assert( m.upstream_hits.size() == n_geometric_constraints_ );

	utility::vector1< Size > rot_vector( n_geometric_constraints_ * 2, 0 );
	for ( Size ii = 1; ii <= n_geometric_constraints_; ++ii ) {
		rot_vector[ ii ] = m.upstream_hits[ ii ].scaffold_build_id();
		rot_vector[ n_geometric_constraints_ + ii ] = m.upstream_hits[ ii ].upstream_conf_id();
	}

	RotamerComboCountMap::const_iterator iter = rotamer_combo_indexer_.find( rot_vector );
	if ( iter == rotamer_combo_indexer_.end() ) {
		Size next_index = rotamer_combo_indexer_.size() + 1;
		rotamer_combo_indexer_[ rot_vector ] = next_index;
		return next_index;
	} else {
		return iter->second;
	}
}

void
SameRotamerComboGrouper::reset()
{
	rotamer_combo_indexer_.clear();
}

void
SameRotamerComboGrouper::set_n_geometric_constraints( Size n_csts )
{
	n_geometric_constraints_ = n_csts;
}


}
}
}
