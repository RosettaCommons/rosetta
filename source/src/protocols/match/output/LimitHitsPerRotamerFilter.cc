// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/output/LimitHitsPerRotamerFilter.hh
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini


// Unit headers
#include <protocols/match/output/LimitHitsPerRotamerFilter.hh>

// Package headers
#include <protocols/match/output/MatchFilter.hh>
#include <protocols/match/Hit.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/OrderedTuple.hh>
// AUTO-REMOVED #include <utility/vector1.hh>

// C++ headers
#include <map>

namespace protocols {
namespace match {
namespace output {

LimitHitsPerRotamerFilter::LimitHitsPerRotamerFilter()
	: StateAccumulatingMatchFilter("LimitHitsPerRotamerFilter")
{}
LimitHitsPerRotamerFilter::LimitHitsPerRotamerFilter( Size n_geometric_constraints ) :
	StateAccumulatingMatchFilter("LimitHitsPerRotamerFilter"),
	n_geometric_constraints_( n_geometric_constraints )
{}

void
LimitHitsPerRotamerFilter::set_n_geometric_constraints( Size n_csts )
{
	n_geometric_constraints_ = n_csts;

}

void
LimitHitsPerRotamerFilter::set_limit_for_rotamer_combo( Size limit )
{
	limit_per_rotamer_combo_ = limit;
}

LimitHitsPerRotamerFilter::~LimitHitsPerRotamerFilter() {}

/// @brief Returns true if the given match passes this filter
bool
LimitHitsPerRotamerFilter::passes_filter(
	match const & m
) const
{
	runtime_assert( m.size() == n_geometric_constraints_ );

	utility::vector1< Size > rot_vector( n_geometric_constraints_ * 2, 0 );
	for ( Size ii = 1; ii <= n_geometric_constraints_; ++ii ) {
		rot_vector[ ii ] = m[ ii ].scaffold_build_id();
		rot_vector[ n_geometric_constraints_ + ii ] = m[ ii ].upstream_conf_id();
	}

	RotamerComboCountMap::const_iterator iter = count_per_rotamer_combo_.find( rot_vector );
	if ( iter == count_per_rotamer_combo_.end() ) {
		return true;
	} else {
		return iter->second < limit_per_rotamer_combo_;
	}
}

/// @brief Note that a particular match has passed all the filters and will be output.
void
LimitHitsPerRotamerFilter::note_match_accepted(
	match const & m
)
{
	runtime_assert( m.size() == n_geometric_constraints_ );

	utility::vector1< Size > rot_vector( n_geometric_constraints_ * 2, 0 );
	for ( Size ii = 1; ii <= n_geometric_constraints_; ++ii ) {
		rot_vector[ ii ] = m[ ii ].scaffold_build_id();
		rot_vector[ n_geometric_constraints_ + ii ] = m[ ii ].upstream_conf_id();
	}

	RotamerComboCountMap::iterator iter = count_per_rotamer_combo_.find( rot_vector );
	if ( iter == count_per_rotamer_combo_.end() ) {
		count_per_rotamer_combo_[ rot_vector ] = 1;
	} else {
		++iter->second;
	}
}

/// @brief Erase all tracking data on which matches have already been output.
void
LimitHitsPerRotamerFilter::reset()
{
	count_per_rotamer_combo_.clear();
}

}
}
}
