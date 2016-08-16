// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/Hit.cc
/// @brief  Constructors for match_dspos1
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

// Unit headers
#include <protocols/match/Hit.hh>

namespace protocols {
namespace match {

bool
upstream_hit::operator < ( upstream_hit const & rhs ) const
{
	if ( first_[1] < rhs.first_[1] ) return true;
	else if ( first_[1] > rhs.first_[1] ) return false;

	if ( first_[2] < rhs.first_[2] ) return true;
	else if ( first_[2] > rhs.first_[2] ) return false;

	if ( first_[3] < rhs.first_[3] ) return true;
	else  if ( first_[3] > rhs.first_[3] ) return false;

	return false;
}

bool
upstream_hit::operator == ( upstream_hit const & rhs ) const
{
	if ( ( first_[1] == rhs.first_[1] )
			&&( first_[2] == rhs.first_[2] )
			&&( first_[3] == rhs.first_[3] ) ) return true;

	return false;
}

bool
downstream_hit::operator < ( downstream_hit const & rhs ) const
{
	if ( downstream_conf_id_ < rhs.downstream_conf_id_ ) return true;
	else if ( downstream_conf_id_ > rhs.downstream_conf_id_ ) return false;

	for ( core::Size i =1; i <= 6; ++i ) {
		if ( second_[i] < rhs.second_[i] ) return true;
		else if (  second_[i] > rhs.second_[i] ) return false;
	}
	return false;
}

bool
downstream_hit::operator == ( downstream_hit const & rhs ) const
{
	if ( ( downstream_conf_id_ == rhs.downstream_conf_id_ )
			&&( second_[1] == rhs.second_[1] )
			&&( second_[2] == rhs.second_[2] )
			&&( second_[3] == rhs.second_[3] )
			&&( second_[4] == rhs.second_[4] )
			&&( second_[5] == rhs.second_[5] )
			&&( second_[6] == rhs.second_[6] ) ) return true;

	return false;
}

match_dspos1::match_dspos1() :
	originating_geom_cst_for_dspos( 0 ),
	downstream_conf_id( 0 ),
	dspos( /* 0 */ )
{}

match_dspos1::match_dspos1( core::Size n_geometric_constraints ) : upstream_hits( n_geometric_constraints ) {}

match_dspos1::match_dspos1(
	match const & m,
	core::Size geomcst_specifying_dspos
) :
	upstream_hits( m.size() ),
	originating_geom_cst_for_dspos( geomcst_specifying_dspos ),
	downstream_conf_id( m[ geomcst_specifying_dspos ].first()[ 4 ] ),
	dspos( m[ geomcst_specifying_dspos ].second() )
{
	for ( core::Size ii = 1; ii <= m.size(); ++ii ) {
		upstream_hits[ ii ] = upstream_hit( m[ ii ] );
	}
}

Hit
fake_hit( upstream_hit const & uhit ) {
	Hit hit;
	hit.first()[ 1 ] = uhit.scaffold_build_id();
	hit.first()[ 2 ] = uhit.upstream_conf_id();
	hit.first()[ 3 ] = uhit.external_geom_id();
	hit.first()[ 4 ] = 0;
	std::fill( hit.second().begin(), hit.second().end(), 0.0 );
	return hit;
}

Hit
fake_hit( downstream_hit const & dhit ) {
	Hit hit;
	hit.first()[ 1 ] = 0;
	hit.first()[ 2 ] = 0;
	hit.first()[ 3 ] = 0;
	hit.first()[ 4 ] = dhit.downstream_conf_id();
	hit.second() = dhit.second();
	return hit;
}


/// @brief Create a hit with the full data from a given match_dspos1 representing
/// the upstream conformation from the originating_geom_cst and its
/// description of the downstream position.
Hit full_hit( match_dspos1 const & m )
{
	core::Size const geomcst_w_dspos = m.originating_geom_cst_for_dspos;
	Hit fullhit;
	fullhit.first()[ 1 ] = m.upstream_hits[ geomcst_w_dspos ].scaffold_build_id();
	fullhit.first()[ 2 ] = m.upstream_hits[ geomcst_w_dspos ].upstream_conf_id();
	fullhit.first()[ 3 ] = m.upstream_hits[ geomcst_w_dspos ].external_geom_id();
	fullhit.first()[ 4 ] = m.downstream_conf_id;
	fullhit.second()     = m.dspos;
	return fullhit;
}

}
}

