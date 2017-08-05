// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pose/rna/leontis_westhof_util.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu

#include <core/pose/rna/leontis_westhof_util.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "core.pose.rna.leontis_westhof_util" );

namespace core {
namespace pose {
namespace rna {

using namespace core::chemical::rna;

/////////////////////////////////////////////////////////////////////////
LW_BaseDoubletOrientation
get_LW_orientation( BaseEdge const & edge1, BaseEdge const & edge2, BaseDoubletOrientation const & orientation )
{
	runtime_assert( edge1 == WATSON_CRICK || edge1 == HOOGSTEEN || edge1 == SUGAR );
	runtime_assert( edge2 == WATSON_CRICK || edge2 == HOOGSTEEN || edge2 == SUGAR );
	runtime_assert( orientation == ANTIPARALLEL || orientation == PARALLEL );

	LW_Table const & lw_table = get_leontis_westhof_table();
	LW_Table::const_iterator it = lw_table.find( std::make_pair( edge1, edge2 ) );
	runtime_assert( it != lw_table.end() );

	std::map< BaseDoubletOrientation, LW_BaseDoubletOrientation > const & orientation_map( it->second );
	return orientation_map.find( orientation )->second;
}

/////////////////////////////////////////////////////////////////////////
BaseDoubletOrientation
get_base_doublet_orientation_from_LW( BaseEdge const & edge1, BaseEdge const & edge2, LW_BaseDoubletOrientation const & lw_orientation )
{
	runtime_assert( edge1 == WATSON_CRICK || edge1 == HOOGSTEEN || edge1 == SUGAR );
	runtime_assert( edge2 == WATSON_CRICK || edge2 == HOOGSTEEN || edge2 == SUGAR );
	runtime_assert( lw_orientation == CIS || lw_orientation == TRANS );

	LW_Table const & lw_table = get_leontis_westhof_table();
	LW_Table::const_iterator it = lw_table.find( std::make_pair( edge1, edge2 ) );
	runtime_assert( it != lw_table.end() );

	std::map< BaseDoubletOrientation, LW_BaseDoubletOrientation > const & orientation_map( it->second );
	if ( orientation_map.find( ANTIPARALLEL )->second == lw_orientation ) return ANTIPARALLEL;
	runtime_assert( orientation_map.find( PARALLEL )->second == lw_orientation );
	return PARALLEL;
}

///////////////////////////////////////////////////////////////////
/// @details
/// From: RNA. 2001 Apr;7(4):499-512. "Geometric nomenclature and classification of RNA base pairs." Leontis NB, Westhof E.
/// Table formatted to match personal communication from N. Leontis to R. Das.
LW_Table const &
get_leontis_westhof_table() {
	static bool init( false );
	static LW_Table lw_table;
	if ( init ) return lw_table;

	lw_table[ std::make_pair( WATSON_CRICK, WATSON_CRICK ) ][ ANTIPARALLEL ] = CIS;
	lw_table[ std::make_pair( WATSON_CRICK, WATSON_CRICK ) ][     PARALLEL ] = TRANS;

	lw_table[ std::make_pair( WATSON_CRICK, HOOGSTEEN )    ][     PARALLEL ] = CIS;
	lw_table[ std::make_pair( WATSON_CRICK, HOOGSTEEN )    ][ ANTIPARALLEL ] = TRANS;

	lw_table[ std::make_pair( WATSON_CRICK, SUGAR )        ][ ANTIPARALLEL ] = CIS;
	lw_table[ std::make_pair( WATSON_CRICK, SUGAR )        ][     PARALLEL ] = TRANS;

	lw_table[ std::make_pair( HOOGSTEEN, HOOGSTEEN )       ][ ANTIPARALLEL ] = CIS;
	lw_table[ std::make_pair( HOOGSTEEN, HOOGSTEEN )       ][     PARALLEL ] = TRANS;

	lw_table[ std::make_pair( HOOGSTEEN, SUGAR )           ][     PARALLEL ] = CIS;
	lw_table[ std::make_pair( HOOGSTEEN, SUGAR )           ][ ANTIPARALLEL ] = TRANS;

	lw_table[ std::make_pair( SUGAR, SUGAR )               ][ ANTIPARALLEL ] = CIS;
	lw_table[ std::make_pair( SUGAR, SUGAR )               ][     PARALLEL ] = TRANS;

	// fill in the rest by symmetry.
	for ( auto const & elem : lw_table )  {
		std::pair< BaseEdge, BaseEdge > const & key = elem.first;
		lw_table[ std::make_pair( key.second, key.first ) ] = elem.second;
	}

	init = true;
	return lw_table;
}

} //rna
} //pose
} //core

