// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ./src/protocols/fldsgn/SheetFoldTypeManager.cc
/// @brief
/// @author Nobuyasu Koga ( nobuyasu@u.washington.edu )

// unit headers
#include <protocols/fldsgn/topology/SheetFoldTypeManager.hh>

// project headers

// utility headers
#include <utility/exit.hh>

// C++ headers
#include <iostream>

namespace protocols {
namespace fldsgn {
namespace topology {

SheetFoldTypeManager::SheetFoldTypeManager():
	initialized_( false )
{
	initialize();
}

void
SheetFoldTypeManager::initialize()
{
	setup_foldtype_names();
	setup_foldtype_strand_pairings();
	initialized_ = true;
}

/// @brief initialize the SheetFoldType name vector and map
void
SheetFoldTypeManager::setup_foldtype_names()
{
	/// 2strands
	name2foldtype_[ "BAB" ] = BABx1;

	/// 3strands
	// parallel
	name2foldtype_[ "RosI" ] = RosI;
	name2foldtype_[ "RosO" ] = RosO;
	name2foldtype_[ "BABx2" ] =BABx2;

	// mixture
	name2foldtype_[ "Thio" ] = Thio;
	name2foldtype_[ "BFr" ] = BFr;
	name2foldtype_[ "EFr" ] = EFr;

	// anti-parallel
	name2foldtype_[ "CFr" ] = CFr;
	name2foldtype_[ "DFr" ] = DFr;

	// 4strands
	// parallel
	name2foldtype_[ "Rsmn2x2" ] = Rsmn2x2;
	name2foldtype_[ "Rsmn3x3_Half" ] = Rsmn3x3_Half;
	name2foldtype_[ "BABx3" ] =BABx3;

	// mixture
	name2foldtype_[ "PG_like" ] = PG_like;
	name2foldtype_[ "Thioredoxin" ] = Thioredoxin;
	name2foldtype_[ "BAB_CFr" ] = BAB_CFr;
	name2foldtype_[ "DFr_BAB" ] = DFr_BAB;
	name2foldtype_[ "BEFr" ] = BEFr;

	// anti-parallel
	name2foldtype_[ "Fd_like" ] = Fd_like;
	name2foldtype_[ "RFd_like" ] = RFd_like;
	name2foldtype_[ "CDFr" ] = CDFr;
	name2foldtype_[ "HPN_CFr" ] = HPN_CFr;
	name2foldtype_[ "DFr_HPN" ] = DFr_HPN;

	// 5strands
	// parallel
	name2foldtype_[ "Flavodoxin" ] = Flavodoxin;
	name2foldtype_[ "Ploop2x3" ] = Ploop2x3;

	// mixture
	name2foldtype_[ "RNAseH" ] = RNAseH;

	// anti-parallel
	name2foldtype_[ "Top7" ] = Top7;

	// 6strands
	name2foldtype_[ "Rsmn3x3" ] = Rsmn3x3;
	name2foldtype_[ "Ploop3x3" ] = Ploop3x3;


	name2foldtype_[ "UNFOLD" ] = UNFOLD;


	name2foldtype_[ "NO_STRANDS" ] = NO_STRANDS;


	name2foldtype_[ "UNKNOWN" ] = UNKNOWN;

	assert( name2foldtype_.size() == n_fold_types );

	foldtype2name_.resize( n_fold_types );

	for ( std::map< String, SheetFoldType >::const_iterator iter = name2foldtype_.begin(),
			iter_end = name2foldtype_.end(); iter != iter_end; ++iter ) {
		foldtype2name_[ iter->second ] = iter->first;
	}

}

/// @brief give a string name of SheetFoldType and return its enum type
SheetFoldType
SheetFoldTypeManager::foldtype_from_name( String const & name )
{
	if ( ! initialized_ ) initialize();
	std::map< String, SheetFoldType >::const_iterator iter( name2foldtype_.find( name ) );
	if ( iter == name2foldtype_.end() ) {
		utility_exit_with_message("unrecognized foldtype type "+name);
	}
	return iter->second;
}

/// @brief give a SheetFoldType and return its string name
std::string
SheetFoldTypeManager::name_from_foldtype( SheetFoldType foldtype )
{
	if ( ! initialized_ ) initialize();
	return foldtype2name_[ foldtype ];
}

/// @brief check whether the string name of strand pairings is in SheetFoldType or not
bool
SheetFoldTypeManager::is_foldtype( String const & name )
{
	if ( ! initialized_ ) initialize();
	std::map< String, SheetFoldType >::const_iterator iter( name2foldtype_.find( name ) );
	return iter != name2foldtype_.end();
}

/// @brief initialize the map of strand pairings and SheetFoldType
void
SheetFoldTypeManager::setup_foldtype_strand_pairings()
{
	// 2strands
	spairs2foldtype_[ "1-2.P" ] = BABx1;

	/// 3strands
	spairs2foldtype_[ "1-2.P;1-3.P" ] = RosI;
	spairs2foldtype_[ "1-3.P;2-3.P" ] = RosO;
	spairs2foldtype_[ "1-2.P;2-3.P" ] = BABx2;

	spairs2foldtype_[ "1-2.P;1-3.A" ] = Thio;
	spairs2foldtype_[ "1-2.P;2-3.A" ] = BFr;
	spairs2foldtype_[ "1-2.A;2-3.P" ] = EFr;

	spairs2foldtype_[ "1-3.A;2-3.A" ] = CFr;
	spairs2foldtype_[ "1-2.A;1-3.A" ] = DFr;

	// 4strands
	spairs2foldtype_[ "1-2.P;1-3.P;3-4.P" ] = Rsmn2x2;
	spairs2foldtype_[ "1-2.P;1-4.P;2-3.P" ] = Rsmn3x3_Half;
	spairs2foldtype_[ "1-2.P;2-3.P;3-4.P" ] = BABx3;

	spairs2foldtype_[ "1-2.A;1-4.P;3-4.A" ] = PG_like;
	spairs2foldtype_[ "1-2.P;1-3.A;3-4.A" ] = Thioredoxin;
	//spairs2foldtype_[ "1-2.P;1-3.A;3-4.A" ] = L30E_like;
	spairs2foldtype_[ "1-2.P;2-4.A;3-4.A" ] = BAB_CFr;
	spairs2foldtype_[ "1-2.A;1-3.A;3-4.P" ] = DFr_BAB;
	spairs2foldtype_[ "1-2.A;2-3.P;3-4.A" ] = BEFr;

	spairs2foldtype_[ "1-3.A;1-4.A;2-3.A" ] = Fd_like;
	spairs2foldtype_[ "1-4.A;2-3.A;2-4.A" ] = RFd_like;
	spairs2foldtype_[ "1-3.A;2-3.A;2-4.A" ] = CDFr;
	spairs2foldtype_[ "1-2.A;2-4.A;3-4.A" ] = HPN_CFr;
	spairs2foldtype_[ "1-2.A;1-3.A;3-4.A" ] = DFr_HPN;

	// 5strands
	spairs2foldtype_[ "1-2.P;1-3.P;3-4.P;4-5.P" ] = Flavodoxin;
	spairs2foldtype_[ "1-3.P;1-4.P;2-3.P;4-5.P" ] = Ploop2x3;

	spairs2foldtype_[ "1-2.A;2-3.A;1-4.P;4-5.P" ] = RNAseH;

	spairs2foldtype_[ "1-2.A;2-4.P;3-5.A;4-5.A" ] = Top7;

	// 6strands
	spairs2foldtype_[ "1-2.P;1-4.P;2-3.P;4-5.P;5-6.P" ] = Rsmn3x3;
	spairs2foldtype_[ "1-4.P;1-5.P;2-3.P;2-4.P;5-6.P" ] = Ploop3x3;


	spairs2foldtype_[ "UNFOLD" ] = UNFOLD;


	spairs2foldtype_[ "" ] = NO_STRANDS;


	spairs2foldtype_[ "UNKNOWN" ] = UNKNOWN;

	assert( spairs2foldtype_.size() == n_fold_types );

	foldtype2spairs_.resize( n_fold_types );

	for ( std::map< String, SheetFoldType >::const_iterator iter = spairs2foldtype_.begin(),
			iter_end = spairs2foldtype_.end(); iter != iter_end; ++iter ) {
		foldtype2spairs_[ iter->second ] = iter->first;

	}

}

/// @brief give a string of strand_pairings and return its enum type
SheetFoldType
SheetFoldTypeManager::foldtype_from_spairs( String const & spairs )
{
	if ( ! initialized_ ) initialize();
	std::map< String, SheetFoldType >::const_iterator iter( spairs2foldtype_.find( spairs ) );
	if ( iter == spairs2foldtype_.end() ) {
		return UNKNOWN;
		//utility_exit_with_message("unrecognized foldtype type "+spairs);
	}
	return iter->second;
}

/// @brief give a SheetFoldType and return its string of strand_pairings

std::string
SheetFoldTypeManager::spairs_from_foldtype( SheetFoldType foldtype )
{
	if ( ! initialized_ ) initialize();
	return foldtype2spairs_[ foldtype ];
}

/// @brief check whether the string spairs of strand_pairings is in SheetFoldType or not
bool
SheetFoldTypeManager::is_sparis_foldtype( String const & spairs )
{
	if ( ! initialized_ ) initialize();
	std::map< String, SheetFoldType >::const_iterator iter( spairs2foldtype_.find( spairs ) );
	return iter != spairs2foldtype_.end();
}

} // namespace topology
} // namespace fldsgn
} // namespace protocols
