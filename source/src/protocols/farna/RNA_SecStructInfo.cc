// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/SecStruct/methods/RNA_BaseBasePotential.cc
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Rhiju Das

// Unit headers
#include <protocols/farna/RNA_SecStructInfo.hh>
#include <protocols/farna/RNA_SecStructInfo.fwd.hh>

// Package headers

// Project headers
#include <core/chemical/AA.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>

#include <utility/vector1.hh>


// Utility headers

// C++

///////////////////////////////////////////////////////
// Keep track of some base geometry that is
// useful for RNA scoring.
///////////////////////////////////////////////////////

namespace protocols {
namespace farna {

/// @details Copy constructors must copy all data, not just some...
RNA_SecStructInfo::RNA_SecStructInfo( RNA_SecStructInfo const & src ) :
	CacheableData()
{
	rna_secstruct_ = src.rna_secstruct_;
}

/// @details Pose must already contain a core::pose::datacache::CacheableDataType::RNA_SCORING_INFO object or this method will fail.
std::string const &
get_rna_secstruct( core::pose::Pose & pose )
{
	//using core::pose::datacache::CacheableDataType::RNA_SECSTRUCT_INFO;

	if ( !pose.data().has( core::pose::datacache::CacheableDataType::RNA_SECSTRUCT_INFO ) ) {
		set_rna_secstruct( pose, std::string( pose.total_residue(), 'X' ) );
	}

	return ( utility::pointer::static_pointer_cast< protocols::farna::RNA_SecStructInfo const > ( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::RNA_SECSTRUCT_INFO ) ) )->get_secstruct();
}

/// @details Either returns a non-const reference to the rna_scoring object already stored
/// in the pose, or creates a new rna scoring info object, places it in the pose, and returns
/// a non-const reference to it.
void
set_rna_secstruct(  core::pose::Pose & pose, std::string const & rna_secstruct_string )
{
	//using core::pose::datacache::CacheableDataType::RNA_SECSTRUCT_INFO;

	if ( pose.data().has( core::pose::datacache::CacheableDataType::RNA_SECSTRUCT_INFO ) ) {
		( utility::pointer::static_pointer_cast< protocols::farna::RNA_SecStructInfo > ( pose.data().get_ptr( core::pose::datacache::CacheableDataType::RNA_SECSTRUCT_INFO ) ) )->set_secstruct( rna_secstruct_string );
	}
	// else
	RNA_SecStructInfoOP rna_secstruct_info( new RNA_SecStructInfo( rna_secstruct_string ) );
	pose.data().set( core::pose::datacache::CacheableDataType::RNA_SECSTRUCT_INFO, rna_secstruct_info );
}


} //farna
} //protocols

