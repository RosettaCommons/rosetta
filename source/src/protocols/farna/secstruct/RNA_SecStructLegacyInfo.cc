// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/SecStruct/methods/RNA_BaseBasePotential.cc
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Rhiju Das

// Unit headers
#include <protocols/farna/secstruct/RNA_SecStructLegacyInfo.hh>
#include <protocols/farna/secstruct/RNA_SecStructLegacyInfo.fwd.hh>

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

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace farna {
namespace secstruct {

/// @details Copy constructors must copy all data, not just some...
RNA_SecStructLegacyInfo::RNA_SecStructLegacyInfo( RNA_SecStructLegacyInfo const & src ) :
	CacheableData()
{
	rna_secstruct_legacy_ = src.rna_secstruct_legacy_;
}

/// @details Pose must already contain a core::pose::datacache::CacheableDataType::RNA_SCORING_INFO object or this method will fail.
std::string const &
get_rna_secstruct_legacy( core::pose::Pose & pose )
{
	if ( !pose.data().has( core::pose::datacache::CacheableDataType::RNA_SECSTRUCT_INFO ) ) {
		set_rna_secstruct_legacy( pose, std::string( pose.size(), 'X' ) );
	}

	return ( utility::pointer::static_pointer_cast< RNA_SecStructLegacyInfo const > ( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::RNA_SECSTRUCT_INFO ) ) )->get_secstruct();
}


/// @details Pose must already contain a core::pose::datacache::CacheableDataType::RNA_SCORING_INFO object or this method will fail.
std::string const &
get_rna_secstruct_legacy_from_const_pose( core::pose::Pose const & pose )
{
	runtime_assert ( pose.data().has( core::pose::datacache::CacheableDataType::RNA_SECSTRUCT_INFO ) );
	return ( utility::pointer::static_pointer_cast< RNA_SecStructLegacyInfo const > ( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::RNA_SECSTRUCT_INFO ) ) )->get_secstruct();
}

/// @details Either returns a non-const reference to the rna_scoring object already stored
/// in the pose, or creates a new rna scoring info object, places it in the pose, and returns
/// a non-const reference to it.
void
set_rna_secstruct_legacy(  core::pose::Pose & pose, std::string const & rna_secstruct_legacy_string )
{
	if ( pose.data().has( core::pose::datacache::CacheableDataType::RNA_SECSTRUCT_INFO ) ) {
		( utility::pointer::static_pointer_cast< RNA_SecStructLegacyInfo > ( pose.data().get_ptr( core::pose::datacache::CacheableDataType::RNA_SECSTRUCT_INFO ) ) )->set_secstruct( rna_secstruct_legacy_string );
	}
	// else
	RNA_SecStructLegacyInfoOP rna_secstruct_legacy_info( new RNA_SecStructLegacyInfo( rna_secstruct_legacy_string ) );
	pose.data().set( core::pose::datacache::CacheableDataType::RNA_SECSTRUCT_INFO, rna_secstruct_legacy_info );
}


/// @details remove RNA scoring info object so that it will be re-initialized if necessary
void
clear_rna_secstruct_legacy_info( core::pose::Pose & pose )
{
	pose.data().clear( core::pose::datacache::CacheableDataType::RNA_SECSTRUCT_INFO );
}

} //secstruct
} //farna
} //protocols


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::farna::secstruct::RNA_SecStructLegacyInfo::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( CEREAL_NVP( rna_secstruct_legacy_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::farna::secstruct::RNA_SecStructLegacyInfo::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( rna_secstruct_legacy_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::farna::secstruct::RNA_SecStructLegacyInfo );
CEREAL_REGISTER_TYPE( protocols::farna::secstruct::RNA_SecStructLegacyInfo )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_farna_secstruct_legacy_RNA_SecStructLegacyInfo )
#endif // SERIALIZATION
