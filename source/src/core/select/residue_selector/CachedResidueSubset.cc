// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/CachedResidueSubset.cc
/// @brief  CacheableData wrapper for ResidueSubset storage
/// @author Tom Linsky ( tlinsky at uw dot edu )

// Unit headers
#include <core/select/residue_selector/CachedResidueSubset.hh>

#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>

#include <basic/Tracer.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>

#include <core/pose/Pose.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/map.hpp>
#include <cereal/types/vector.hpp>
#endif // SERIALIZATION

namespace core {
namespace select {
namespace residue_selector {

static basic::Tracer TR( "core.select.residue_selector.CachedResidueSubset" );

CachedResidueSubset & CachedResidueSubset::from_pose_datacache(core::pose::Pose & pose) {

	using namespace core::pose::datacache;

	// If the pose doesn't have cached subset data, add blank data
	if ( !pose.data().has( CacheableDataType::STORED_RESIDUE_SUBSET ) ) {
		CachedResidueSubsetOP blank_stored_subsets( new CachedResidueSubset );
		pose.data().set( CacheableDataType::STORED_RESIDUE_SUBSET, blank_stored_subsets );
	}

	// Grab a reference to the data
	debug_assert( utility::pointer::dynamic_pointer_cast< CachedResidueSubset >( pose.data().get_ptr( CacheableDataType::STORED_RESIDUE_SUBSET ) ) );
	return
		*( utility::pointer::static_pointer_cast< CachedResidueSubset >( pose.data().get_ptr( CacheableDataType::STORED_RESIDUE_SUBSET ) ) );
}

// @brief default constructor
CachedResidueSubset::CachedResidueSubset():
	basic::datacache::CacheableData(),
	subsets_()
{}

void
CachedResidueSubset::set_subset( ResidueSubsetCOP const subset, std::string const & name )
{
	debug_assert( subset );
	subsets_[ name ] = utility::pointer::make_shared< BoolVector >( subset->begin(), subset->end() );
}

ResidueSubsetCOP
CachedResidueSubset::get_subset( std::string const & name ) const
{
	auto subset = subsets_.find( name );
	if ( subset == subsets_.end() ) {
		std::stringstream msg;
		msg << "CachedResidueSubset: No cached residue subset named " << name << " was found in the pose.  Available subsets are: ";
		for ( auto const & subset2 : subsets_ ) {
			msg << subset2.first << " ";
		}
		msg << std::endl;
		throw CREATE_EXCEPTION(utility::excn::BadInput,  msg.str() );
	}
	return utility::pointer::make_shared< ResidueSubset >( subset->second->begin(), subset->second->end() );
}

// @brief check to see if the task you're interested in is in the object
bool
CachedResidueSubset::has_subset( std::string const & name ) const
{
	return ( subsets_.find( name ) != subsets_.end() );
}

basic::datacache::CacheableDataOP
CachedResidueSubset::clone() const
{
	return utility::pointer::make_shared< CachedResidueSubset > ( *this );
}

basic::datacache::CacheableDataOP
CachedResidueSubset::fresh_instance() const
{
	return utility::pointer::make_shared< CachedResidueSubset >();
}

} // residue_selector
} // select
} // core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::select::residue_selector::CachedResidueSubset::save( Archive & arc ) const
{
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( CEREAL_NVP( subsets_ ) ); // std::map< std::string, ResidueSubsetOP >
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::select::residue_selector::CachedResidueSubset::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( subsets_ );
}

SAVE_AND_LOAD_SERIALIZABLE( core::select::residue_selector::CachedResidueSubset );
CEREAL_REGISTER_TYPE( core::select::residue_selector::CachedResidueSubset )
CEREAL_REGISTER_DYNAMIC_INIT( core_select_residue_selector_CachedResidueSubset )
#endif // SERIALIZATION
