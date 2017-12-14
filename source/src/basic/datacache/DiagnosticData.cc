// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/datacache/DiagnosticData.hh
/// @brief
/// @author Phil Bradley

// unit headers
#include <basic/datacache/DiagnosticData.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION

namespace basic {
namespace datacache {

DiagnosticData::DiagnosticData( std::map < std::string, double > const & data_in ) : CacheableData(), data_(data_in) {}
DiagnosticData::~DiagnosticData() = default;
CacheableDataOP DiagnosticData::clone() const { return CacheableDataOP( new DiagnosticData(*this) ); }
std::map < std::string, double > const & DiagnosticData::data() const { return data_; }

DiagnosticDataOP DiagnosticData::shared_from_this() { return utility::pointer::static_pointer_cast<DiagnosticData>( CacheableData::shared_from_this() ); }

} // namespace datacache
} // namespace basic


#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
basic::datacache::DiagnosticData::DiagnosticData() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
basic::datacache::DiagnosticData::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( CEREAL_NVP( data_ ) ); // std::map<std::string, double>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
basic::datacache::DiagnosticData::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( data_ ); // std::map<std::string, double>
}

SAVE_AND_LOAD_SERIALIZABLE( basic::datacache::DiagnosticData );
CEREAL_REGISTER_TYPE( basic::datacache::DiagnosticData )

CEREAL_REGISTER_DYNAMIC_INIT( basic_datacache_DiagnosticData )
#endif // SERIALIZATION
