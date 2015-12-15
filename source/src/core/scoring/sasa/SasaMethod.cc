// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody_design/SasaMethod.cc
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <core/scoring/sasa/SasaMethod.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {
namespace sasa {

SasaMethod::SasaMethod(Real probe_radius, SasaRadii radii_set):
	utility::pointer::ReferenceCount(),
	probe_radius_(probe_radius),
	radii_set_(radii_set),
	include_probe_radius_(true),
	use_big_polar_H_(false)
{

}

SasaMethod::~SasaMethod(){}

void
SasaMethod::set_include_probe_radius_in_calc(bool include_probe_radius) {
	include_probe_radius_ = include_probe_radius;
}

void
SasaMethod::set_probe_radius(Real probe_radius){
	probe_radius_ = probe_radius;
}

void
SasaMethod::set_radii_set(SasaRadii radii_set) {
	radii_set_ = radii_set;
}

void
SasaMethod::set_use_big_polar_hydrogen(bool big_polar_h){
	use_big_polar_H_ = big_polar_h;
}


}
} //scoring
} //core

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::scoring::sasa::SasaMethod::SasaMethod() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::sasa::SasaMethod::save( Archive & arc ) const {
	arc( CEREAL_NVP( probe_radius_ ) ); // Real
	arc( CEREAL_NVP( radii_set_ ) ); // enum core::scoring::sasa::SasaRadii
	arc( CEREAL_NVP( include_probe_radius_ ) ); // _Bool
	arc( CEREAL_NVP( use_big_polar_H_ ) ); // _Bool
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::sasa::SasaMethod::load( Archive & arc ) {
	arc( probe_radius_ ); // Real
	arc( radii_set_ ); // enum core::scoring::sasa::SasaRadii
	arc( include_probe_radius_ ); // _Bool
	arc( use_big_polar_H_ ); // _Bool
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::sasa::SasaMethod );
CEREAL_REGISTER_TYPE( core::scoring::sasa::SasaMethod )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_sasa_SasaMethod )
#endif // SERIALIZATION
