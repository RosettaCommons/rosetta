// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody_design/SasaMethod.cc
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <core/scoring/sasa/SasaMethod.hh>

//Core headers:
#include <core/chemical/ResidueType.hh>
#include <core/chemical/Atom.hh>
#include <core/conformation/Residue.hh>

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
{}

SasaMethod::~SasaMethod()= default;

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

/// @brief Given the name of the SasaMethodHPMode, get the mode.
/// @returns SasaMethodHPMode::INVALID_MODE if the string can't be parsed; the correct
/// SasaMethodHPMode if it can.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
SasaMethodHPMode
SasaMethod::sasa_metric_mode_from_name(
	std::string const &mode_name
) {
	for ( core::Size i(1); i < static_cast<core::Size>( SasaMethodHPMode::END_OF_LIST ); ++i ) {
		if ( sasa_metric_name_from_mode( static_cast< SasaMethodHPMode >(i) ) == mode_name ) {
			return static_cast< SasaMethodHPMode >( i );
		}
	}
	return SasaMethodHPMode::INVALID_MODE;
}

/// @brief Given the SasaMethodHPMode, get the name.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
std::string
SasaMethod::sasa_metric_name_from_mode(
	SasaMethodHPMode const mode
) {
	switch( mode ) {
	case SasaMethodHPMode::ALL_SASA :
		return std::string( "all_sasa" );
	case SasaMethodHPMode::POLAR_SASA :
		return std::string( "polar_sasa" );
	case SasaMethodHPMode::HYDROPHOBIC_SASA :
		return std::string( "hydrophobic_sasa" );
	default :
		utility_exit_with_message( "Error in SasaMetric::sasa_metric_name_from_mode(): Invalid SasaMethodHPMode provided." );
	};
	return std::string(""); //Keep compiler happy.
}

/// @brief Construct a comma-separeted string listing all of the sasa metric modes.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
std::string
SasaMethod::list_sasa_method_hp_modes() {
	std::stringstream ostream;
	for ( core::Size i(1); i < static_cast<core::Size>(SasaMethodHPMode::END_OF_LIST); ++i ) {
		if ( i > 1 ) {
			ostream << ", ";
		}
		ostream << sasa_metric_name_from_mode( static_cast<SasaMethodHPMode>(i) );
	}
	return ostream.str();
}

/// @brief Set whether we're counting all SASA (default), polar SASA, or hydrophobic SASA.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
void
SasaMethod::set_sasa_method_hp_mode(
	SasaMethodHPMode const mode_in
) {
	runtime_assert_string_msg( static_cast< core::Size >(mode_in) > 0 && mode_in < SasaMethodHPMode::END_OF_LIST, "Error in SasaMethod::set_sasa_method_hp_mode():  Mode not recognized!" );
	sasa_method_hp_mode_ = mode_in;
}

void
SasaMethod::set_use_big_polar_hydrogen(bool big_polar_h){
	use_big_polar_H_ = big_polar_h;
}

/// @brief Given a residue, an atom index, and a SasaMethodHPMode, determine whether the atom is one to skip (returns true)
/// or count (returns false).
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
bool
SasaMethod::skip_atom(
	core::conformation::Residue const & rsd,
	core::Size const atom_index,
	SasaMethodHPMode const hp_mode
) {
	switch( hp_mode ) {
	case SasaMethodHPMode::ALL_SASA :
		return false;
	case SasaMethodHPMode::POLAR_SASA :
		{
		core::chemical::Atom const & curatom( rsd.type().atom(atom_index) );
		if ( curatom.heavyatom_has_polar_hydrogens() || curatom.is_acceptor() || curatom.is_polar_hydrogen() ) {
			return false;
		} else {
			return true;
		}
	}
	case SasaMethodHPMode::HYDROPHOBIC_SASA :
		{
		core::chemical::Atom const & curatom( rsd.type().atom(atom_index) );
		if ( curatom.heavyatom_has_polar_hydrogens() || curatom.is_acceptor() || curatom.is_polar_hydrogen() ) {
			return true;
		} else {
			return false;
		}
	}
	default :
		utility_exit_with_message( "Error in SasaMethod::skip_atom(): Unknown SasaMethodHPMode!" );
	};
	return false; //Keep the compiler happy.
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
	arc( CEREAL_NVP( sasa_method_hp_mode_ ) ); //SasaMethodHPMode enum class
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::sasa::SasaMethod::load( Archive & arc ) {
	arc( probe_radius_ ); // Real
	arc( radii_set_ ); // enum core::scoring::sasa::SasaRadii
	arc( include_probe_radius_ ); // _Bool
	arc( use_big_polar_H_ ); // _Bool
	arc( sasa_method_hp_mode_ ); //SasaMethodHPMode enum class
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::sasa::SasaMethod );
CEREAL_REGISTER_TYPE( core::scoring::sasa::SasaMethod )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_sasa_SasaMethod )
#endif // SERIALIZATION
