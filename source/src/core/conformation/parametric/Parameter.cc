// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/conformation/parametric/Parameter.cc
/// @brief  A class for holding a single parameter for parametric backbone generation.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Unit header
#include <core/conformation/parametric/Parameter.hh>

// Package headers

// Project headers

// Basic headers
#include <basic/Tracer.hh>

// Numeric headers

// Utility Headers
#include <utility/vector1.hh>


#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace conformation {
namespace parametric {

static basic::Tracer TR( "core.conformation.parametric.Parameter" );

/// @brief Constructor.
///
Parameter::Parameter() :
	parameter_name_(""),
	parameter_description_(""),
	short_parameter_description_(""),
	parameter_units_("dimensionless"),
	parameter_type_(PT_invalid_type),
	can_be_set_(true),
	can_be_copied_(true),
	can_be_sampled_(true),
	can_be_perturbed_(true),
	global_for_parameters_set_(true),
	copy_suffix_("copied_from"),
	copy_from_parameters_index_(0),
	value_set_(false)
	//TODO -- initialize variables here.
{
}

Parameter::Parameter( Parameter const & src ) :
	utility::pointer::ReferenceCount(),
	utility::pointer::enable_shared_from_this< Parameter >(),
	parameter_name_(src.parameter_name_),
	parameter_description_( src.parameter_description_ ),
	short_parameter_description_(src.short_parameter_description_),
	parameter_units_(src.parameter_units_),
	parameter_type_(src.parameter_type_),
	can_be_set_(src.can_be_set_),
	can_be_copied_(src.can_be_copied_),
	can_be_sampled_(src.can_be_sampled_),
	can_be_perturbed_(src.can_be_perturbed_),
	global_for_parameters_set_(src.global_for_parameters_set_),
	copy_suffix_(src.copy_suffix_),
	copy_from_parameters_index_(src.copy_from_parameters_index_),
	value_set_(src.value_set_)
{
}

/// @brief Pure virtual destructor.
/// @details Counter-intuitively, C++ requires pure virtual destructors to be implemented.
Parameter::~Parameter() {}

/// @brief Set parameter name.
void
Parameter::set_parameter_name(
	std::string const &name_in
) {
	debug_assert(!name_in.empty());
	parameter_name_ = name_in;
}

/// @brief Set parameter description.
/// @details This is a lay-language description used for annotating user output.  It should be a short phrase with no capitals.
void
Parameter::set_parameter_description(
	std::string const & description_in
) {
	debug_assert(!description_in.empty());
	parameter_description_ = description_in;
}

/// @brief Set short parameter description.
/// @details This is a one or two word lay-language description used for annotating user output.  The first word should be capitalized.
void
Parameter::set_short_parameter_description(
	std::string const & short_description_in
) {
	debug_assert(!short_description_in.empty());
	short_parameter_description_ = short_description_in;
}

/// @brief Set parameter units.
void
Parameter::set_parameter_units(
	std::string const & units_in
) {
	debug_assert(!units_in.empty());
	parameter_units_ = units_in;
}

/// @brief Set the parameter type.
void
Parameter::set_parameter_type(
	ParameterType const type_in
) {
	parameter_type_ = type_in;
}

/// @brief Set whether this parameter is one that can be set directly by the user, whether it can be sampled, and whether it can be perturbed.
/// @details Default true in all cases.  If setting is false, the parameter value must be read from the database, e.g. from a Crick parameters file.
void
Parameter::set_can_be_set_sampled_perturbed_copied(
	bool const can_be_set,
	bool const can_be_copied,
	bool const can_be_sampled,
	bool const can_be_perturbed
) {
	can_be_set_ = can_be_set;
	can_be_copied_ = can_be_copied;
	can_be_sampled_ = can_be_sampled;
	can_be_perturbed_ = can_be_perturbed;
}

/// @brief Set whether this property is meant to be global for a parameters set.
void
Parameter::set_global_for_parameters_set(
	bool const setting
) {
	global_for_parameters_set_ = setting;
}

/// @brief Set the suffix for the copy tag (e.g. "copies_helix" in the case of helical bundles).
/// @details The leading underscore should be omitted.
void
Parameter::set_copy_suffix(
	std::string const suffix_in
) {
#ifndef NDEBUG
	debug_assert( !suffix_in.empty() );
	std::string const suffix_in_stripped( utility::strip( suffix_in ) );
	debug_assert( suffix_in_stripped == suffix_in );
#endif
	copy_suffix_ = suffix_in;
}


/// @brief Set the index of the Parameters object from which we will copy this parameter's value.
/// @details A setting of zero means no copying.
void
Parameter::set_copy_from_parameters_index(
	core::Size const setting
) {
	reset_sampling_and_perturbation_settings();
	reset_value_settings();
	copy_from_parameters_index_ = setting;
}

/// @brief Reset the copying options.
void
Parameter::reset_copying_settings() {
	copy_from_parameters_index_ = 0;
}


} // namespace parametric
} // namespace conformation
} // namespace core


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::conformation::parametric::Parameter::save( Archive & arc ) const {
	arc( CEREAL_NVP( parameter_name_ ) ); //std::string
	arc( CEREAL_NVP( parameter_description_ ) ); //std::string
	arc( CEREAL_NVP( short_parameter_description_ ) ); //std::string
	arc( CEREAL_NVP( parameter_units_ ) ); //std::string
	arc( CEREAL_NVP( parameter_type_ ) ); //ParameterType enum
	arc( CEREAL_NVP( can_be_set_ ) ); //bool
	arc( CEREAL_NVP( can_be_copied_ ) ); //bool
	arc( CEREAL_NVP( can_be_sampled_ ) ); //bool
	arc( CEREAL_NVP( can_be_perturbed_ ) ); //bool
	arc( CEREAL_NVP( global_for_parameters_set_ ) ); //bool
	arc( CEREAL_NVP( copy_suffix_ ) ); //std::string
	arc( CEREAL_NVP( copy_from_parameters_index_ ) ); //core::Size
	arc( CEREAL_NVP( value_set_ ) ); //bool
	//TODO -- archive private member variables here.
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::conformation::parametric::Parameter::load( Archive & arc ) {
	arc( parameter_name_ ); //std::string
	arc( parameter_description_ ); //std::string
	arc( short_parameter_description_ ); //std::string
	arc( parameter_units_ ); //std::string
	arc( parameter_type_ ); //ParameterType enum
	arc( can_be_set_ ); //bool
	arc( can_be_copied_ ); //bool
	arc( can_be_sampled_ ); //bool
	arc( can_be_perturbed_ ); //bool
	arc( global_for_parameters_set_ ); //bool
	arc( copy_suffix_ ); //std::string
	arc( copy_from_parameters_index_ ); //core::Size
	arc( value_set_ ); //bool
	//TODO -- de-archive private member variables here.
}

SAVE_AND_LOAD_SERIALIZABLE( core::conformation::parametric::Parameter );
CEREAL_REGISTER_TYPE( core::conformation::parametric::Parameter )

CEREAL_REGISTER_DYNAMIC_INIT( core_conformation_parametric_Parameter )
#endif // SERIALIZATION
