// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/helical_bundle/parameters/OmegaBundleParameter.cc
/// @brief A class for the omega paremeter, derived from the generic RealValuedParameter class.  Omega has a few additional options that can be configured, for pitch-copying vs. value-copying.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

// Associated headers:
#include <protocols/helical_bundle/parameters/OmegaBundleParameter.hh>

// Core headers:
#include <core/conformation/parametric/Parameters.hh>
#include <core/conformation/parametric/Parameter.hh>
#include <core/conformation/parametric/RealValuedParameter.hh>

// Protocols headers:
#include <protocols/helical_bundle/BundleParametrizationCalculator.hh>

// Utility headers:
#include <utility/tag/Tag.hh>

// Basic headers:
#include <basic/Tracer.hh>

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


static basic::Tracer TR( "protocols.helical_bundle.parameters.OmegaBundleParameter" );


namespace protocols {
namespace helical_bundle {
namespace parameters {

/// @brief Default constructor.
OmegaBundleParameter::OmegaBundleParameter() :
	core::conformation::parametric::RealValuedParameter(),
	pitch_copying_mode_(false)
{}

/// @brief Copy constructor.
OmegaBundleParameter::OmegaBundleParameter(
	OmegaBundleParameter const &src
) :
	core::conformation::parametric::RealValuedParameter( src ),
	pitch_copying_mode_(src.pitch_copying_mode_)
{}

/// @details Destructor
OmegaBundleParameter::~OmegaBundleParameter() = default;

/// @brief Clone operator.
/// @details Make a copy of this object and return a smart pointer to the copy.
core::conformation::parametric::ParameterOP
OmegaBundleParameter::clone() const {
	return core::conformation::parametric::ParameterOP( OmegaBundleParameterOP( new OmegaBundleParameter( *this ) ) );
}

/// @brief Given another parameter of the same type, copy its value.  This does *not* set value_set_ to true.
/// @details Performs type checking in debug mode.  Note that this version has options for copying pitch or
/// copying the value directly.
/// @returns Returns TRUE for failure, FALSE for success.
bool
OmegaBundleParameter::copy_value_from_parameter(
	core::conformation::parametric::ParameterCOP other_parameter,
	core::conformation::parametric::ParametersCOP other_parameter_collection,
	core::conformation::parametric::ParametersCOP this_parameter_collection
) {
#ifndef NDEBUG
	OmegaBundleParameterCOP other_parameter_omega( utility::pointer::dynamic_pointer_cast< OmegaBundleParameter const >( other_parameter ) );
	debug_assert( other_parameter_omega != nullptr );
#endif

	if ( pitch_copying_mode_ ) {
		return copy_pitch_from_parameter( other_parameter, other_parameter_collection, this_parameter_collection );
	}
	return core::conformation::parametric::RealValuedParameter::copy_value_from_parameter( other_parameter, other_parameter_collection, this_parameter_collection );
}

//////////////// PARSE FUNCTIONS ///////////////////

/// @brief Return the XSD information for this parameter.
/// @details Calls the equivalent function for RealValuedParameter, then adds additional settings for pitch copying.
/// @note Must be implemented by derived classes.
void
OmegaBundleParameter::provide_xsd_information(
	utility::tag::AttributeList & xsd,
	bool const provide_setting,
	bool const provide_copying,
	bool const provide_grid_sampling,
	bool const provide_perturbation
) const {
	core::conformation::parametric::RealValuedParameter::provide_xsd_information( xsd, provide_setting, provide_copying, provide_grid_sampling, provide_perturbation );
	if ( provide_copying ) {
		provide_xsd_pitch_copying_information( xsd );
	}
}

/// @brief Reset the copying options.
/// @details This version resets pitch_copying_mode_, and calls the parent reset_copying_settings() function.
void
OmegaBundleParameter::reset_copying_settings() {
	Parameter::reset_copying_settings();
	pitch_copying_mode_ = false;
}

/////////////////// PROTECTED FUNCTIONS: ///////////////////

/// @brief Parse a setting for this parameter (e.g. r0="12.5").
/// @details Returns false if this parameter isn't provided.  Also parses copying.
/// @note May be overridden by derived classes.
bool
OmegaBundleParameter::parse_setting_only(
	utility::tag::TagCOP tag,
	bool const parse_copying_too
) {
	if ( parse_copying_too && ( tag->hasOption( parameter_name() + "_" + copy_suffix() ) || tag->hasOption( "pitch_from_helix" ) ) ) {
		runtime_assert_string_msg( !tag->hasOption( parameter_name() ), "Error in RealValuedParameter::parse_setting_only(): Options were provided for both copying parameter values and setting parameter values for the " + parameter_name() + " parameter." );
		if ( tag->hasOption( parameter_name() + "_" + copy_suffix() ) ) {
			set_copy_information( tag );
			set_copies_pitch(false);
		} else {
			set_pitch_copying_information( tag );
		}
		return true;
	}
	runtime_assert_string_msg( !tag->hasOption( parameter_name() + "_" + copy_suffix() ) && !tag->hasOption( "pitch_from_helix" ), "Error in RealValuedParameter::parse_setting_only(): Options were provided for copying parameter values for the " + parameter_name() + " parameter, and should not have been."  );
	if ( !tag->hasOption( parameter_name() ) ) return false;
	set_value( tag->getOption<core::Real>( parameter_name() ) );
	return true;
}

///////////// PRIVATE FUNCTIONS ////////////////////////

/// @brief Given another parameter of the same type, copy its pitch.  This does *not* set value_set_ to true.
/// @details Performs type checking in debug mode.
/// @returns Returns TRUE for failure, FALSE for success.
bool
OmegaBundleParameter::copy_pitch_from_parameter(
	core::conformation::parametric::ParameterCOP /*other_parameter*/,
	core::conformation::parametric::ParametersCOP other_parameter_collection,
	core::conformation::parametric::ParametersCOP this_parameter_collection
) {
	core::Real const other_r0( utility::pointer::static_pointer_cast< core::conformation::parametric::RealValuedParameter const >( other_parameter_collection->parameter_cop( static_cast<core::Size>( protocols::helical_bundle::BPC_r0 ) ) )->value() );
	core::Real const other_omega0( utility::pointer::static_pointer_cast< core::conformation::parametric::RealValuedParameter const >( other_parameter_collection->parameter_cop( static_cast<core::Size>( protocols::helical_bundle::BPC_omega0 ) ) )->value() );
	core::Real const other_z1( utility::pointer::static_pointer_cast< core::conformation::parametric::RealValuedParameter const >( other_parameter_collection->parameter_cop( static_cast<core::Size>( protocols::helical_bundle::BPC_z1 ) ) )->value() );
	core::Real const other_sinalpha( other_r0 * other_omega0 / other_z1 );
	if ( other_sinalpha > 1 || other_sinalpha < -1 ) {
		if ( TR.visible() ) TR << "Failed to copy pitch angle.  Current parameters do not generate a sensible pitch angle for a copied helix." << std::endl;
		return true;
	}

	//If we've got a good pitch angle, then continue:
	core::Real const other_alpha( asin(other_sinalpha) );

	core::Real const this_r0( utility::pointer::static_pointer_cast< core::conformation::parametric::RealValuedParameter const >( this_parameter_collection->parameter_cop( static_cast<core::Size>( protocols::helical_bundle::BPC_r0 ) ) )->value() ); //Already set, if sampled or if copied.
	core::Real const this_z1( utility::pointer::static_pointer_cast< core::conformation::parametric::RealValuedParameter const >( this_parameter_collection->parameter_cop( static_cast<core::Size>( protocols::helical_bundle::BPC_z1 ) ) )->value() ); //Cannot be sampled or copied.

	/********************
	We know: tan(alpha)=2*PI*R0/P, where alpha is the pitch angle, P is the pitch (rise per turn about major helix), and R0 is the major radius.
	sin(alpha)=R0*omega0/z1
	We want: P' = P
	2*PI*RO'/tan(alpha') = 2*PI*R0/tan(alpha)
	tan(alpha) = R0/R0'*tan(alpha')
	alpha = atan(R0/R0'*tan(alpha')
	R0*omega0/z1 = sin(atan(R0/R0'*tan(asin(R0'*omega0'/z1'))))
	omega0 = z1/R0*sin(atan(R0/R0'*tan(asin(R0'*omega0'/z1'))))
	********************/

	set_value_sneakily( this_z1/this_r0 * sin(atan(this_r0/other_r0*tan(other_alpha) ) ) ); //The "sneaky" version does not set value_set=true.
	return false;
}

/// @brief Parse copy information from the tag, and set it.
void
OmegaBundleParameter::set_pitch_copying_information(
	utility::tag::TagCOP tag
) {
	static const std::string copy_tag( "pitch_from_helix" );
	runtime_assert( tag->hasOption( copy_tag ) ); //Should already have been confirmed.
	reset_sampling_and_perturbation_settings();
	reset_value_settings();
	set_copy_from_parameters_index( tag->getOption< core::Size >( copy_tag ) );
	set_copies_pitch(true);
}

/// @brief Return the XSD information for copying this parameter from another.
void
OmegaBundleParameter::provide_xsd_pitch_copying_information(
	utility::tag::AttributeList & xsd
) const {
	using namespace utility::tag;
	xsd + XMLSchemaAttribute::attribute_w_default( "pitch_from_helix", xsct_non_negative_integer, "The index of the parametric object (e.g. the helix, in the case of a helical bundle) from which pitch value should be copied in order to set omega0, the twist per residue.  An alternative to \"omega0_copies_helix\".", "0" );
}


} //protocols
} //helical_bundle
} //parameters

#ifdef SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::helical_bundle::parameters::OmegaBundleParameter::save( Archive & arc ) const {
	arc( cereal::base_class< core::conformation::parametric::RealValuedParameter >( this ) );
	arc( CEREAL_NVP( pitch_copying_mode_ ) ); //bool
	//TODO -- archive private member variables here.
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::helical_bundle::parameters::OmegaBundleParameter::load( Archive & arc ) {
	arc( cereal::base_class< core::conformation::parametric::RealValuedParameter >( this ) );
	arc( pitch_copying_mode_ ); //bool
	//TODO -- de-archive private member variables here.
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::helical_bundle::parameters::OmegaBundleParameter );
CEREAL_REGISTER_TYPE( protocols::helical_bundle::parameters::OmegaBundleParameter )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_helical_bundle_parameters_OmegaBundleParameter )
#endif // SERIALIZATION
