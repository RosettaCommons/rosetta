// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/conformation/parametric/RealValuedParameter.cc
/// @brief  A class for holding a single real-valued parameter for parametric backbone generation.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Unit header
#include <core/conformation/parametric/RealValuedParameter.hh>

// Package headers

// Project headers

// Basic headers
#include <basic/Tracer.hh>

// Numeric headers
#include <numeric/angle.functions.hh>
#include <numeric/random/random.hh>
#include <numeric/random/random.functions.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>
#include <utility/string_util.hh>

// Cereal headers
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace conformation {
namespace parametric {

static basic::Tracer TR( "core.conformation.parametric.RealValuedParameter" );

/// @brief Constructor.
///
RealValuedParameter::RealValuedParameter() :
	Parameter(),
	value_(0.0),
	default_value_(0.0),
	default_value_was_set_(false),
	value_min_(0.0),
	value_max_(0.0),
	value_samples_(1),
	perturbation_magnitude_(0.0),
	perturbation_type_(RPT_gaussian),
	sampling_set_(false),
	perturbation_set_(false),
	input_is_angle_in_degrees_(false)
	//TODO -- initialize variables here.
{
	Parameter::set_parameter_type( PT_generic_real );
	set_can_be_set_sampled_perturbed_copied( true, true, true, true );
}

RealValuedParameter::RealValuedParameter( RealValuedParameter const & src ) :
	Parameter(src),
	value_(src.value_),
	default_value_(src.default_value_),
	default_value_was_set_(src.default_value_was_set_),
	value_min_(src.value_min_),
	value_max_(src.value_max_),
	value_samples_(src.value_samples_),
	perturbation_magnitude_(src.perturbation_magnitude_),
	perturbation_type_(src.perturbation_type_),
	sampling_set_(src.sampling_set_),
	perturbation_set_(src.perturbation_set_),
	input_is_angle_in_degrees_(src.input_is_angle_in_degrees_)
{
}

RealValuedParameter::~RealValuedParameter() {}


/// @brief Make a copy of this object ( allocate actual memory for it )
///
ParameterOP
RealValuedParameter::clone() const
{
	return ParameterOP( RealValuedParameterOP( new RealValuedParameter( *this ) ) );
}

/// @brief Set whether this parameter expects that its input is an angle in degrees.
/// @details Must be set before any values are set.
void
RealValuedParameter::set_input_is_angle_in_degrees(
	bool const setting
) {
	static std::string const errmsg( "Error in RealValuedParameter::set_input_is_angle_in_degrees(): " );
	runtime_assert_string_msg( !value_was_set(), errmsg + "The input format cannot be changed after the value is set." );
	runtime_assert_string_msg( !sampling_set_, errmsg + "The input format cannot be changed after the sampling parameters are set."  );
	runtime_assert_string_msg( !perturbation_set_, errmsg + "The input format cannot be changed after the perturbaton parameters are set."  );
	runtime_assert_string_msg( !default_value_was_set_, errmsg + "The input format cannot be changed after the default value has been set."  );
	input_is_angle_in_degrees_ = setting;
}

/// @brief Set the value of this parameter.
/// @note This makes the value_set_ boolean true.
void
RealValuedParameter::set_value(
	core::Real const &value_in,
	bool const ignore_use_degrees/*=false*/
) {
	if ( ignore_use_degrees ) {
		value_ = value_in;
	} else {
		value_ = convert_angle( value_in );
	}
	correct_range();
	set_value_was_set();
}

/// @brief Sets the initial value of this parameter, keeping the value_set_ boolean false.
/// @details Also sets default_value_.
void
RealValuedParameter::set_default_value(
	core::Real const &value_in
) {
	runtime_assert_string_msg(!value_was_set(), "Error in RealValuedParameter::set_default_value(): The default value cannot be set after the value is set!");
	value_ = convert_angle( value_in );
	correct_range();
	default_value_ = value_;
	default_value_was_set_ = true;
}

/// @brief Sets the parameter type.
/// @details Override is limited to real types.
void
RealValuedParameter::set_parameter_type(
	ParameterType const type_in
) {
	runtime_assert( type_in == PT_angle || type_in == PT_generic_real || type_in == PT_generic_nonnegative_valued_real || type_in == PT_generic_positive_valued_real );
	Parameter::set_parameter_type(type_in);
	correct_range();
}

/// @brief Set the options for sampling this parameter.
void
RealValuedParameter::set_sampling_options(
	core::Real const &value_min_in,
	core::Real const &value_max_in,
	core::Size const samples_in
) {
	reset_sampling_and_perturbation_settings();
	reset_copying_settings();
	reset_value_settings();
	value_min_ = convert_angle( value_min_in );
	value_max_ = convert_angle( value_max_in );
	value_samples_ = samples_in;
	correct_range();
	sampling_set_ = true;
}

/// @brief Set the perturbation type, by string.
/// @details Throws an error if not a valid type.
void
RealValuedParameter::set_perturbation_type(
	std::string const &pert_type_in
) {
	set_perturbation_type( perturbation_type_enum_from_string( pert_type_in ) );
}

/// @brief Set the perturbation type, by enum.
/// @details Throws an error if invalid enum.
void
RealValuedParameter::set_perturbation_type(
	RealPerturbationType const pert_type_in
) {
	runtime_assert_string_msg( pert_type_in > 0 && pert_type_in < RPT_unknown_type,
		"Error in core::conformation::parametric::RealValuedParameter::set_perturbation_type(): Unknown perturbation type specified.");
	perturbation_type_ = pert_type_in;
	// Note: it is possible that no perturbation was set if we're only setting the perturbation type.  So
	// we DON'T set perturbation_set_ = true.
}

/// @brief Set the options for perturbing this parameter (with perturbation type set by string).
void
RealValuedParameter::set_perturbation_options(
	core::Real const &pert_magnitude_in,
	std::string const &pert_type_in
) {
	set_perturbation_options( pert_magnitude_in, perturbation_type_enum_from_string( pert_type_in ) );
}

/// @brief Set the options for perturbing this parameter (with perturbation type set by enum).
void
RealValuedParameter::set_perturbation_options(
	core::Real const &pert_magnitude_in,
	RealPerturbationType const pert_type_in
) {
	runtime_assert_string_msg( pert_type_in > 0 && pert_type_in < RPT_unknown_type,
		"Error in core::conformation::parametric::RealValuedParameter::set_perturbation_options(): Unknown perturbation type specified.");
	reset_sampling_and_perturbation_settings();
	reset_copying_settings();
	reset_value_settings();
	perturbation_magnitude_ = convert_angle( pert_magnitude_in );
	perturbation_type_ = pert_type_in;
	correct_range();
	perturbation_set_ = true;
}

/// @brief Given another parameter of the same type, copy its value.  This does *not* set value_set_ to true.
/// @details Performs type checking in debug mode.
/// @returns Returns TRUE for failure, FALSE for success.
bool
RealValuedParameter::copy_value_from_parameter(
	ParameterCOP other_parameter,
	ParametersCOP /*other_parameter_collection*/,
	ParametersCOP /*this_parameter_collection*/
) {
#ifndef NDEBUG
	debug_assert( utility::pointer::dynamic_pointer_cast< RealValuedParameter const >(other_parameter) != nullptr );
#endif
	RealValuedParameterCOP other_parameter_cast( utility::pointer::static_pointer_cast< RealValuedParameter const >( other_parameter ) );
	value_ = other_parameter_cast->value();
	return false;
}

//////////////// PARSE FUNCTIONS ///////////////////

/// @brief Given a tag, parse out the setting for this parameter.
/// @details If parse_setting is true, this tries to parse the value for this parameter.  If parse_grid_sampling is true, this tries to parse
/// a range of values to sample, and a number of samples.  If parse_perturbation is true, this parses options for perturbing the value of the
/// parameter.  Note that, if multiple options are true, grid sampling or perturbation options are given priority over a flat setting (and an error
/// is thrown if more than one of these is provided).
/// @note Must be implemented by derived classes.
void
RealValuedParameter::parse_setting(
	utility::tag::TagCOP tag,
	bool const parse_setting,
	bool const parse_grid_sampling,
	bool const parse_perturbation,
	bool const parse_copying
) {
	debug_assert( !( parse_perturbation && parse_grid_sampling ) ); //Cannot parse both perturbation and grid sampling options.
	debug_assert(parse_setting); //Should always be parsing setting, even if parsing sampling.

	if ( parse_grid_sampling ) {
		if ( !parse_grid_sampling_options(tag) ) parse_setting_only(tag, parse_copying);
		return;
	} else if ( parse_perturbation ) {
		if ( !parse_perturbation_options(tag) ) parse_setting_only(tag, parse_copying);
		return;
	}

	parse_setting_only(tag, parse_copying);
}

/// @brief Return the XSD information for this parameter.
/// @note Must be implemented by derived classes.
void
RealValuedParameter::provide_xsd_information(
	utility::tag::AttributeList & xsd,
	bool const provide_setting,
	bool const provide_copying,
	bool const provide_grid_sampling,
	bool const provide_perturbation
) const {
	debug_assert( !(provide_grid_sampling && provide_perturbation) ); //Can't provide both.
	debug_assert( provide_setting ); //Should always be true.  Silly to make this an option, I guess.  Meh.

	provide_xsd_setting_information( xsd );

	if ( provide_grid_sampling ) {
		provide_xsd_grid_sampling_information( xsd );
	}
	if ( provide_perturbation ) {
		provide_xsd_perturbation_information( xsd );
	}
	if ( provide_copying ) {
		provide_xsd_copying_information( xsd );
	}

}


/// @brief Reset the sampling and perturbation options before storing this Parameter object in a parametric Conformation.
/// @details Pure virtual.  Must be implemented by derived classes.
void
RealValuedParameter::reset_sampling_and_perturbation_settings() {
	value_min_ = 0;
	value_max_ = 0;
	value_samples_ = 1;
	perturbation_magnitude_ = 0;
	perturbation_type_ = RPT_gaussian;
	sampling_set_ = false;
	perturbation_set_ = false;
}

/// @brief Reset the value settings for this Parameter object.
void
RealValuedParameter::reset_value_settings() {
	set_value( default_value(), true );
	reset_value_was_set();
}

//////////////// PROTECTED FUNCTIONS ////////////////

/// @brief Set the value without setting value_was_set_ = true.
void
RealValuedParameter::set_value_sneakily(
	core::Real const &value_in
) {
	value_ = value_in;
}

/// @brief Parse copy information from the tag, and set it.
void
RealValuedParameter::set_copy_information(
	utility::tag::TagCOP tag
) {
	std::string const copy_tag( parameter_name() + "_" + copy_suffix() );
	runtime_assert( tag->hasOption( copy_tag ) ); //Should already have been confirmed.
	set_copy_from_parameters_index( tag->getOption< core::Size >( copy_tag ) );
}

//////////////// PRIVATE FUNCTIONS ////////////////

/// @brief Return the XSD information for grid-sampling this parameter.
void
RealValuedParameter::provide_xsd_grid_sampling_information(
	utility::tag::AttributeList & xsd
) const {
	using namespace utility::tag;
	xsd + XMLSchemaAttribute::attribute_w_default( parameter_name() + "_min" , xsct_real, "Minimum value of sampling range for " + parameter_name() + ".", "0" );
	xsd + XMLSchemaAttribute::attribute_w_default( parameter_name() + "_max" , xsct_real, "Maximum value of sampling range for " + parameter_name() + ".", "0" );
	xsd + XMLSchemaAttribute::attribute_w_default( parameter_name() + "_samples" , xsct_positive_integer, "Number of samples when sampling " + parameter_name() + ".  Must be greater than zero.", "1" );
}

/// @brief Return the XSD information for perturbing this parameter.
void
RealValuedParameter::provide_xsd_perturbation_information(
	utility::tag::AttributeList & xsd
) const {
	using namespace utility::tag;
	xsd + XMLSchemaAttribute::attribute_w_default( parameter_name() + "_perturbation" , xsct_real, "Perturbation magnitude for perturbing " + parameter_name() + ".", "0.0" );
	xsd + XMLSchemaAttribute::attribute_w_default( parameter_name() + "_perturbation_type" , xs_string, "Perturbation type for perturbing " + parameter_name() + ".  Can be \"gaussian\" or \"uniform\".", "gaussian" );
	//TODO: FIGURE OUT HOW TO ADD RESTRICTION THAT THE PERTURBATION TYPE CAN ONLY BE "gaussian" OR "uniform"
}

/// @brief Return the XSD information for setting this parameter.
void
RealValuedParameter::provide_xsd_setting_information(
	utility::tag::AttributeList & xsd
) const {
	using namespace utility::tag;
	xsd + XMLSchemaAttribute::attribute_w_default( parameter_name(), xsct_real, parameter_description(), std::to_string( default_value() ) );
}

/// @brief Return the XSD information for copying this parameter from another.
void
RealValuedParameter::provide_xsd_copying_information(
	utility::tag::AttributeList & xsd
) const {
	using namespace utility::tag;
	xsd + XMLSchemaAttribute::attribute_w_default( parameter_name() + "_" + copy_suffix(), xsct_non_negative_integer, "The index of the parametric object (e.g. the helix, in the case of a helical bundle) from which the value for " + parameter_name() + " should be copied.", "0" );
}

///////////////// THE FOLLOWING TWO FUNCTIONS ARE PUBLIC AND STATIC: /////////////////

/// @brief Given a perturbation type enum, get the string for that enum.
/// @details Throws an error for an invalid enum.
std::string const &
RealValuedParameter::perturbation_type_string_from_enum(
	RealPerturbationType const type_enum
) {
	static utility::vector1< std::string > const types { //Must match enum in RealValuedParameter.hh.
		"gaussian",
		"uniform"
		};
	runtime_assert(type_enum > 0 && type_enum < RPT_end_of_list);
	return types[type_enum];
}

/// @brief Given a perturbation type string, get the enum for that string.
/// @details Returns RPT_unknown_type if unknown string provided.
RealPerturbationType
RealValuedParameter::perturbation_type_enum_from_string(
	std::string const &type_string
) {
	for ( core::Size i(1); static_cast<RealPerturbationType>(i)<RPT_end_of_list; ++i ) {
		if ( !type_string.compare( perturbation_type_string_from_enum(static_cast<RealPerturbationType>(i)) ) ) return static_cast<RealPerturbationType>(i);
	}
	return RPT_unknown_type;
}

/// @brief Given the current value of this parameter and the current perturbation
/// settings, return a perturbed value.
/// @note Does *not* alter the value of this parameter.
core::Real
RealValuedParameter::generate_perturbed_value( core::Real const &current_value ) const {
	runtime_assert_string_msg( perturbation_type_ == RPT_uniform || perturbation_type_ == RPT_gaussian, "Error in RealValuedParameter::generate_perturbed_value(): An unknown perturbation type was set." );
	core::Real pert;
	if ( perturbation_type_ == RPT_uniform ) {
		pert = (numeric::random::uniform() * 2.0 - 1.0) * perturbation_magnitude_ + current_value;
	} else { // RPT_gaussian case
		pert = numeric::random::gaussian() * perturbation_magnitude_ + current_value;
	}
	if ( parameter_type() == PT_angle ) {
		pert = numeric::principal_angle_radians( pert );
	} else {
		if ( pert < 0 ) {
			if ( parameter_type() == PT_generic_nonnegative_valued_real ) pert = 0;
			else if ( parameter_type() == PT_generic_positive_valued_real ) pert = 1e-12;
		}
	}

	return pert;
}

/////////////////// PROTECTED FUNCTIONS: ///////////////////

/// @brief Parse a setting for this parameter (e.g. r0="12.5").
/// @details Returns false if this parameter isn't provided.  Also parses copying.
/// @note May be overridden by derived classes.
bool
RealValuedParameter::parse_setting_only(
	utility::tag::TagCOP tag,
	bool const parse_copying_too
) {
	if ( parse_copying_too && tag->hasOption( parameter_name() + "_" + copy_suffix() ) ) {
		runtime_assert_string_msg( !tag->hasOption( parameter_name() ), "Error in RealValuedParameter::parse_setting_only(): Options were provided for both copying parameter values and setting parameter values for the " + parameter_name() + " parameter." );
		set_copy_information( tag );
		return true;
	}
	runtime_assert_string_msg( !tag->hasOption( parameter_name() + "_" + copy_suffix() ), "Error in RealValuedParameter::parse_setting_only(): Options were provided for copying parameter values for the " + parameter_name() + " parameter, and should not have been."  );
	if ( !tag->hasOption( parameter_name() ) ) return false;
	reset_sampling_and_perturbation_settings();
	reset_copying_settings();
	reset_value_settings();
	set_value( tag->getOption<core::Real>( parameter_name() ) );
	return true;
}

///////////////// MORE PRIVATE FUNCTIONS: /////////////////

/// @brief If this is a parameter storing an angle, correct the value to be (-Pi:Pi].  Also, check that the value is positive if it is meant to be, or
/// nonnegative if it is meant to be.
/// @details Also checks sampling and perturbing options to ensure that these are reasonable.
void
RealValuedParameter::correct_range() {
	static const std::string errmsg( "Error in core::conformation::parametric::RealValuedParameter::correct_range(): " );
	if ( parameter_type() == PT_angle ) {
		value_ = numeric::principal_angle_radians(value_);
		//value_min_ = numeric::principal_angle_radians(value_min_);
		//value_max_ = numeric::principal_angle_radians(value_max_);
	} else if ( parameter_type() == PT_generic_nonnegative_valued_real ) {
		runtime_assert_string_msg( value_ >= 0, errmsg + "The stored value is less than zero!" );
		runtime_assert_string_msg( value_min_ >= 0 && value_max_ >= 0, errmsg + "The stored value sampling range includes values less than zero!" );
		if ( value_min_ > value_max_ ) std::swap( value_min_, value_max_ );
	} else if ( parameter_type() == PT_generic_positive_valued_real ) {
		runtime_assert_string_msg( value_ > 0, errmsg + "The stored value is not greater than zero!" );
		runtime_assert_string_msg( value_min_ > 0 && value_max_ > 0, errmsg + "The stored value sampling range includes values that are not greater than zero!" );
		if ( value_min_ > value_max_ ) std::swap( value_min_, value_max_ );
	}
	runtime_assert_string_msg( perturbation_magnitude_ >= 0, errmsg + "The perturbation magnitude must be greater than or equal to zero." );
	runtime_assert_string_msg( perturbation_type_ > 0 && perturbation_type_ < RPT_end_of_list, errmsg + "The perturbation type is unknown." );
}

/// @brief Parse a sampling range for this parameter (e.g. r0_min="12.3" r0_max="14.3" r0_samples="5").
/// @details Returns false if none of these parameters is provided.
bool
RealValuedParameter::parse_grid_sampling_options(
	utility::tag::TagCOP tag
) {
	static std::string const errmsg( "Error in core::conformation::parametric::RealValuedParameter::parse_grid_sampling(): " );

	bool const hasmin( tag->hasOption( parameter_name() + "_min" ) );
	bool const hasmax( tag->hasOption( parameter_name() + "_max" ) );
	bool const hassamples( tag->hasOption( parameter_name() + "_samples" ) );

	if ( !hasmin && !hasmax && !hassamples ) {
		return false; //No options specified for sampling.
	}

	runtime_assert_string_msg( !tag->hasOption( parameter_name() ), errmsg + "Sampling options were specified alongside an absolute value for parameter " + parameter_name() + "." ); //Ought not to have an option for the absolute value if we are sampling.
	runtime_assert_string_msg( !tag->hasOption( parameter_name() + "_perturbation" ) && !tag->hasOption( parameter_name() + "_perturbation_type" ), errmsg + "Sampling options were specified alongside perturbation options for parameter " + parameter_name() + "." ); //Ought not to have an option for perturbing we are sampling.
	runtime_assert_string_msg( hasmin && hasmax && hassamples, errmsg + "All of \"" + parameter_name() + "_min\", \"" + parameter_name() + "_max\", and \"" + parameter_name() + "_samples\" must be specified for sampling." );

	set_sampling_options( tag->getOption<core::Real>(parameter_name() + "_min"), tag->getOption<core::Real>(parameter_name() + "_max"), tag->getOption<core::Size>(parameter_name() + "_samples") );
	return true;
}

/// @brief Parse perturbation options for this parameter (e.g. r0_perturbation="5.0" r0_perturbation_type="gaussian").
/// @details Returns false if the _perturbation option is not provided.
bool
RealValuedParameter::parse_perturbation_options(
	utility::tag::TagCOP tag
) {
	static std::string const errmsg( "Error in core::conformation::parametric::RealValuedParameter::parse_perturbation(): " );

	bool const haspert( tag->hasOption( parameter_name() + "_perturbation" ) );
	bool const hastype( tag->hasOption( parameter_name() + "_perturbation_type" ) );

	if ( !haspert && !hastype ) {
		return false; //No options specified for perturbing.
	}

	runtime_assert_string_msg( !tag->hasOption( parameter_name() ), errmsg + "Perturbation options were specified alongside an absolute value for parameter " + parameter_name() + "." ); //Ought not to have an option for the absolute value if we are perturbing.
	runtime_assert_string_msg( !tag->hasOption( parameter_name() + "_min" ) && !tag->hasOption( parameter_name() + "_max") && !tag->hasOption( parameter_name() + "_samples" ), errmsg + "Sampling options were specified alongside perturbation options for parameter " + parameter_name() + "." ); //Ought not to have an option for sampling we are perturbing.
	runtime_assert_string_msg( haspert, errmsg + "A perturbation type was specified for parameter " + parameter_name() + ", but no perturbation magnitude was provided." );

	set_perturbation_options( tag->getOption<core::Real>(parameter_name() + "_perturbation"), tag->getOption<std::string>(parameter_name() + "_perturbation_type", perturbation_type_string_from_enum( perturbation_type_ ) ) );
	return true;
}

} // namespace parametric
} // namespace conformation
} // namespace core


#ifdef SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::conformation::parametric::RealValuedParameter::save( Archive & arc ) const {
	arc( cereal::base_class< core::conformation::parametric::Parameter >( this ) );
	arc( CEREAL_NVP( value_ ) ); //core::Real
	arc( CEREAL_NVP( default_value_ ) ); //core::Real
	arc( CEREAL_NVP( default_value_was_set_ ) ); //bool
	arc( CEREAL_NVP( value_min_ ) ); //core::Real
	arc( CEREAL_NVP( value_max_ ) ); //core::Real
	arc( CEREAL_NVP( value_samples_ ) ); //core::Size
	arc( CEREAL_NVP( perturbation_magnitude_ ) ); //core::Real
	arc( CEREAL_NVP( perturbation_type_ ) ); //RealPerturbationType enum
	arc( CEREAL_NVP( sampling_set_ ) ); //bool
	arc( CEREAL_NVP( perturbation_set_ ) ); //bool
	arc( CEREAL_NVP( input_is_angle_in_degrees_ ) ); //bool
	//TODO -- archive private member variables here.
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::conformation::parametric::RealValuedParameter::load( Archive & arc ) {
	arc( cereal::base_class< core::conformation::parametric::Parameter >( this ) );
	arc( value_ ); //core::Real
	arc( default_value_ ); //core::Real
	arc( default_value_was_set_ ); //bool
	arc( value_min_ ); //core::Real
	arc( value_max_ ); //core::Real
	arc( value_samples_ ); //core::Size
	arc( perturbation_magnitude_ ); //core::Real
	arc( perturbation_type_ ); //RealPerturbationType enum
	arc( sampling_set_ ); //bool
	arc( perturbation_set_ ); //bool
	arc( input_is_angle_in_degrees_ ); //bool
	//TODO -- de-archive private member variables here.
}

SAVE_AND_LOAD_SERIALIZABLE( core::conformation::parametric::RealValuedParameter );
CEREAL_REGISTER_TYPE( core::conformation::parametric::RealValuedParameter )

CEREAL_REGISTER_DYNAMIC_INIT( core_conformation_parametric_RealValuedParameter )
#endif // SERIALIZATION
