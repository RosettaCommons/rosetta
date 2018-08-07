// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/conformation/parametric/SizeVectorValuedParameter.cc
/// @brief  A class for holding a single utility::vector1<core::Size>-valued parameter for parametric backbone generation.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Unit header
#include <core/conformation/parametric/SizeVectorValuedParameter.hh>

// Package headers

// Project headers

// Basic headers
#include <basic/Tracer.hh>

// Numeric headers

// Utility Headers
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>

// C++ headers
#include <sstream>

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

static basic::Tracer TR( "core.conformation.parametric.SizeVectorValuedParameter" );

/// @brief Constructor.
///
SizeVectorValuedParameter::SizeVectorValuedParameter() :
	Parameter(),
	values_( utility::vector1< core::Size >( 1, 0 ) ),
	default_value_(values_)
	//TODO -- initialize variables here.
{
	Parameter::set_parameter_type( PT_generic_whole_number_vector );
	set_can_be_set_sampled_perturbed_copied( true, false, false, false );
}

SizeVectorValuedParameter::SizeVectorValuedParameter( SizeVectorValuedParameter const & src ) :
	Parameter(src),
	values_(src.values_),
	default_value_(src.default_value_)
{
}

SizeVectorValuedParameter::~SizeVectorValuedParameter() {}

/// @brief Make a copy of this object ( allocate actual memory for it )
///
ParameterOP
SizeVectorValuedParameter::clone() const
{
	return ParameterOP( SizeVectorValuedParameterOP( new SizeVectorValuedParameter( *this ) ) );
}

/// @brief Set the value of this parameter.
void
SizeVectorValuedParameter::set_value(
	utility::vector1< core::Size > const & values_in
) {
	values_ = values_in;
	correct_range();
	set_value_was_set();
}

/// @brief Sets the parameter type.
/// @details Override is limited to nonnegative or positive integer types.
void
SizeVectorValuedParameter::set_parameter_type(
	ParameterType const type_in
) {
	runtime_assert( type_in == PT_generic_natural_number_vector || type_in == PT_generic_whole_number_vector );
	Parameter::set_parameter_type(type_in);
	correct_range();
}

/// @brief Sets the default value.
void
SizeVectorValuedParameter::set_default_value(
	utility::vector1< core::Size > const & values_in
) {
	runtime_assert_string_msg( !value_was_set(), "Error in SizeVectorValuedParameter::set_default_value(): The default value cannot be set after the value has been set!");
	default_value_ = values_in;
	values_ = values_in;
}

/// @brief Given another parameter of the same type, copy its value.  This does *not* set value_set_ to true.
/// @details Performs type checking in debug mode.
bool
SizeVectorValuedParameter::copy_value_from_parameter(
	ParameterCOP other_parameter,
	ParametersCOP /*other_parameter_collection*/,
	ParametersCOP /*this_parameter_collection*/
) {
#ifndef NDEBUG
	debug_assert( utility::pointer::dynamic_pointer_cast< SizeVectorValuedParameter const >(other_parameter) != nullptr );
#endif
	SizeVectorValuedParameterCOP other_parameter_cast( utility::pointer::static_pointer_cast< SizeVectorValuedParameter const >( other_parameter ) );
	values_ = other_parameter_cast->value();
	return false;
}

//////////////////////////// PARSE FUNCTIONS ////////////////////////////////////////////////////////

/// @brief Given a tag, parse out the setting for this parameter.
/// @details If parse_setting is true, this tries to parse the value for this parameter.  If parse_grid_sampling is true, this tries to parse
/// a range of values to sample, and a number of samples.  If parse_perturbation is true, this parses options for perturbing the value of the
/// parameter.  Note that, if multiple options are true, grid sampling or perturbation options are given priority over a flat setting (and an error
/// is thrown if more than one of these is provided).
/// @note A vector of integers must be provided as a whitespace-separated string.
void
SizeVectorValuedParameter::parse_setting(
	utility::tag::TagCOP tag,
	bool const parse_setting,
	bool const parse_grid_sampling,
	bool const parse_perturbation,
	bool const parse_copying
) {
	std::string const errmsg("Error in core::conformation::parametric::SizeValuedParameter::parse_setting(): ");
	runtime_assert_string_msg( parse_setting && !parse_grid_sampling && !parse_perturbation && !parse_copying, errmsg + "Grid sampling, perturbations, and parameter copying are not currently supported for integer vector-valued parameters." );

	if ( tag->hasOption(parameter_name()) ) {
		std::string const vectorstring( tag->getOption< std::string >(parameter_name()) );
		std::string const vectorstring2( utility::strip( vectorstring, " \t\n" ) ); //Remove trailing whitespace
		std::istringstream vectorstream( vectorstring2 );
		utility::vector1< core::Size > values_new;
		core::Size val(0);
		while ( !vectorstream.eof() ) {
			vectorstream >> val;
			runtime_assert_string_msg( !vectorstream.fail(), errmsg + "Could not parse integer vector.  The vector must be provided as a series of whitespace-separated integers." );
			values_new.push_back( val );
		}
		set_value( values_new );
	}
}

/// @brief Return the XSD information for this parameter.
/// @note Currently, SizeVector-valued parameters can't be sampled or perturbed, so this only provides setting information.
void
SizeVectorValuedParameter::provide_xsd_information(
	utility::tag::AttributeList & xsd,
	bool const provide_setting,
	bool const /*provide_copying*/,
	bool const /*provide_grid_sampling*/,
	bool const /*provide_perturbation*/
) const {
	using namespace utility::tag;

	if ( provide_setting ) {
		std::stringstream default_stream;
		for ( core::Size i(1), imax( default_value_.size() ); i<=imax; ++i ) {
			default_stream << default_value_[i];
			if ( i<imax ) default_stream << " ";
		}
		xsd + XMLSchemaAttribute::attribute_w_default( parameter_name(), parameter_type() == PT_generic_natural_number_vector ? xsct_positive_integer_wsslist : xsct_nnegative_int_wsslist, parameter_description(), default_stream.str() );
	}
}

/// @brief Reset the sampling and perturbation options before storing this Parameter object in a parametric Conformation.
/// @details Pure virtual.  Must be implemented by derived classes.
void
SizeVectorValuedParameter::reset_sampling_and_perturbation_settings() { /*GNDN -- Does nothing for this derived class, since there's nothing to do.*/ }

/// @brief Reset the value settings for this Parameter object.
void
SizeVectorValuedParameter::reset_value_settings() {
	set_value( default_value() );
	reset_value_was_set();
}

//////////////////////////// PRIVATE FUNCTIONS ////////////////////////////////////////////////////////

/// @brief Check that the value is positive if it is meant to be, or nonnegative if it is meant to be.
void
SizeVectorValuedParameter::correct_range() {
	static const std::string errmsg( "Error in core::conformation::parametric::SizeVectorValuedParameter::correct_range(): " );
	runtime_assert_string_msg( !values_.empty(), errmsg + "The stored vector is empty!" );
	if ( parameter_type() == PT_generic_natural_number_vector ) {
		for ( core::Size i(1), imax(values_.size()); i<=imax; ++i ) {
			runtime_assert_string_msg( values_[i] > 0, errmsg + "A natural number's value must be greater than zero!" );
		}
	}
}

} // namespace parametric
} // namespace conformation
} // namespace core


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::conformation::parametric::SizeVectorValuedParameter::save( Archive & arc ) const {
	arc( cereal::base_class< core::conformation::parametric::Parameter >( this ) );
	arc( CEREAL_NVP( values_ ) ); //utility::vector1<core::Size>
	arc( CEREAL_NVP( default_value_ ) ); //utility::vector1<core::Size>
	//TODO -- archive private member variables here.
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::conformation::parametric::SizeVectorValuedParameter::load( Archive & arc ) {
	arc( cereal::base_class< core::conformation::parametric::Parameter >( this ) );
	arc( values_ ); //utility::vector1<core::Size>
	arc( default_value_ ); //utility::vector1<core::Size>
	//TODO -- de-archive private member variables here.
}

SAVE_AND_LOAD_SERIALIZABLE( core::conformation::parametric::SizeVectorValuedParameter );
CEREAL_REGISTER_TYPE( core::conformation::parametric::SizeVectorValuedParameter )

CEREAL_REGISTER_DYNAMIC_INIT( core_conformation_parametric_SizeVectorValuedParameter )
#endif // SERIALIZATION
