// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/conformation/parametric/ParametrizationCalculator.cc
/// @brief  Function definitions for the ParametrizationCalculator class, a base class
/// from which classes that calculate particular parametrizations will be derived.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Unit header
#include <core/conformation/parametric/ParametrizationCalculator.hh>
#include <core/conformation/parametric/RealValuedParameter.hh>
#include <core/conformation/parametric/RealVectorValuedParameter.hh>
#include <core/conformation/parametric/SizeValuedParameter.hh>
#include <core/conformation/parametric/SizeVectorValuedParameter.hh>
#include <core/conformation/parametric/BooleanValuedParameter.hh>

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

static basic::Tracer TR( "core.conformation.parametric.ParametrizationCalculator" );

/// @brief Constructor.
///
ParametrizationCalculator::ParametrizationCalculator() :
	parameters_( new Parameters )
	//TODO -- initialize variables here.
{
}

/// @brief Constructor with input Parameters object.
///
ParametrizationCalculator::ParametrizationCalculator( ParametersOP parameters ) :
	parameters_( parameters )
	//TODO -- initialize variables here.
{
}

/// @brief Copy constructor.
/// @details Deep-copies the stored parameters.
ParametrizationCalculator::ParametrizationCalculator( ParametrizationCalculator const & src ) :
	utility::pointer::ReferenceCount(),
	utility::pointer::enable_shared_from_this< ParametrizationCalculator >(),
	parameters_( src.parameters_->clone() )
{
}

/// @brief Pure virtual destructor means don't instantiate this base class!.
/// @details Counter-intuitively, C++ requires pure virtual destructors to be implemented.
ParametrizationCalculator::~ParametrizationCalculator() {}

////////////////////// PUBLIC FUNCTIONS ///////////////////////////

/// @details Access a parameter.
///
ParameterOP
ParametrizationCalculator::parameter(
	core::Size const index
) {
	debug_assert( index > 0 );
	debug_assert( index <= parameters_->num_parameters() );
	return parameters_->parameter_op(index);
}

/// @brief Access a parameter.
/// @details Const-access
ParameterCOP
ParametrizationCalculator::parameter_cop(
	core::Size const index
) const {
	debug_assert( index > 0 );
	debug_assert( index <= parameters_->num_parameters() );
	return parameters_->parameter_cop(index);
}

/// @brief Accesses a real-valued parameter.
/// @details Returns nullptr if the parameter with the given index is not a RealValuedParameter.
RealValuedParameterOP
ParametrizationCalculator::real_parameter( core::Size const index ) {
	return utility::pointer::dynamic_pointer_cast<RealValuedParameter>( parameter(index) );
}

/// @brief Accesses a vector-valued parameter.
/// @details Returns nullptr if the parameter with the given index is not a RealVectorValuedParameter.
RealVectorValuedParameterOP
ParametrizationCalculator::realvector_parameter( core::Size const index ) {
	return utility::pointer::dynamic_pointer_cast<RealVectorValuedParameter>( parameter(index) );
}

/// @brief Accesses a size-valued parameter.
/// @details Returns nullptr if the parameter with the given index is not a SizeValuedParameter.
SizeValuedParameterOP
ParametrizationCalculator::size_parameter( core::Size const index ) {
	return utility::pointer::dynamic_pointer_cast<SizeValuedParameter>( parameter(index) );
}

/// @brief Accesses a sizevector-valued parameter.
/// @details Returns nullptr if the parameter with the given index is not a SizeVectorValuedParameter.
SizeVectorValuedParameterOP
ParametrizationCalculator::sizevector_parameter( core::Size const index ) {
	return utility::pointer::dynamic_pointer_cast<SizeVectorValuedParameter>( parameter(index) );
}

/// @brief Accesses a boolean-valued parameter.
/// @details Returns nullptr if the parameter with the given index is not a BooleanValuedParameter.
BooleanValuedParameterOP
ParametrizationCalculator::boolean_parameter( core::Size const index ) {
	return utility::pointer::dynamic_pointer_cast<BooleanValuedParameter>( parameter(index) );
}

/// @brief Accesses a real-valued parameter (const-access).
/// @details Returns nullptr if the parameter with the given index is not a RealValuedParameter.
RealValuedParameterCOP
ParametrizationCalculator::real_parameter_cop( core::Size const index ) const {
	return utility::pointer::dynamic_pointer_cast<RealValuedParameter const>( parameter_cop(index) );
}

/// @brief Accesses a realvector-valued parameter (const-access).
/// @details Returns nullptr if the parameter with the given index is not a RealVectorValuedParameter.
RealVectorValuedParameterCOP
ParametrizationCalculator::realvector_parameter_cop( core::Size const index ) const {
	return utility::pointer::dynamic_pointer_cast<RealVectorValuedParameter const>( parameter_cop(index) );
}


/// @brief Accesses a size-valued parameter (const-access).
/// @details Returns nullptr if the parameter with the given index is not a SizeValuedParameter.
SizeValuedParameterCOP
ParametrizationCalculator::size_parameter_cop( core::Size const index ) const {
	return utility::pointer::dynamic_pointer_cast<SizeValuedParameter const>( parameter_cop(index) );
}

/// @brief Accesses a sizevector-valued parameter (const-access).
/// @details Returns nullptr if the parameter with the given index is not a SizeVectorValuedParameter.
SizeVectorValuedParameterCOP
ParametrizationCalculator::sizevector_parameter_cop( core::Size const index ) const {
	return utility::pointer::dynamic_pointer_cast<SizeVectorValuedParameter const>( parameter_cop(index) );
}


/// @brief Accesses a boolean-valued parameter (const-access).
/// @details Returns nullptr if the parameter with the given index is not a BooleanValuedParameter.
BooleanValuedParameterCOP
ParametrizationCalculator::boolean_parameter_cop( core::Size const index ) const  {
	return utility::pointer::dynamic_pointer_cast<BooleanValuedParameter const>( parameter_cop(index) );
}


////////////////////// PROTECTED FUNCTIONS ///////////////////////////

/// @brief Add a parameter to this calculator, automatically determining the type.
///
void
ParametrizationCalculator::add_parameter(
	std::string const &parameter_name,
	ParameterType type,
	std::string const &description,
	std::string const &short_description,
	std::string const &units,
	ParameterizationCalculatorProperties const &properties
) {
	if ( type == PT_angle || type == PT_generic_real || type == PT_generic_positive_valued_real || type == PT_generic_nonnegative_valued_real ) {
		add_real_parameter(parameter_name, type, description, short_description, units, properties);
	} else if ( type == PT_angle_vector || type == PT_generic_real_vector || type == PT_generic_positive_valued_real_vector || type == PT_generic_nonnegative_valued_real_vector ) {
		add_realvector_parameter(parameter_name, type, description, short_description, units, properties);
	} else if ( type == PT_generic_natural_number || type == PT_generic_whole_number ) {
		add_size_parameter(parameter_name, type, description, short_description, units, properties);
	} else if ( type == PT_generic_natural_number_vector || type == PT_generic_whole_number_vector ) {
		add_sizevector_parameter(parameter_name, type, description, short_description, units, properties);
	} else if ( type == PT_boolean ) {
		add_boolean_parameter(parameter_name, type, description, short_description, units, properties);
	} else {
		utility_exit_with_message( "Error in core::conformation::parametric::ParametrizationCalculator::add_parameter(): Invalid parameter type!" );
	}
}

/// @brief Add a real-valued parameter to this calculator.
///
void
ParametrizationCalculator::add_real_parameter(
	std::string const &parameter_name,
	ParameterType type,
	std::string const &description,
	std::string const &short_description,
	std::string const &units,
	ParameterizationCalculatorProperties const &properties
) {
	RealValuedParameterOP param(  new RealValuedParameter );
	param->set_parameter_name(parameter_name);
	param->set_parameter_description(description);
	param->set_short_parameter_description(short_description);
	param->set_parameter_units(units);
	if ( type == PT_generic_positive_valued_real ) param->set_value(1.0);
	param->set_parameter_type(type); //Checks compatible type.
	param->set_can_be_set_sampled_perturbed_copied( properties.can_be_set, properties.can_be_copied, properties.can_be_sampled, properties.can_be_perturbed );
	param->set_global_for_parameters_set( properties.global_for_parameters_set );
	add_parameter( param );
}

/// @brief Add a vector-valued parameter to this calculator.
void
ParametrizationCalculator::add_realvector_parameter(
	std::string const &parameter_name,
	ParameterType type,
	std::string const &description,
	std::string const &short_description,
	std::string const &units,
	ParameterizationCalculatorProperties const &properties
) {
	RealVectorValuedParameterOP param(  new RealVectorValuedParameter );
	param->set_parameter_name(parameter_name);
	param->set_parameter_description(description);
	param->set_short_parameter_description(short_description);
	param->set_parameter_units(units);
	if ( type == PT_generic_positive_valued_real ) param->set_value(utility::vector1<core::Real>(1, 1.0));
	param->set_parameter_type(type); //Checks compatible type.
	param->set_can_be_set_sampled_perturbed_copied( properties.can_be_set, properties.can_be_copied, properties.can_be_sampled, properties.can_be_perturbed );
	param->set_global_for_parameters_set( properties.global_for_parameters_set );
	add_parameter( param );
}

/// @brief Add a core::Size-valued parameter to this calculator.
///
void
ParametrizationCalculator::add_size_parameter(
	std::string const &parameter_name,
	ParameterType type,
	std::string const &description,
	std::string const &short_description,
	std::string const &units,
	ParameterizationCalculatorProperties const &properties
) {
	SizeValuedParameterOP param( new SizeValuedParameter );
	param->set_parameter_name(parameter_name);
	param->set_parameter_description(description);
	param->set_short_parameter_description(short_description);
	param->set_parameter_units(units);
	if ( type == PT_generic_natural_number ) param->set_value(1);
	param->set_parameter_type(type); //Checks compatible type.
	param->set_can_be_set_sampled_perturbed_copied( properties.can_be_set, properties.can_be_copied, properties.can_be_sampled, properties.can_be_perturbed );
	param->set_global_for_parameters_set( properties.global_for_parameters_set );
	add_parameter( param );
}

/// @brief Add a integer vector-valued parameter to this calculator.
///
void
ParametrizationCalculator::add_sizevector_parameter(
	std::string const &parameter_name,
	ParameterType type,
	std::string const &description,
	std::string const &short_description,
	std::string const &units,
	ParameterizationCalculatorProperties const &properties
) {
	SizeVectorValuedParameterOP param( new SizeVectorValuedParameter );
	param->set_parameter_name(parameter_name);
	param->set_parameter_description(description);
	param->set_short_parameter_description(short_description);
	param->set_parameter_units(units);
	if ( type == PT_generic_natural_number_vector ) param->set_value(utility::vector1< core::Size >(1,1) );
	param->set_parameter_type(type); //Checks compatible type.
	param->set_can_be_set_sampled_perturbed_copied( properties.can_be_set, properties.can_be_copied, properties.can_be_sampled, properties.can_be_perturbed );
	param->set_global_for_parameters_set( properties.global_for_parameters_set );
	add_parameter( param );
}

/// @brief Add a Boolean-valued parameter to this calculator.
///
void
ParametrizationCalculator::add_boolean_parameter(
	std::string const &parameter_name,
	ParameterType type,
	std::string const &description,
	std::string const &short_description,
	std::string const &units,
	ParameterizationCalculatorProperties const &properties
) {
	BooleanValuedParameterOP param( new BooleanValuedParameter );
	param->set_parameter_name(parameter_name);
	param->set_parameter_description(description);
	param->set_short_parameter_description(short_description);
	param->set_parameter_units(units);
	param->set_parameter_type(type); //Checks compatible type.
	param->set_can_be_set_sampled_perturbed_copied( properties.can_be_set, properties.can_be_copied, properties.can_be_sampled, properties.can_be_perturbed );
	param->set_global_for_parameters_set( properties.global_for_parameters_set );
	add_parameter( param );
}

/// @brief Add a custom parameter
void
ParametrizationCalculator::add_custom_parameter(
	std::string const &parameter_name,
	ParameterType const type,
	std::string const &description,
	std::string const &short_description,
	std::string const &units,
	ParameterizationCalculatorProperties const &properties,
	ParameterOP parameter_in
) {
	parameter_in->set_parameter_name(parameter_name);
	parameter_in->set_parameter_description(description);
	parameter_in->set_short_parameter_description(short_description);
	parameter_in->set_parameter_units(units);
	parameter_in->set_parameter_type(type); //Checks compatible type.
	parameter_in->set_can_be_set_sampled_perturbed_copied( properties.can_be_set, properties.can_be_copied, properties.can_be_sampled, properties.can_be_perturbed );
	parameter_in->set_global_for_parameters_set( properties.global_for_parameters_set );
	add_parameter( parameter_in );
}

////////////////////// PRIVATE FUNCTIONS ///////////////////////////

/// @brief Add a parameter to this calculator.
///
void
ParametrizationCalculator::add_parameter(
	ParameterOP parameter_in
) {
	parameters_->add_parameter( parameter_in );
}

} // namespace parametric
} // namespace conformation
} // namespace core


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::conformation::parametric::ParametrizationCalculator::save( Archive & arc ) const {
	arc( CEREAL_NVP( parameters_ ) );
	//TODO -- archive private member variables here.
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::conformation::parametric::ParametrizationCalculator::load( Archive & arc ) {
	arc( parameters_ );
	//TODO -- de-archive private member variables here.
}

SAVE_AND_LOAD_SERIALIZABLE( core::conformation::parametric::ParametrizationCalculator );
CEREAL_REGISTER_TYPE( core::conformation::parametric::ParametrizationCalculator )

CEREAL_REGISTER_DYNAMIC_INIT( core_conformation_parametric_ParametrizationCalculator )
#endif // SERIALIZATION
