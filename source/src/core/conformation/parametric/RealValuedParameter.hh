// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/conformation/parameters/RealValuedParameter.hh
/// @brief  Prototypes and method declarations for the RealValuedParameter class, a class for holding a single real-valued parameter for parametric backbone generation.
/// @author Vikram K. Mulligan (vmullig@uw.edu)


#ifndef INCLUDED_core_conformation_parametric_RealValuedParameter_hh
#define INCLUDED_core_conformation_parametric_RealValuedParameter_hh


// Unit headers
#include <core/conformation/parametric/Parameter.hh>
#include <core/conformation/parametric/RealValuedParameter.fwd.hh>

// Package headers
#include <core/conformation/Residue.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

// Numeric headers
#include <numeric/constants.hh>

// C++ headers


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace conformation {
namespace parametric {

/// @brief The types of perturbations.
/// @details If you add to this list, be sure to update the perturbation_type_string_from_enum() function.
enum RealPerturbationType {
	RPT_gaussian=1, //Keep first
	RPT_uniform,
	RPT_unknown_type, //Keep second-to-last
	RPT_end_of_list=RPT_unknown_type //Keep last
};

/// @brief  RealValuedParameter class, used to store a single real-valued parameter for parametric backbone generation.
///
class RealValuedParameter : public Parameter
{
public:

	/// @brief constructors
	///
	RealValuedParameter();

	RealValuedParameter( RealValuedParameter const & src );

	~RealValuedParameter() override;

	/// @brief Make a copy of this object ( allocate actual memory for it )
	///
	ParameterOP clone() const override;

public: //Getters

	/// @brief Get the value of this parameter.
	inline core::Real const & value() const { return value_; }

	/// @brief Get the default value of this parameter.
	inline core::Real const & default_value() const { return default_value_; }

	/// @brief Get the minimum value of this parameter.
	inline core::Real const & value_min() const { return value_min_; }

	/// @brief Get the maximum value of this parameter.
	inline core::Real const & value_max() const { return value_max_; }

	/// @brief Get the samples.
	inline core::Size value_samples() const { return value_samples_; }

	/// @brief Get the perturbation type.
	inline RealPerturbationType perturbation_type() const { return perturbation_type_; }

	/// @brief Get the perturbation magnitude of this parameter.
	inline core::Real const & perturbation_magnitude() const { return perturbation_magnitude_; }

public: //Setters

	/// @brief Set whether this parameter expects that its input is an angle in degrees.
	/// @details Must be set before any values are set.
	void set_input_is_angle_in_degrees( bool const setting );

	/// @brief Set the value of this parameter.
	/// @details If ignore_use_degrees is true, this will just set the value to the input value.  Otherwise, it convert from degrees
	/// to radians if input_is_angle_in_degrees_ is true.
	/// @note This makes the value_set_ boolean true.
	void set_value( core::Real const &value_in, bool const ignore_use_degrees=false );

	/// @brief Sets the initial value of this parameter, keeping the value_set_ boolean false.
	/// @details Also sets default_value_.
	void set_default_value( core::Real const &value_in );

	/// @brief Sets the parameter type.
	/// @details Override is limited to real types.
	void set_parameter_type( ParameterType const type_in ) override;

	/// @brief Set the options for sampling this parameter.
	void set_sampling_options( core::Real const &value_min_in, core::Real const &value_max_in, core::Size const samples_in );

	/// @brief Set the perturbation type, by string.
	/// @details Throws an error if not a valid type.
	void set_perturbation_type( std::string const &pert_type_in );

	/// @brief Set the perturbation type, by enum.
	/// @details Throws an error if invalid enum.
	void set_perturbation_type( RealPerturbationType const pert_type_in );

	/// @brief Set the options for perturbing this parameter (with perturbation type set by string).
	void set_perturbation_options( core::Real const &pert_magnitude_in, std::string const &pert_type_in );

	/// @brief Set the options for perturbing this parameter (with perturbation type set by enum).
	void set_perturbation_options( core::Real const &pert_magnitude_in, RealPerturbationType const pert_type_in );

	/// @brief Given another parameter of the same type, copy its value.  This does *not* set value_set_ to true.
	/// @details Performs type checking in debug mode.
	/// @returns Returns TRUE for failure, FALSE for success.
	bool copy_value_from_parameter( ParameterCOP other_parameter, ParametersCOP other_parameter_collection, ParametersCOP this_parameter_collection ) override;

public: //Parse functions

	/// @brief Given a tag, parse out the setting for this parameter.
	/// @details If parse_setting is true, this tries to parse the value for this parameter.  If parse_grid_sampling is true, this tries to parse
	/// a range of values to sample, and a number of samples.  If parse_perturbation is true, this parses options for perturbing the value of the
	/// parameter.  Note that, if multiple options are true, grid sampling or perturbation options are given priority over a flat setting (and an error
	/// is thrown if more than one of these is provided).
	/// @note Must be implemented by derived classes.
	void parse_setting( utility::tag::TagCOP tag, bool const parse_setting, bool const parse_grid_sampling, bool const parse_perturbation, bool const parse_copying ) override;

	/// @brief Return the XSD information for this parameter.
	/// @note Must be implemented by derived classes.
	void provide_xsd_information( utility::tag::AttributeList & xsd, bool const provide_setting, bool const provide_copying, bool const provide_grid_sampling, bool const provide_perturbation ) const override;

	/// @brief Reset the sampling and perturbation options before storing this Parameter object in a parametric Conformation.
	/// @details Pure virtual.  Must be implemented by derived classes.
	void reset_sampling_and_perturbation_settings() override;

	/// @brief Reset the value settings for this Parameter object.
	void reset_value_settings() override;

protected: //Functions that can be called by derived classes.

	/// @brief Set the value without setting value_was_set_ = true.
	void set_value_sneakily( core::Real const &value_in );

	/// @brief Parse copy information from the tag, and set it.
	void set_copy_information( utility::tag::TagCOP tag );

private: //Functions

	/// @brief Ensure that an angle value is in radians.
	/// @details  Checks the use_degrees_ boolean.  If true, converts degrees to radians; if false, returns input value.
	inline core::Real convert_angle( core::Real const &val ) const {
		if ( input_is_angle_in_degrees_ ) return ( val / 180.0 * numeric::constants::d::pi );
		return val; //Default case -- don't alter the value.
	}

	/// @brief Return the XSD information for grid-sampling this parameter.
	void provide_xsd_grid_sampling_information( utility::tag::AttributeList & xsd ) const;

	/// @brief Return the XSD information for perturbing this parameter.
	void provide_xsd_perturbation_information( utility::tag::AttributeList & xsd ) const;

	/// @brief Return the XSD information for setting this parameter.
	void provide_xsd_setting_information( utility::tag::AttributeList & xsd ) const;

	/// @brief Return the XSD information for copying this parameter from another.
	void provide_xsd_copying_information( utility::tag::AttributeList & xsd ) const;

public:

	/// @brief Given a perturbation type enum, get the string for that enum.
	/// @details Throws an error for an invalid enum.
	static std::string const & perturbation_type_string_from_enum( RealPerturbationType const type_enum );

	/// @brief Given a perturbation type string, get the enum for that string.
	/// @details Returns RPT_unknown_type if unknown string provided.
	static RealPerturbationType perturbation_type_enum_from_string( std::string const &type_string );

	/// @brief Has the perturbation information for this parameter been set?
	inline bool perturbation_set() const { return perturbation_set_; }

	/// @brief Has the sampling information for this parameter been set?
	inline bool sampling_set() const { return sampling_set_; }

	/// @brief Given the current value of this parameter and the current perturbation
	/// settings, return a perturbed value.
	/// @note Does *not* alter the value of this parameter.
	core::Real generate_perturbed_value( core::Real const &current_value ) const;

protected: //Functions that can be overridden by derived classes.

	/// @brief Parse a setting for this parameter (e.g. r0="12.5").
	/// @details Returns false if this parameter isn't provided.  Also parses copying.
	/// @note May be overridden by derived classes.
	virtual bool parse_setting_only( utility::tag::TagCOP tag, bool const parse_copying_too );

private:

	/// @brief If this is a parameter storing an angle, correct the value to be (-Pi:Pi].  Also, check that the value is positive if it is meant to be, or
	/// nonnegative if it is meant to be.
	/// @details Also checks sampling and perturbing options to ensure that these are reasonable.
	void correct_range();

	/// @brief Parse a sampling range for this parameter (e.g. r0_min="12.3" r0_max="14.3" r0_samples="5").
	/// @details Returns false if none of these parameters is provided.
	bool parse_grid_sampling_options( utility::tag::TagCOP tag );

	/// @brief Parse perturbation options for this parameter (e.g. r0_perturbation="5.0" r0_perturbation_type="gaussian").
	/// @details Returns false if the _perturbation option is not provided.
	bool parse_perturbation_options( utility::tag::TagCOP tag );

private:

	/********************************************************************************
	PRIVATE DATA
	*********************************************************************************/

	/// @brief The value of this parameter.
	core::Real value_;

	/// @brief The default value of this parameter.
	core::Real default_value_;

	/// @brief Has the default value been set?
	bool default_value_was_set_;

	/// @brief The min value to be sampled for this parameter.
	core::Real value_min_;

	/// @brief The max value to be sampled for this parameter.
	core::Real value_max_;

	/// @brief The number of samples for this parameter.
	core::Size value_samples_;

	/// @brief The magnitude of the perturbation of this parameter.
	core::Real perturbation_magnitude_;

	/// @brief The perturbation type for this parameter.
	RealPerturbationType perturbation_type_;

	/// @brief Was the sampling range set?
	bool sampling_set_;

	/// @brief Were perturbation options set?
	bool perturbation_set_;

	/// @brief Does this parameter expect that its input will be an angle in degrees?
	/// @details Will be converted to radians internally, in this case.
	bool input_is_angle_in_degrees_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; //class RealValuedParameter

} // namespace parametric
} // namespace conformation
} // namespace core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_conformation_parametric_RealValuedParameter )
#endif // SERIALIZATION


#endif
