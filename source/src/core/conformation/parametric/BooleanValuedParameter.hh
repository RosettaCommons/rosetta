// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/conformation/parameters/BooleanValuedParameter.hh
/// @brief  Prototypes and method declarations for the BooleanValuedParameter class, a class for holding a single Boolean-valued parameter for parametric backbone generation.
/// @author Vikram K. Mulligan (vmullig@uw.edu)


#ifndef INCLUDED_core_conformation_parametric_BooleanValuedParameter_hh
#define INCLUDED_core_conformation_parametric_BooleanValuedParameter_hh


// Unit headers
#include <core/conformation/parametric/Parameter.hh>
#include <core/conformation/parametric/BooleanValuedParameter.fwd.hh>

// Package headers
#include <core/conformation/Residue.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.fwd.hh>

// Numeric headers

// C++ headers


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace conformation {
namespace parametric {

/// @brief  BooleanValuedParameter class, used to store a single Boolean-valued parameter for parametric backbone generation.
///
class BooleanValuedParameter : public Parameter
{
public:

	/// @brief constructors
	///
	BooleanValuedParameter();

	BooleanValuedParameter( BooleanValuedParameter const & src );

	~BooleanValuedParameter() override;

	/// @brief Make a copy of this object ( allocate actual memory for it )
	///
	ParameterOP clone() const override;

public: //Getters

	/// @brief Get the value of this parameter.
	inline bool value() const { return value_; }

	/// @brief Get the default value of this parameter.
	inline bool default_value() const { return default_value_; }

public: //Setters

	/// @brief Set the value of this parameter.
	void set_value( bool const value_in );

	/// @brief Sets the parameter type.
	/// @details Override is limited to Boolean types.
	void set_parameter_type( ParameterType const type_in ) override;

	/// @brief Set the default value of this parameter.
	void set_default_value( bool const value_in );

	/// @brief Given another parameter of the same type, copy its value.  This does *not* set value_set_ to true.
	/// @details Performs type checking in debug mode.
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
	/// @note Currently, boolean values can't be sampled or perturbed, so this only provides setting information.
	void provide_xsd_information( utility::tag::AttributeList & xsd, bool const provide_setting, bool const provide_copying, bool const provide_grid_sampling, bool const provide_perturbation ) const override;

	/// @brief Reset the sampling and perturbation options before storing this Parameter object in a parametric Conformation.
	/// @details Pure virtual.  Must be implemented by derived classes.
	void reset_sampling_and_perturbation_settings() override;

	/// @brief Reset the value settings for this Parameter object.
	void reset_value_settings() override;

private:

	/********************************************************************************
	PRIVATE DATA
	*********************************************************************************/

	/// @brief The value of this parameter.
	bool value_;

	/// @brief The default value of this parameter.
	bool default_value_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; //class BooleanValuedParameter

} // namespace parametric
} // namespace conformation
} // namespace core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_conformation_parametric_BooleanValuedParameter )
#endif // SERIALIZATION


#endif
