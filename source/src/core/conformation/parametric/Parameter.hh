// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/conformation/parameters/Parameter.hh
/// @brief  Prototypes and method declarations for the Parameter class, a class for holding a single parameter for parametric backbone generation.
/// @author Vikram K. Mulligan (vmullig@uw.edu)


#ifndef INCLUDED_core_conformation_parametric_Parameter_hh
#define INCLUDED_core_conformation_parametric_Parameter_hh


// Unit headers
#include <core/conformation/parametric/Parameter.fwd.hh>
#include <core/conformation/parametric/Parameters.fwd.hh>

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
#include <utility/tag/XMLSchemaGeneration.hh>

// Numeric headers

// C++ headers


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace conformation {
namespace parametric {

/// @brief The types of parameters possible.
///
enum ParameterType {
	PT_generic_real = 1, //Keep this first.  Can be zero, negative, or positive.
	PT_generic_nonnegative_valued_real, //Can be zero or positive.
	PT_generic_positive_valued_real, //Must be greater than zero.
	PT_angle, // Will be stored in radians, in the range (-Pi, Pi], but might be entered in degrees by the user.
	PT_boolean, //True or false
	PT_generic_integer, // Negative or positive.
	PT_generic_whole_number, // 0 or positive.
	PT_generic_natural_number, // Strictly positive (nonzero).
	PT_generic_integer_vector, // Negative or positive.
	PT_generic_whole_number_vector, // 0 or positive.
	PT_generic_natural_number_vector, // Strictly positive (nonzero).
	PT_generic_real_vector,
	PT_generic_nonnegative_valued_real_vector,
	PT_generic_positive_valued_real_vector,
	PT_angle_vector,
	PT_invalid_type, //Keep this second-to-last.
	PT_end_of_list = PT_invalid_type //Keep this last.
};

/// @brief  Parameter class, used to store a single parameter for parametric backbone generation.
///
class Parameter : public utility::pointer::ReferenceCount, public utility::pointer::enable_shared_from_this< Parameter >
{
public:

	/// @brief constructors
	///
	Parameter();

	Parameter( Parameter const & src );

	/// @brief Pure virtual destructor.
	/// @details Counter-intuitively, C++ requires pure virtual destructors to be implemented.
	~Parameter() override=0;

	/// @brief Make a copy of this object ( allocate actual memory for it )
	///
	virtual
	ParameterOP clone() const=0;

	/// @brief Self pointers (const shared pointer).
	inline ParameterCOP get_self_ptr() const { return shared_from_this(); }

	/// @brief Self pointers (nonconst shared pointer).
	inline ParameterOP get_self_ptr() { return shared_from_this(); }

	/// @brief Self pointers (const weak pointer).
	inline ParameterCAP get_self_weak_ptr() const { return ParameterCAP( shared_from_this() ); }

	/// @brief Self pointers (nonconst weak pointer).
	inline ParameterAP get_self_weak_ptr() { return ParameterAP( shared_from_this() ); }

public: //Getters

	/// @brief Get parameter name.
	inline std::string const & parameter_name() const { return parameter_name_; }

	/// @brief Get parameter description.
	/// @details This is a lay-language description used for annotating user output.  It should be a short phrase with no capitals.
	inline std::string const & parameter_description() const { return parameter_description_; }

	/// @brief Get short parameter description.
	/// @details This is a one or two word lay-language description used for annotating user output.  The first word should be capitalized.
	inline std::string const & short_parameter_description() const { return short_parameter_description_; }

	/// @brief Get units.
	inline std::string const & parameter_units() const { return parameter_units_; }

	/// @brief Get the parameter type.
	inline ParameterType parameter_type() const { return parameter_type_; }

	/// @brief Is this parameter one that can be set directly by the user?
	/// @details Default true.  If false, it must be read from the database, e.g. from a Crick parameters file.
	inline bool can_be_set() const { return can_be_set_; }

	/// @brief Is this parameter one that can be copied from another parameter?
	/// @details Default true.
	inline bool can_be_copied() const { return can_be_copied_; }

	/// @brief Is this parameter one that can be sampled?
	/// @details Default true.  If false, it must be read from the database, e.g. from a Crick parameters file.
	inline bool can_be_sampled() const { return can_be_sampled_; }

	/// @brief Is this parameter one that can be perturbed?
	/// @details Default true.  If false, it must be read from the database, e.g. from a Crick parameters file.
	inline bool can_be_perturbed() const { return can_be_perturbed_; }

	/// @brief Get whether this property is meant to be global for a parameters set.
	inline bool global_for_parameters_set() const { return global_for_parameters_set_; }

	/// @brief Get the suffix for the copy tag (e.g. "copies_helix" in the case of helical bundles).
	/// @details Does not include the leading underscore.
	inline std::string const &copy_suffix() const { return copy_suffix_; }

	/// @brief Get the index of the Parameters object from which we will copy this parameter's value.
	/// @details A setting of zero means no copying.
	inline core::Size copy_from_parameters_index() const { return copy_from_parameters_index_; }

	/// @brief Has the value been set?
	inline bool value_was_set() const { return value_set_; }

	/// @brief Has the copying information been set?
	inline bool copying_information_was_set() const { return (copy_from_parameters_index_ != 0); }

public: //Setters

	/// @brief Set parameter name.
	void set_parameter_name( std::string const &name_in );

	/// @brief Set parameter description.
	/// @details This is a lay-language description used for annotating user output.  It should be a short phrase with no capitals.
	void set_parameter_description( std::string const & description_in );

	/// @brief Set short parameter description.
	/// @details This is a one or two word lay-language description used for annotating user output.  The first word should be capitalized.
	void set_short_parameter_description( std::string const & short_description_in );

	/// @brief Set parameter units.
	void set_parameter_units( std::string const & units_in );

	/// @brief Set the parameter type.
	virtual void set_parameter_type( ParameterType const type_in );

	/// @brief Set whether this parameter is one that can be set directly by the user, whether it can be sampled, and whether it can be perturbed.
	/// @details Default true in all cases.  If setting is false, the parameter value must be read from the database, e.g. from a Crick parameters file.
	void set_can_be_set_sampled_perturbed_copied( bool const can_be_set, bool const can_be_copied, bool const can_be_sampled, bool const can_be_perturbed );

	/// @brief Set whether this property is meant to be global for a parameters set.
	void set_global_for_parameters_set( bool const setting );

	/// @brief Set the suffix for the copy tag (e.g. "copies_helix" in the case of helical bundles).
	/// @details The leading underscore should be omitted.
	void set_copy_suffix( std::string const suffix_in );

	/// @brief Set the index of the Parameters object from which we will copy this parameter's value.
	/// @details A setting of zero means no copying.
	void set_copy_from_parameters_index( core::Size const setting );

	/// @brief Given another parameter of the same type, copy its value.  This does *not* set value_set_ to true.
	/// @details Must be implemented by derived classes.  Type checking can be performed by derived classes.
	/// @details Returns TRUE for failure and FALSE for success.
	virtual bool copy_value_from_parameter( ParameterCOP other_parameter, ParametersCOP other_parameter_collection, ParametersCOP this_parameter_collection ) = 0;

public: //Parse functions

	/// @brief Given a tag, parse out the setting for this parameter.
	/// @details If parse_setting is true, this tries to parse the value for this parameter.  If parse_grid_sampling is true, this tries to parse
	/// a range of values to sample, and a number of samples.  If parse_perturbation is true, this parses options for perturbing the value of the
	/// parameter.  If parse_copying is true, this parses options for copying.  Note that, if multiple options are true an error should be thrown
	/// if more than one of these is provided.
	/// @note Must be implemented by derived classes.
	virtual void parse_setting( utility::tag::TagCOP tag, bool const parse_setting, bool const parse_grid_sampling, bool const parse_perturbation, bool const parse_copying )=0;

	/// @brief Return the XSD information for this parameter.
	/// @note Must be implemented by derived classes.
	virtual void provide_xsd_information( utility::tag::AttributeList & xsd, bool const provide_setting, bool const provide_copying, bool const provide_grid_sampling, bool const provide_perturbation ) const = 0;

	/// @brief Reset the sampling and perturbation options before storing this Parameter object in a parametric Conformation.
	/// @details Pure virtual.  Must be implemented by derived classes.
	virtual void reset_sampling_and_perturbation_settings()=0;

	/// @brief Reset the copying options.
	virtual void reset_copying_settings();

	/// @brief Reset the value settings for this Parameter object.
	/// @details Pure virtual.  Must be implemented by derived classes.
	virtual void reset_value_settings()=0;

protected:

	/// @brief Allows derived classes to indicate that the value has been set.
	inline void set_value_was_set() { value_set_ = true; }

	/// @brief Allows derived classes to indicate that the value has not been set.
	inline void reset_value_was_set() { value_set_ = false; }

private:

	/********************************************************************************
	PRIVATE DATA
	*********************************************************************************/

	/// @brief The name of the parameter (e.g. "r0", "omega0", etc.)
	///
	std::string parameter_name_;

	/// @brief Parameter description.
	/// @details This is a lay-language description used for annotating user output.  It should be a short phrase with no capitals.
	std::string parameter_description_;

	/// @brief A shorter, lay-language descriptoin used for annotating user output.  Should be a word or two, with the first word captialized.
	std::string short_parameter_description_;

	/// @brief A word or two for the units.  Defaults to "dimensionless".
	std::string parameter_units_;

	/// @brief The type of data associated with the parameter (e.g. generic real, angle, integer, boolean).
	///
	ParameterType parameter_type_;

	/// @brief Is this parameter one that can be set directly by the user?
	/// @details Default true.  If false, it must be read from the database, e.g. from a Crick parameters file.
	bool can_be_set_;

	/// @brief Is this parameter one that can be copied from another parameter?
	/// @details Default true.
	bool can_be_copied_;

	/// @brief Is this parameter one that can be sampled?
	/// @details Default true.
	bool can_be_sampled_;

	/// @brief Is this parameter one that can be perturbed?
	/// @details Default true.
	bool can_be_perturbed_;

	/// @brief Is this a parameter that is intended to be global for a ParametersSet?
	/// @details Default false.
	bool global_for_parameters_set_;

	/// @brief A string appended to the parameter name for options related to copying the parameter value.
	/// @details For example, "copies_helix" in the case of parametric helices.
	std::string copy_suffix_;

	/// @brief Index of the Parameters object from which we should copy this parameter's value.
	core::Size copy_from_parameters_index_;

	/// @brief Was the value set?
	bool value_set_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; //class Parameter

} // namespace parametric
} // namespace conformation
} // namespace core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_conformation_parametric_Parameter )
#endif // SERIALIZATION

#endif
