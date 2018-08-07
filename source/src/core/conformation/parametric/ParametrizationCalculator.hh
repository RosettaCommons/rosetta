// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/conformation/parameters/ParametrizationCalculator.hh
/// @brief  Prototypes and method declarations for the ParametrizationCalculator class, a base class
/// from which classes that calculate particular parametrizations will be derived.
/// @author Vikram K. Mulligan (vmullig@uw.edu)


#ifndef INCLUDED_core_conformation_parametric_ParametrizationCalculator_hh
#define INCLUDED_core_conformation_parametric_ParametrizationCalculator_hh


// Unit headers
#include <core/conformation/parametric/ParametrizationCalculator.fwd.hh>
#include <core/conformation/parametric/Parameters.hh>
#include <core/conformation/parametric/Parameter.hh>
#include <core/conformation/parametric/RealValuedParameter.fwd.hh>
#include <core/conformation/parametric/RealVectorValuedParameter.fwd.hh>
#include <core/conformation/parametric/SizeValuedParameter.fwd.hh>
#include <core/conformation/parametric/SizeVectorValuedParameter.fwd.hh>
#include <core/conformation/parametric/BooleanValuedParameter.fwd.hh>

// Package headers

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// Numeric headers

// C++ headers


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace conformation {
namespace parametric {

struct ParameterizationCalculatorProperties {

	bool can_be_set;
	bool can_be_copied;
	bool can_be_sampled;
	bool can_be_perturbed;
	bool global_for_parameters_set;

	/// @brief Struct default constructor.
	ParameterizationCalculatorProperties() :
		can_be_set(true),
		can_be_copied(true),
		can_be_sampled(true),
		can_be_perturbed(true),
		global_for_parameters_set(true)
	{}

	/// @brief Struct options constructor.
	ParameterizationCalculatorProperties( bool const setting_allowed, bool const copying_allowed, bool const sampling_allowed, bool const perturbing_allowed, bool const is_global ) :
		can_be_set(setting_allowed),
		can_be_copied(copying_allowed),
		can_be_sampled(sampling_allowed),
		can_be_perturbed(perturbing_allowed),
		global_for_parameters_set(is_global)
	{}

};

/// @brief  ParametrizationCalculator class, used for parametric backbone generation.
///
class ParametrizationCalculator : public utility::pointer::ReferenceCount, public utility::pointer::enable_shared_from_this< ParametrizationCalculator >
{
public:

	/// @brief Constructor.
	///
	ParametrizationCalculator();

	/// @brief Constructor with input Parameters object.
	///
	ParametrizationCalculator( ParametersOP parameters );

	/// @brief Copy constructor.
	/// @details Deep-copies the stored parameters.
	ParametrizationCalculator( ParametrizationCalculator const & src );

	/// @brief Pure virtual destructor means don't instantiate this base class!.
	/// @details Counter-intuitively, C++ requires pure virtual destructors to be implemented.
	~ParametrizationCalculator() override=0;

	/// @brief Copy this object ( allocate actual memory for it )
	/// @details Must be implemented by derived classes.
	virtual
	ParametrizationCalculatorOP clone() const=0;

	/// self pointers
	inline ParametrizationCalculatorCOP get_self_ptr() const { return shared_from_this(); }
	inline ParametrizationCalculatorOP get_self_ptr() { return shared_from_this(); }
	inline ParametrizationCalculatorCAP get_self_weak_ptr() const { return ParametrizationCalculatorCAP( shared_from_this() ); }
	inline ParametrizationCalculatorAP get_self_weak_ptr() { return ParametrizationCalculatorAP( shared_from_this() ); }

public:

	/// @brief Get a const-owning pointer to the Parameters object.
	inline ParametersCOP parameters_cop() const { return ParametersCOP( parameters_ ); }

	/// @brief Access a parameter.
	///
	ParameterOP parameter( core::Size const index );

	/// @brief Access a parameter.
	/// @details Const-access
	ParameterCOP parameter_cop( core::Size const index ) const;

	/// @brief Accesses a real-valued parameter.
	/// @details Returns nullptr if the parameter with the given index is not a RealValuedParameter.
	RealValuedParameterOP real_parameter( core::Size const index );

	/// @brief Accesses a realvector-valued parameter.
	/// @details Returns nullptr if the parameter with the given index is not a RealVectorValuedParameter.
	RealVectorValuedParameterOP realvector_parameter( core::Size const index );

	/// @brief Accesses a size-valued parameter.
	/// @details Returns nullptr if the parameter with the given index is not a SizeValuedParameter.
	SizeValuedParameterOP size_parameter( core::Size const index );

	/// @brief Accesses a sizevector-valued parameter.
	/// @details Returns nullptr if the parameter with the given index is not a SizeVectorValuedParameter.
	SizeVectorValuedParameterOP sizevector_parameter( core::Size const index );

	/// @brief Accesses a boolean-valued parameter.
	/// @details Returns nullptr if the parameter with the given index is not a BooleanValuedParameter.
	BooleanValuedParameterOP boolean_parameter( core::Size const index );

	/// @brief Accesses a real-valued parameter (const-access).
	/// @details Returns nullptr if the parameter with the given index is not a RealValuedParameter.
	RealValuedParameterCOP real_parameter_cop( core::Size const index ) const;

	/// @brief Accesses a realvector-valued parameter (const-access).
	/// @details Returns nullptr if the parameter with the given index is not a RealVectorValuedParameter.
	RealVectorValuedParameterCOP realvector_parameter_cop( core::Size const index ) const;

	/// @brief Accesses a size-valued parameter (const-access).
	/// @details Returns nullptr if the parameter with the given index is not a SizeValuedParameter.
	SizeValuedParameterCOP size_parameter_cop( core::Size const index ) const;

	/// @brief Accesses a sizevector-valued parameter (const-access).
	/// @details Returns nullptr if the parameter with the given index is not a SizeVectorValuedParameter.
	SizeVectorValuedParameterCOP sizevector_parameter_cop( core::Size const index ) const;

	/// @brief Accesses a boolean-valued parameter (const-access).
	/// @details Returns nullptr if the parameter with the given index is not a BooleanValuedParameter.
	BooleanValuedParameterCOP boolean_parameter_cop( core::Size const index ) const;

protected: //Functions

	/// @brief Add a parameter to this calculator, automatically determining the type.
	///
	void add_parameter( std::string const &parameter_name, ParameterType type, std::string const &description, std::string const &short_description, std::string const &units, ParameterizationCalculatorProperties const &properties );

	/// @brief Add a real-valued parameter to this calculator.
	///
	void add_real_parameter( std::string const &parameter_name, ParameterType type, std::string const &description, std::string const &short_description, std::string const &units, ParameterizationCalculatorProperties const &properties );

	/// @brief Add a vector-valued parameter to this calculator.
	void add_realvector_parameter( std::string const &parameter_name, ParameterType type, std::string const &description, std::string const &short_description, std::string const &units, ParameterizationCalculatorProperties const &properties );

	/// @brief Add a core::Size-valued parameter to this calculator.
	///
	void add_size_parameter( std::string const &parameter_name, ParameterType type, std::string const &description, std::string const &short_description, std::string const &units, ParameterizationCalculatorProperties const &properties );

	/// @brief Add a integer vector-valued parameter to this calculator.
	///
	void add_sizevector_parameter( std::string const &parameter_name, ParameterType type, std::string const &description, std::string const &short_description, std::string const &units, ParameterizationCalculatorProperties const &properties );

	/// @brief Add a Boolean-valued parameter to this calculator.
	///
	void add_boolean_parameter( std::string const &parameter_name, ParameterType type, std::string const &description, std::string const &short_description, std::string const &units, ParameterizationCalculatorProperties const &properties );

	/// @brief Add a custom parameter
	void add_custom_parameter( std::string const &parameter_name, ParameterType const type, std::string const &description, std::string const &short_description, std::string const &units, ParameterizationCalculatorProperties const &properties, ParameterOP parameter_in );

private: //Functions

	/// @brief Add a parameter to this calculator.
	///
	void add_parameter( ParameterOP parameter_in );

private:

	/********************************************************************************
	PRIVATE DATA
	*********************************************************************************/

	/// @brief The container of parameters that the user can twiddle given a particular
	/// parametrization.
	ParametersOP parameters_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; //class ParametrizationCalculator

} // namespace parametric
} // namespace conformation
} // namespace core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_conformation_parametric_ParametrizationCalculator )
#endif // SERIALIZATION


#endif
