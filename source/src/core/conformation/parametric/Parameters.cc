// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/conformation/parametric/Parameters.cc
/// @brief  A class for holding sets of parameters for parametric backbone generation.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Unit header
#include <core/conformation/parametric/Parameters.hh>
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

static basic::Tracer TR( "core.conformation.parametric.Parameters" );

/// @brief Constructor.
///
Parameters::Parameters() :
	residue_list_(),
	parameter_list_()
{
}

/// @brief Copy constructor.
/// @details Deep-copies the residue list and the parameters list.
Parameters::Parameters( Parameters const & src ) :
	utility::pointer::ReferenceCount(),
	utility::pointer::enable_shared_from_this< Parameters >()
{
	residue_list_.clear();
	if ( src.residue_list_.size()>0 ) {
		residue_list_.reserve(src.residue_list_.size());
		for ( core::Size i(1), imax(src.residue_list_.size()); i<=imax; ++i ) {
			residue_list_.push_back( src.residue_list_[i]->clone() ); //This copies the residue that was being pointed at.
			//Note that when copying a Conformation, I need to add logic that will ensure that the Parameters objects that result have owning pointers to the residues in the Conformation,
			//rather than to residues that only exist in the Parameters object.
		}
	}

	parameter_list_.clear();
	if ( src.parameter_list_.size() > 0 ) {
		parameter_list_.reserve(src.parameter_list_.size());
		for ( core::Size i(1), imax(src.parameter_list_.size()); i<=imax; ++i ) {
			parameter_list_.push_back( src.parameter_list_[i]->clone() );
		}
	}
}

Parameters::~Parameters() = default;


/// @brief make a copy of this residue( allocate actual memory for it )
///
ParametersOP
Parameters::clone() const
{
	return ParametersOP( new Parameters( *this ) );
}

/// @brief Clears the sampling and perturbing information in the individual parameters.
void
Parameters::reset_sampling_and_perturbing_info() const {
	for ( core::Size i(1), imax(parameter_list_.size()); i<=imax; ++i ) {
		parameter_list_[i]->reset_sampling_and_perturbation_settings();
	}
}

/// @brief Add a parameter.
/// @details Does NOT clone the parameter, but stores the OP directly.
void
Parameters::add_parameter(
	ParameterOP parameter
) {
	debug_assert( parameter != nullptr );
	parameter_list_.push_back( parameter );
}

/// @brief Access a parameter, by index.
/// @details Non-const access.
ParameterOP
Parameters::parameter_op(
	core::Size const index
) const {
	debug_assert( index > 0 );
	debug_assert( index <= parameter_list_.size() );
	return parameter_list_[index];
}

/// @brief Access a parameter, by index.
/// @details Const access.
ParameterCOP
Parameters::parameter_cop(
	core::Size const index
) const {
	debug_assert( index > 0 );
	debug_assert( index <= parameter_list_.size() );
	return ParameterCOP(parameter_list_[index]);
}

/// @brief Replace one of the contained parameter objects with
/// a copy of an input parameter object.
void
Parameters::replace_parameter_via_clone(
	core::Size const param_index,
	ParameterCOP new_parameter,
	bool const reset_sampling_copying_perturbing/*=true*/
) {
	debug_assert( param_index > 0 );
	debug_assert( param_index <= parameter_list_.size() );
	parameter_list_[param_index] = new_parameter->clone();
	if ( reset_sampling_copying_perturbing ) {
		parameter_list_[param_index]->reset_copying_settings();
		parameter_list_[param_index]->reset_sampling_and_perturbation_settings();
	}
}

} // namespace parametric
} // namespace conformation
} // namespace core


#ifdef SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::conformation::parametric::Parameters::save( Archive & arc ) const {
	arc( CEREAL_NVP( residue_list_ ) ); // utility::vector1<core::conformation::ResidueCOP>
	arc( CEREAL_NVP( parameter_list_ ) ); // utility::vector1<core::conformation::ResidueCOP>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::conformation::parametric::Parameters::load( Archive & arc ) {
	utility::vector1< std::shared_ptr< core::conformation::Residue > > local_residue_list;
	arc( local_residue_list ); // utility::vector1<core::conformation::ResidueCOP>
	residue_list_ = local_residue_list; // copy the non-const pointer(s) into the const pointer(s)
	arc( parameter_list_ );
}

SAVE_AND_LOAD_SERIALIZABLE( core::conformation::parametric::Parameters );
CEREAL_REGISTER_TYPE( core::conformation::parametric::Parameters )

CEREAL_REGISTER_DYNAMIC_INIT( core_conformation_parametric_Parameters )
#endif // SERIALIZATION
