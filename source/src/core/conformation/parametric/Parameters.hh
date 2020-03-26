// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/conformation/parameters/Parameters.hh
/// @brief  Prototypes and method declarations for the Parameters class, a class for holding parameters for parametric backbone generation.
/// @author Vikram K. Mulligan (vmullig@uw.edu)


#ifndef INCLUDED_core_conformation_parametric_Parameters_hh
#define INCLUDED_core_conformation_parametric_Parameters_hh


// Unit headers
#include <core/conformation/parametric/Parameters.fwd.hh>
#include <core/conformation/parametric/Parameter.fwd.hh>

// Package headers
#include <core/conformation/Residue.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>
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

/// @brief  Parameters class, used to store sets of parameters for parametric backbone generation.
///
class Parameters : public utility::VirtualBase, public utility::pointer::enable_shared_from_this< Parameters >
{
public:

	/// @brief constructors
	///
	Parameters();

	/// @brief Copy constructor.
	/// @details Deep-copies the residue list and the parameters list.
	Parameters( Parameters const & src );

	~Parameters() override;

	/// @brief Copy this residue( allocate actual memory for it )
	///
	virtual
	ParametersOP clone() const;

	/// self pointers
	inline ParametersCOP get_self_ptr() const { return shared_from_this(); }
	inline ParametersOP get_self_ptr() { return shared_from_this(); }
	inline ParametersCAP get_self_weak_ptr() const { return ParametersCAP( shared_from_this() ); }
	inline ParametersAP get_self_weak_ptr() { return ParametersAP( shared_from_this() ); }

	/// @brief Add a residue to the list of residues that these parameters describe.
	void add_residue( core::conformation::ResidueCOP residue );

	/// @brief Get an owning pointer to a residue in the list of residues that these parameters describe.
	core::conformation::ResidueCOP residue( core::Size const index ) const;

	/// @brief Get a const owning pointer to the first residue in the list of residues that these parameters describe.
	/// @details  Note that this might not be the first residue in linear sequence, if the residues were put in in non-
	/// sequential order or the residue numbering has changed.
	core::conformation::ResidueCOP first_residue() const;

	/// @brief Get a const owning pointer to the last residue in the list of residues that these parameters describe.
	/// @details  Note that this might not be the last residue in linear sequence, if the residues were put in in non-
	/// sequential order or the residue numbering has changed.
	core::conformation::ResidueCOP last_residue() const;

	/// @brief Get the index of the first residue.
	/// @details MUST BE REWRITTEN when the OP issue is resolved.
	core::Size first_residue_index() const;

	/// @brief Get the index of the last residue.
	/// @details MUST BE REWRITTEN when the OP issue is resolved.
	core::Size last_residue_index() const;

	/// @brief Returns the number (count) of residues that these parameters describe.
	core::Size n_residue() const;

	/// @brief Clears the list of residues that these parameters describe.
	///
	void reset_residue_list();

	/// @brief Clears the sampling and perturbing information in the individual parameters.
	void reset_sampling_and_perturbing_info() const;

	/// @brief Assign an element in the residue list to be an owning pointer to an existing residue.
	void set_residue( core::Size const index, core::conformation::ResidueCOP existing_residue );

	/// @brief Get the number of parameters stored in this Parameters object.
	core::Size num_parameters() const;

	/// @brief Add a parameter.
	/// @details Does NOT clone the parameter, but stores the OP directly.
	void add_parameter( ParameterOP parameter );

	/// @brief Access a parameter, by index.
	/// @details Non-const access.
	ParameterOP parameter_op( core::Size const index ) const;

	/// @brief Access a parameter, by index.
	/// @details Const access.
	ParameterCOP parameter_cop( core::Size const index ) const;

	/// @brief Replace one of the contained parameter objects with
	/// a copy of an input parameter object.
	/// @details If reset_sampling_copying_perturbing is true, then
	/// the sampling, copying, and perturbing information of the
	/// copied parameter object are all reset.  True by default.
	void replace_parameter_via_clone( core::Size const param_index, ParameterCOP new_parameter, bool const reset_sampling_copying_perturbing=true );

private:

	/********************************************************************************
	PRIVATE DATA
	*********************************************************************************/

	/// @brief A list of the residues in this conformation that are described by the parameters in this Parameters object.
	/// @details This is a list of owning pointers so that the residue indices don't mess things up.  (That is, as residue
	/// indices change, these should still point to the proper residues).
	utility::vector1 < core::conformation::ResidueCOP > residue_list_;

	/// @brief A list of the parameters that describe these residues.
	utility::vector1< ParameterOP > parameter_list_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; //class Parameters

} // namespace parametric
} // namespace conformation
} // namespace core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_conformation_parametric_Parameters )
#endif // SERIALIZATION

#endif
